"""Read one structure with one HTTP range request.

Each member of a bulk tar is a complete gzip file, so the byte range
``[offset, offset + length)`` of the tar is the ``.cif.gz`` or ``.pdb.gz`` itself.

What a read checks, in order:

1. The answer must be ``206 Partial Content`` for exactly the range asked for. A ``200``
   means the server ignored ``Range`` and is about to send the whole tar (up to 40 GB):
   the connection is closed without reading the body and :class:`FetchError` is raised.
2. The tar's total size, ``ETag`` and ``Last-Modified`` must equal what the index
   recorded when it was built. A difference means the tar was replaced and the offsets
   can no longer be trusted: :class:`StaleIndexError`, raised before the body is read.
3. The bytes must gunzip. If they do not: :class:`StaleIndexError`.

Retries: up to four attempts with exponential backoff, on timeouts, dropped connections
and 5xx answers only. A 4xx answer is never retried.
"""

from __future__ import annotations

import gzip
import http.client
import os
import random
import re
import tempfile
import time
import urllib.error
import urllib.request
import zlib
from collections.abc import Iterable
from concurrent.futures import ThreadPoolExecutor
from email.message import Message
from pathlib import Path

from . import _index
from ._errors import FetchError, FetchManyError, StaleIndexError
from ._index import Index, Member, Tar
from ._version import USER_AGENT

DEFAULT_TIMEOUT = 60.0
DEFAULT_WORKERS = 4

_MAX_ATTEMPTS = 4
_BACKOFF_BASE = 0.5
_BACKOFF_MAX = 8.0
_sleep = time.sleep

_CONTENT_RANGE_RE = re.compile(r"^bytes\s+(\d+)-(\d+)/(\d+|\*)$")
_REFRESH = "The index is out of date: call scigantic_alphafold.refresh_index() (or run `scigantic-alphafold refresh`)."


class _Retryable(Exception):
    """A failure worth another attempt: a timeout, a dropped connection or a 5xx answer."""


def _same_etag(indexed: str, served: str) -> bool:
    def strong(value: str) -> str:
        value = value.strip()
        return value[2:] if value.startswith("W/") else value

    return strong(indexed) == strong(served)


def _check_validators(member: Member, tar: Tar, headers: Message, total: str) -> None:
    served_etag = headers.get("ETag")
    served_modified = headers.get("Last-Modified")
    problem: str | None = None
    if tar.etag and served_etag and not _same_etag(tar.etag, served_etag):
        problem = f"its ETag is now {served_etag.strip()}, the index recorded {tar.etag}"
    elif tar.last_modified and served_modified and tar.last_modified.strip() != served_modified.strip():
        problem = f"its Last-Modified is now {served_modified.strip()}, the index recorded {tar.last_modified}"
    elif total != "*" and int(total) != tar.size_bytes:
        problem = f"its size is now {total} bytes, the index recorded {tar.size_bytes}"
    if problem is not None:
        raise StaleIndexError(f"{tar.name} changed on the server since the index was built: {problem}. {_REFRESH}")


def _request_once(member: Member, tar: Tar, timeout: float) -> bytes:
    first = member.offset
    last = member.offset + member.length - 1
    request = urllib.request.Request(
        member.url,
        headers={"Range": f"bytes={first}-{last}", "User-Agent": USER_AGENT, "Accept-Encoding": "identity"},
    )
    try:
        response = urllib.request.urlopen(request, timeout=timeout)
    except urllib.error.HTTPError as exc:
        code = exc.code
        exc.close()
        if code >= 500:
            raise _Retryable(f"HTTP {code}") from None
        if code == 416:
            raise StaleIndexError(
                f"{tar.name} is shorter on the server than the index expects (HTTP 416 for bytes "
                f"{first}-{last}). {_REFRESH}"
            ) from None
        raise FetchError(f"HTTP {code} for {member.url} (not retried)", status=code) from None
    except (OSError, http.client.HTTPException) as exc:
        raise _Retryable(str(exc) or type(exc).__name__) from None

    # From here the connection is open. Closing it without reading is what keeps a
    # full-tar 200 answer from being downloaded.
    try:
        status = response.status
        if status != 206:
            raise FetchError(
                f"{member.url} answered HTTP {status} instead of 206 to a range request, so it would have "
                f"sent the whole tar ({tar.size_bytes} bytes). The connection was closed without reading it.",
                status=status,
            )
        content_range = _CONTENT_RANGE_RE.match((response.headers.get("Content-Range") or "").strip())
        if content_range is None:
            raise FetchError(f"{member.url} answered 206 without a usable Content-Range header", status=status)
        _check_validators(member, tar, response.headers, content_range.group(3))
        if int(content_range.group(1)) != first or int(content_range.group(2)) != last:
            raise FetchError(
                f"{member.url} answered with bytes {content_range.group(1)}-{content_range.group(2)}, "
                f"not the requested {first}-{last}",
                status=status,
            )
        declared = response.headers.get("Content-Length")
        if declared is not None and declared.strip() != str(member.length):
            raise FetchError(
                f"{member.url} answered with Content-Length {declared.strip()}, expected {member.length}",
                status=status,
            )
        try:
            data: bytes = response.read(member.length)
        except (OSError, http.client.HTTPException) as exc:
            raise _Retryable(str(exc) or type(exc).__name__) from None
        if len(data) != member.length:
            raise _Retryable(f"connection closed after {len(data)} of {member.length} bytes")
        return data
    finally:
        response.close()


def _gunzip(member: Member, data: bytes) -> bytes:
    try:
        return gzip.decompress(data)
    except (OSError, EOFError, zlib.error) as exc:
        raise StaleIndexError(
            f"The bytes at offset {member.offset} of {member.tar} are not the gzip file of "
            f"{member.entry_id} ({exc}). {_REFRESH}"
        ) from None


def read_member(member: Member, tar: Tar, *, timeout: float = DEFAULT_TIMEOUT) -> tuple[bytes, bytes]:
    """Range-read one member. Returns (the stored gzip bytes, the decompressed bytes)."""
    failure = ""
    for attempt in range(_MAX_ATTEMPTS):
        if attempt:
            _sleep(min(_BACKOFF_BASE * 2 ** (attempt - 1), _BACKOFF_MAX) * (1.0 + random.random() / 4))
        try:
            data = _request_once(member, tar, timeout)
        except _Retryable as exc:
            failure = str(exc)
            continue
        return data, _gunzip(member, data)
    raise FetchError(f"Gave up on {member.url} after {_MAX_ATTEMPTS} attempts: {failure}")


def _read(index: Index, id: str, fmt: str, timeout: float) -> tuple[Member, bytes, bytes]:
    member = index.resolve(id, fmt)
    raw, text = read_member(member, index.tar(member.tar), timeout=timeout)
    return member, raw, text


def fetch(id: str, fmt: str = "cif", *, timeout: float = DEFAULT_TIMEOUT) -> str:
    """The structure file as text (mmCIF by default, ``fmt="pdb"`` for PDB).

    ``id`` is a UniProt accession (``P69905``), an isoform accession (``P12345-2``) or an
    entry id (``AF-P69905-F1``), in any letter case. A bare accession means its single
    ``F1`` entry; :class:`AmbiguousEntryError` lists the candidates when there are several
    fragments. :class:`NotInBulkSetError` when the accession is not in the bulk tars.
    """
    return _read(_index.get_index(), id, fmt, timeout)[2].decode("utf-8")


def fetch_bytes(id: str, fmt: str = "cif", *, timeout: float = DEFAULT_TIMEOUT) -> bytes:
    """The ``.cif.gz`` or ``.pdb.gz`` exactly as stored in the tar (checked to gunzip)."""
    return _read(_index.get_index(), id, fmt, timeout)[1]


def download(
    id: str,
    out_dir: str | os.PathLike[str] = ".",
    fmt: str = "cif",
    decompress: bool = True,
    *,
    timeout: float = DEFAULT_TIMEOUT,
) -> Path:
    """Write ``AF-<accession>-F1-model_v6.cif`` (``.cif.gz`` with ``decompress=False``) and return its path.

    The file is written under a temporary name in ``out_dir`` and renamed into place, so
    a reader never sees a partial file. ``out_dir`` is created if it does not exist.
    """
    index = _index.get_index()
    member, raw, text = _read(index, id, fmt, timeout)
    version = index.meta().get("afdb_version", "6")
    name = f"{member.entry_id}-model_v{version}.{member.fmt}" + ("" if decompress else ".gz")
    folder = Path(out_dir)
    folder.mkdir(parents=True, exist_ok=True)
    target = folder / name
    fd, temp = tempfile.mkstemp(dir=folder, prefix=f".{name}.", suffix=".part")
    try:
        with os.fdopen(fd, "wb") as out:
            out.write(text if decompress else raw)
        os.chmod(temp, 0o644)
        os.replace(temp, target)
    except BaseException:
        try:
            os.unlink(temp)
        except FileNotFoundError:
            pass
        raise
    return target


def fetch_many(
    ids: Iterable[str],
    fmt: str = "cif",
    workers: int = DEFAULT_WORKERS,
    *,
    on_error: str = "raise",
    timeout: float = DEFAULT_TIMEOUT,
) -> dict[str, str]:
    """Fetch several structures in parallel. Returns ``{id as given: text}``.

    Every id is attempted, whatever happens to the others. When all succeed the dict is
    complete, in the order given (duplicates collapse to one key). When some fail:

    * ``on_error="raise"`` (the default) raises :class:`FetchManyError` once everything
      has finished. Its ``.results`` holds the successes, exactly as the return value
      would have, and its ``.errors`` maps each failed id to its exception.
    * ``on_error="skip"`` returns only the successes.

    ``workers`` is the number of simultaneous connections (default 4).
    """
    if on_error not in ("raise", "skip"):
        raise ValueError(f"on_error must be 'raise' or 'skip', got {on_error!r}")
    _index.check_fmt(fmt)
    wanted = list(dict.fromkeys([ids] if isinstance(ids, str) else ids))
    if not wanted:
        return {}
    index = _index.get_index()

    def one(id: str) -> str | Exception:
        try:
            return _read(index, id, fmt, timeout)[2].decode("utf-8")
        except Exception as exc:
            return exc

    with ThreadPoolExecutor(max_workers=max(1, min(int(workers), len(wanted)))) as pool:
        outcomes = list(pool.map(one, wanted))
    results: dict[str, str] = {}
    errors: dict[str, Exception] = {}
    for id, outcome in zip(wanted, outcomes):
        if isinstance(outcome, Exception):
            errors[id] = outcome
        else:
            results[id] = outcome
    if errors and on_error == "raise":
        raise FetchManyError(results, errors)
    return results
