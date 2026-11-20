"""Build the byte-offset index for the AlphaFold DB v6 bulk tars.

Nothing in this module is needed to read structures. It produces the index
that the reader uses:

* :func:`build_part` streams one tar over HTTP (never to disk), walks it with
  ``tarfile`` in stream mode and writes ``<tar name>.tsv.gz`` with the offset
  and length of every member's data bytes.
* :func:`merge_parts` turns a directory of part files into ``index-v6.sqlite``
  and ``index-v6.sqlite.gz`` and refuses to do so if anything is missing or a
  count does not match EBI's ``download_metadata.json``.
* :func:`verify_index` reads random members back from the live server with the
  recorded offsets and checks that each one is the model it claims to be.
* :func:`write_summary` and :func:`needs_rebuild` support the scheduled
  rebuild in ``.github/workflows/build-index.yml``.

Run ``python -m scigantic_alphafold._build --help`` for the command line used
by the workflow. Standard library only.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import gzip
import hashlib
import http.client
import json
import os
import random
import re
import shutil
import sqlite3
import sys
import tarfile
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import IO, Any, Iterable, NamedTuple, Sequence, cast

from ._version import USER_AGENT

__all__ = [
    "V6_BASE_URL",
    "METADATA_URL",
    "SCHEMA_VERSION",
    "AFDB_VERSION",
    "BuildError",
    "TarChangedError",
    "MergeError",
    "VerifyError",
    "build_part",
    "merge_parts",
    "verify_index",
    "write_summary",
    "needs_rebuild",
    "plan_jobs",
    "segment_count",
    "segment_bounds",
    "part_file_name",
    "SEGMENT_BYTES",
    "SEGMENT_THRESHOLD_BYTES",
    "fetch_download_metadata",
    "main",
]

V6_BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/alphafold/v6"
METADATA_URL = "https://ftp.ebi.ac.uk/pub/databases/alphafold/download_metadata.json"
SCHEMA_VERSION = "1"
AFDB_VERSION = "6"

#: Member names that are indexed. Anything else lands in a ``#skipped`` line.
#: Verified 2026-11-20 against proteome, global health and both Swiss-Prot
#: tars: all use ``AF-<accession>-F<n>-model_v6.<cif|pdb>.gz`` at the tar root.
_MEMBER_RE = re.compile(r"(AF-[A-Z0-9]+(?:-[0-9]+)?-F[0-9]+)-model_v6\.(cif|pdb)\.gz")
_ENTRY_RE = re.compile(r"AF-([A-Z0-9]+)(?:-[0-9]+)?-F[0-9]+")
_PROTEOME_TAR_RE = re.compile(r"(UP[0-9]+)_([0-9]+)_[A-Za-z0-9]+_v6\.tar")
_SWISSPROT_TAR_RE = re.compile(r"swissprot_(cif|pdb)_v6\.tar")
_TAR_NAME_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]*\.tar")

_PART_SUFFIX = ".tsv.gz"
_FORMATS = ("cif", "pdb")

#: Block size handed to ``tarfile`` for the streaming walk.
_BUFSIZE = 64 * 1024
#: Upper bound on reconnects for one tar, whatever progress is being made.
_MAX_TOTAL_RECONNECTS = 500
_MAX_BACKOFF_SECONDS = 60.0

#: Tars larger than the threshold are built in parallel byte-range segments of
#: at most SEGMENT_BYTES each. A hosted runner measured about 1.2 MB/s per
#: stream from EBI (2026-11-20), so 4 GiB is roughly an hour per job.
SEGMENT_BYTES = 4 * 1024**3
SEGMENT_THRESHOLD_BYTES = 6 * 1024**3
#: A segment looks at most this far past its end for the next member header,
#: and does not accept a resynchronization candidate larger than this.
_RESYNC_LIMIT = 64 * 1024 * 1024

#: Tars whose member count is known to differ from ``num_predicted_structures``
#: in download_metadata.json, as ``{tar name: {fmt: expected member count}}``.
#: Empty on purpose: every tar measured so far has exactly one member per
#: format per predicted structure. Add an entry here, with a comment saying
#: what was measured and when, rather than loosening the check in merge_parts.
KNOWN_COUNT_EXCEPTIONS: dict[str, dict[str, int]] = {}

_SCHEMA = """
CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT NOT NULL);
CREATE TABLE tars(
  tar_id INTEGER PRIMARY KEY,
  name TEXT NOT NULL UNIQUE,
  url TEXT NOT NULL,
  size_bytes INTEGER NOT NULL,
  last_modified TEXT,
  etag TEXT,
  type TEXT,
  reference_proteome TEXT,
  tax_id INTEGER,
  species TEXT,
  common_name TEXT,
  n_members INTEGER NOT NULL
);
CREATE TABLE members(
  entry_id TEXT NOT NULL,
  fmt TEXT NOT NULL,
  tar_id INTEGER NOT NULL,
  offset INTEGER NOT NULL,
  length INTEGER NOT NULL,
  PRIMARY KEY(entry_id, fmt, tar_id)
) WITHOUT ROWID;
"""


class BuildError(RuntimeError):
    """A tar could not be streamed or parsed into a part file."""


class TarChangedError(BuildError):
    """The tar changed on the server between two requests of one build."""


class MergeError(RuntimeError):
    """The part files do not add up to a complete, consistent index."""


class VerifyError(RuntimeError):
    """A spot check read from the live server did not return a valid model."""


def _log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def _sleep(seconds: float) -> None:
    time.sleep(seconds)


def _is_retryable(exc: BaseException) -> bool:
    """Server errors, timeouts and dropped connections. Never 4xx."""
    if isinstance(exc, urllib.error.HTTPError):
        return exc.code >= 500
    return isinstance(exc, (OSError, http.client.HTTPException))


def _request(
    url: str, *, headers: dict[str, str] | None = None, method: str = "GET"
) -> urllib.request.Request:
    merged = {"User-Agent": USER_AGENT, "Accept-Encoding": "identity"}
    if headers:
        merged.update(headers)
    return urllib.request.Request(url, headers=merged, method=method)


def _open(request: urllib.request.Request, timeout: float) -> http.client.HTTPResponse:
    return cast(http.client.HTTPResponse, urllib.request.urlopen(request, timeout=timeout))


def _with_retries(what: str, attempt: Any, *, max_retries: int, backoff: float) -> Any:
    failures = 0
    while True:
        try:
            return attempt()
        except Exception as exc:
            if not _is_retryable(exc):
                raise
            failures += 1
            if failures > max_retries:
                raise
            delay = min(backoff * 2 ** (failures - 1), _MAX_BACKOFF_SECONDS)
            _log(f"{what}: {exc!r}; retry {failures}/{max_retries} in {delay:.1f}s")
            _sleep(delay)


# --------------------------------------------------------------------------
# Streaming one tar


class _HTTPTarStream:
    """Sequential, read-only view of one remote file that survives disconnects.

    ``read`` hands out the file's bytes in order. When the connection drops,
    times out or the server answers 5xx, the stream reconnects with
    ``Range: bytes=<pos>-`` and carries on from the exact byte it had reached,
    so the tar parser above it never notices. Every reconnect checks that the
    ETag (or Last-Modified when there is no ETag) and the total size are
    unchanged; if they are not, the offsets recorded so far describe a file
    that no longer exists and the build stops with :class:`TarChangedError`.
    """

    def __init__(
        self,
        url: str,
        *,
        start: int = 0,
        timeout: float = 60.0,
        max_retries: int = 8,
        backoff: float = 2.0,
    ) -> None:
        self.url = url
        self.timeout = timeout
        self.max_retries = max_retries
        self.backoff = backoff
        #: First byte this stream serves; a segment starts in the middle.
        self.start = start
        self.pos = start
        self.size = -1
        self.etag: str | None = None
        self.last_modified: str | None = None
        self.reconnects = 0
        #: End offset of the last non-zero byte handed out so far.
        self.last_nonzero_end = 0
        self._resp: http.client.HTTPResponse | None = None
        self._failures = 0
        self._retrying(self._connect)

    # -- connection handling

    def _connect(self) -> None:
        headers: dict[str, str] = {}
        if self.pos > 0:
            headers["Range"] = f"bytes={self.pos}-"
            if self.size >= 0 and self.etag and not self.etag.startswith("W/"):
                headers["If-Range"] = self.etag
        resp = _open(_request(self.url, headers=headers), self.timeout)
        try:
            self._check(resp)
        except BaseException:
            resp.close()
            raise
        self._resp = resp

    def _check(self, resp: http.client.HTTPResponse) -> None:
        etag = resp.headers.get("ETag")
        last_modified = resp.headers.get("Last-Modified")
        first = self.size < 0
        if not first and (etag != self.etag or last_modified != self.last_modified):
            raise TarChangedError(
                f"{self.url} changed during the build: ETag {self.etag!r} -> {etag!r}, "
                f"Last-Modified {self.last_modified!r} -> {last_modified!r}"
            )
        problem: type[BuildError] = BuildError if first else TarChangedError
        if self.pos == 0:
            if resp.status != 200:
                raise problem(f"{self.url}: expected HTTP 200, got {resp.status}")
            length = resp.headers.get("Content-Length")
            if length is None:
                raise BuildError(f"{self.url}: response has no Content-Length")
            total = int(length)
        else:
            if resp.status != 206:
                raise problem(
                    f"{self.url}: asked for bytes from {self.pos} and got HTTP {resp.status} "
                    "instead of 206; the server does not serve this file by range"
                )
            content_range = resp.headers.get("Content-Range", "").strip()
            match = re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)", content_range)
            if (
                match is None
                or int(match.group(1)) != self.pos
                or int(match.group(2)) != int(match.group(3)) - 1
            ):
                raise problem(
                    f"{self.url}: Content-Range {content_range!r}, "
                    f"expected bytes {self.pos}-<last>/<size>"
                )
            total = int(match.group(3))
        if first:
            self.size = total
            self.etag = etag
            self.last_modified = last_modified
        elif total != self.size:
            raise TarChangedError(
                f"{self.url}: size changed during the build, {self.size} -> {total}"
            )

    def _disconnect(self) -> None:
        if self._resp is not None:
            try:
                self._resp.close()
            except Exception:
                pass
            self._resp = None

    def _retrying(self, action: Any) -> Any:
        while True:
            try:
                result = action()
            except Exception as exc:
                self._disconnect()
                if not _is_retryable(exc):
                    raise
                self._failures += 1
                self.reconnects += 1
                if self._failures > self.max_retries:
                    raise BuildError(
                        f"{self.url}: gave up at byte {self.pos} of {self.size} after "
                        f"{self.max_retries} consecutive failed attempts: {exc!r}"
                    ) from exc
                if self.reconnects > _MAX_TOTAL_RECONNECTS:
                    raise BuildError(
                        f"{self.url}: gave up at byte {self.pos} of {self.size} after "
                        f"{_MAX_TOTAL_RECONNECTS} reconnects"
                    ) from exc
                delay = min(self.backoff * 2 ** (self._failures - 1), _MAX_BACKOFF_SECONDS)
                _log(
                    f"{self.url}: {exc!r} at byte {self.pos}; "
                    f"reconnect {self._failures}/{self.max_retries} in {delay:.1f}s"
                )
                _sleep(delay)
            else:
                return result

    # -- file-like interface used by tarfile

    def _read_once(self, n: int) -> bytes:
        if self._resp is None:
            self._connect()
        assert self._resp is not None
        try:
            chunk = self._resp.read(n)
        except http.client.IncompleteRead as exc:
            chunk = exc.partial
        if not chunk:
            raise ConnectionError(f"connection closed at byte {self.pos} of {self.size}")
        return chunk

    def read(self, n: int = -1) -> bytes:
        remaining = self.size - self.pos
        if remaining <= 0:
            return b""
        want = remaining if n is None or n < 0 else min(n, remaining)
        if want == 0:
            return b""
        chunk = cast(bytes, self._retrying(lambda: self._read_once(want)))
        self._failures = 0
        start = self.pos
        self.pos += len(chunk)
        if chunk[-1]:
            self.last_nonzero_end = self.pos
        else:
            stripped = len(chunk.rstrip(b"\0"))
            if stripped:
                self.last_nonzero_end = start + stripped
        return chunk

    def tell(self) -> int:
        return self.pos

    def close(self) -> None:
        self._disconnect()


def _tar_name(tar: str) -> str:
    if not _TAR_NAME_RE.fullmatch(tar):
        raise ValueError(
            f"not a tar file name: {tar!r} (expected e.g. UP000000429_85962_HELPY_v6.tar)"
        )
    return tar


def _clean(name: str) -> str:
    return name.replace("\t", " ").replace("\r", " ").replace("\n", " ")


def _align_up(n: int) -> int:
    return (n + 511) // 512 * 512


def segment_count(size: int) -> int:
    """Number of segments a tar of ``size`` bytes is built in (1: not split)."""
    if size <= SEGMENT_THRESHOLD_BYTES:
        return 1
    return -(-size // SEGMENT_BYTES)


def segment_bounds(size: int, index: int, count: int) -> tuple[int, int]:
    """Byte range of segment ``index`` (1-based) of ``count`` equal segments."""
    if count < 1 or not 1 <= index <= count:
        raise ValueError(f"no segment {index} of {count}")
    start = min(_align_up(size * (index - 1) // count), size)
    end = size if index == count else min(_align_up(size * index // count), size)
    return start, end


def part_file_name(tar: str, segment: tuple[int, int] | None = None) -> str:
    """``<tar>.tsv.gz`` or, for a segment, ``<tar>.seg03of10.tsv.gz``."""
    if segment is None:
        return tar + _PART_SUFFIX
    index, count = segment
    width = max(2, len(str(count)))
    return f"{tar}.seg{index:0{width}d}of{count:0{width}d}{_PART_SUFFIX}"


def _parse_header(block: bytes) -> tarfile.TarInfo | None:
    """A tar header block with ustar magic and a valid checksum, or None."""
    if len(block) != 512 or block[257:262] != b"ustar":
        return None
    try:
        return tarfile.TarInfo.frombuf(block, "utf-8", "surrogateescape")
    except tarfile.TarError:
        return None


def _is_model(info: tarfile.TarInfo) -> bool:
    return info.isreg() and info.size > 0 and _MEMBER_RE.fullmatch(info.name) is not None


class _Lookahead:
    """Reader over a stream that can look at bytes before handing them on.

    The resynchronization needs to inspect a candidate header, the first bytes
    of its data and the block after the member before the tar parser is
    allowed to start there. Whatever was looked at stays buffered and is
    served first by :meth:`read`.
    """

    def __init__(self, stream: _HTTPTarStream) -> None:
        self.stream = stream
        self.base = stream.pos
        self.buf = bytearray()

    def fill(self, upto: int) -> bool:
        """Buffer up to absolute offset ``upto``. False if the file ends first."""
        while self.base + len(self.buf) < upto:
            chunk = self.stream.read(_BUFSIZE)
            if not chunk:
                return False
            self.buf += chunk
        return True

    def peek(self, offset: int, n: int) -> bytes:
        self.fill(offset + n)
        at = offset - self.base
        return bytes(self.buf[at : at + n])

    def drop_before(self, offset: int) -> None:
        del self.buf[: offset - self.base]
        self.base = offset

    def read(self, n: int = -1) -> bytes:
        if self.buf:
            out = bytes(self.buf if n is None or n < 0 else self.buf[:n])
            del self.buf[: len(out)]
            self.base += len(out)
            return out
        chunk = self.stream.read(n)
        self.base += len(chunk)
        return chunk


def _resync(look: _Lookahead, start: int, end: int, total: int, what: str) -> int | None:
    """Offset of the first model member header at or after ``start``.

    ``start`` is 512-aligned but normally lands inside some member's data.
    Tar headers sit on 512-byte boundaries, so every following block is tried
    as a header. A block is accepted only if all of this holds:

    * ustar magic and a valid header checksum;
    * a regular file whose name matches the model naming;
    * the member's data starts with the gzip magic;
    * the block right after the member is another valid header, the zero
      block that ends the archive, or the end of the file.

    Magic alone would be easy to hit inside compressed data; the checksum,
    the name pattern and the follow-on block together are not. The merge
    still does not take this on trust: it requires the previous segment to
    have stopped at exactly the offset returned here.

    Returns None when the file ends before any such header (the segment
    starts in the last member or in the end-of-archive padding). Gives up
    with :class:`BuildError` if no header shows up within ``_RESYNC_LIMIT``
    bytes past ``end``, so a segment never reads far beyond its range.
    """
    zero = bytes(512)
    pos = start
    while pos + 512 <= total:
        if pos >= end + _RESYNC_LIMIT:
            raise BuildError(
                f"{what}: no member header found between byte {start} and byte {pos}, "
                f"{_RESYNC_LIMIT} bytes past the end of the segment"
            )
        look.drop_before(pos)
        info = _parse_header(look.peek(pos, 512))
        if info is not None and _is_model(info) and info.size <= _RESYNC_LIMIT:
            data = pos + 512
            after = data + _align_up(info.size)
            if after == total or after + 512 <= total:
                follow = look.peek(after, 512) if after < total else zero
                if look.peek(data, 2) == b"\x1f\x8b" and (
                    follow == zero or _parse_header(follow) is not None
                ):
                    return pos
        pos += 512
    return None


def build_part(
    tar: str,
    out_dir: str | os.PathLike[str],
    *,
    base_url: str = V6_BASE_URL,
    segment: tuple[int, int] | None = None,
    bounds: tuple[int, int] | None = None,
    timeout: float = 60.0,
    max_retries: int = 8,
    backoff: float = 2.0,
    progress_seconds: float = 30.0,
) -> Path:
    """Stream one bulk tar, or one segment of it, and write its part file.

    The tar is never written to disk. For every member the part file records
    the offset of its data bytes (``TarInfo.offset_data``) and its size, which
    is exactly the byte range a reader has to request to get that one
    ``.cif.gz`` or ``.pdb.gz``.

    A dropped connection is resumed with a Range request from the byte that
    had been reached (see :class:`_HTTPTarStream`), with at most
    ``max_retries`` consecutive failures and exponential ``backoff``. The part
    is written under a temporary name and renamed at the end, so a file with
    the final name is always complete.

    ``segment=(index, count)`` (1-based) indexes one of ``count`` equal byte
    ranges of the tar instead, so that a 40 GB tar can be built by several
    jobs at once. The output is ``<tar>.seg<index>of<count>.tsv.gz``. A segment
    starts at a 512-aligned offset that is usually inside a member; it finds
    the first model member header at or after that offset (see
    :func:`_resync`), parses from there and stops at the first model member
    header at or past its end. A member therefore belongs to the segment in
    which its header starts, and a segment reads at most one member (plus one
    read buffer) beyond its end. ``bounds=(start, end)`` overrides the byte
    range computed from ``index`` and ``count``; both are rounded up to a
    multiple of 512.
    """
    name = _tar_name(tar)
    url = base_url.rstrip("/") + "/" + name
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    if bounds is not None and segment is None:
        raise ValueError("bounds needs segment=(index, count)")
    final = out_path / part_file_name(name, segment)
    tmp = out_path / f".{final.name}.{os.getpid()}.tmp"
    label = name if segment is None else f"{name} segment {segment[0]}/{segment[1]}"

    start = 0
    seg_end: int | None = None
    if segment is not None:
        known_size = _head(url, timeout=timeout)["size_bytes"]
        if not isinstance(known_size, int):
            raise BuildError(f"{url}: HEAD response has no Content-Length")
        lo, hi = bounds if bounds is not None else segment_bounds(known_size, *segment)
        start, seg_end = _align_up(lo), min(_align_up(hi), known_size)
        if not 0 <= start < seg_end:
            raise ValueError(
                f"{label}: empty byte range {start}-{seg_end} in a tar of {known_size} bytes"
            )

    stream = _HTTPTarStream(
        url, start=start, timeout=timeout, max_retries=max_retries, backoff=backoff
    )
    counts = {fmt: 0 for fmt in _FORMATS}
    n_skipped = 0
    started = time.monotonic()
    last_report = started
    try:
        total = stream.size
        if seg_end is None:
            seg_end = total
        elif total != known_size:
            raise TarChangedError(f"{url}: size changed during the build, {known_size} -> {total}")
        span = seg_end - start
        header: dict[str, Any] = {
            "name": name,
            "url": url,
            "size_bytes": total,
            "last_modified": stream.last_modified,
            "etag": stream.etag,
        }
        if segment is not None:
            header["segment"] = {
                "index": segment[0], "count": segment[1], "start": start, "end": seg_end
            }
        look = _Lookahead(stream)
        with open(tmp, "wb") as raw:
            with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as out:
                out.write(("#" + json.dumps(header, sort_keys=True) + "\n").encode())
                first_header = 0 if start == 0 else _resync(look, start, seg_end, total, label)
                stop: int | None = None
                end_of_members: int | None = None
                if first_header is not None and first_header >= seg_end:
                    # No member header starts inside this segment.
                    stop = first_header
                elif first_header is not None:
                    with tarfile.open(
                        fileobj=cast(IO[bytes], look), mode="r|", bufsize=_BUFSIZE
                    ) as archive:
                        seen: list[tarfile.TarInfo] = getattr(archive, "members")
                        while True:
                            member = archive.next()
                            if member is None:
                                end_of_members = first_header + archive.offset
                                break
                            # tarfile keeps every TarInfo it has seen; with half
                            # a million members that is memory we do not need.
                            seen.clear()
                            # Offsets from tarfile count from where it started.
                            header_at = first_header + member.offset
                            data_at = first_header + member.offset_data
                            match = _MEMBER_RE.fullmatch(member.name) if member.isreg() else None
                            if match is not None and member.size > 0 and header_at >= seg_end:
                                stop = header_at
                                break
                            if match is None:
                                n_skipped += 1
                                out.write(f"#skipped\t{_clean(member.name)}\n".encode())
                            else:
                                if member.size <= 0:
                                    raise BuildError(f"{name}: member {member.name} is empty")
                                if data_at + member.size > total:
                                    raise BuildError(
                                        f"{name}: member {member.name} runs past the end of the tar"
                                    )
                                entry_id, fmt = match.group(1), match.group(2)
                                counts[fmt] += 1
                                out.write(f"{entry_id}\t{fmt}\t{data_at}\t{member.size}\n".encode())
                            now = time.monotonic()
                            if now - last_report >= progress_seconds:
                                last_report = now
                                done = stream.pos - start
                                rate = done / max(now - started, 1e-9) / 1e6
                                _log(
                                    f"{label}: {done / 1e6:,.0f} of {span / 1e6:,.0f} MB "
                                    f"({rate:.1f} MB/s), {sum(counts.values()):,} members, "
                                    f"{stream.reconnects} reconnects"
                                )
                if end_of_members is not None:
                    # tarfile stops quietly at the first block it cannot parse,
                    # so prove that nothing but zero padding follows.
                    while stream.last_nonzero_end <= end_of_members and look.read(1 << 20):
                        pass
                    if stream.last_nonzero_end > end_of_members:
                        raise BuildError(
                            f"{name}: parsing stopped at byte {end_of_members} but the tar has "
                            f"data up to byte {stream.last_nonzero_end} "
                            "(corrupt header or concatenated tar)"
                        )
                    if stream.pos != total:
                        raise BuildError(f"{name}: read up to byte {stream.pos} of {total}")
                trailer: dict[str, Any] = {
                    "n_cif": counts["cif"],
                    "n_pdb": counts["pdb"],
                    "n_skipped": n_skipped,
                    "bytes_read": stream.pos - start,
                    "reconnects": stream.reconnects,
                }
                if segment is None:
                    trailer["end_of_members"] = end_of_members
                else:
                    # first_header: where this segment began parsing (None if
                    # the file ended before any header). stop: the header of
                    # the next segment's first member, or where the archive
                    # ends when reached_end is true.
                    trailer["first_header"] = first_header
                    trailer["stop"] = stop if stop is not None else end_of_members
                    trailer["reached_end"] = end_of_members is not None
                out.write(("#end\t" + json.dumps(trailer, sort_keys=True) + "\n").encode())
            raw.flush()
            os.fsync(raw.fileno())
        os.replace(tmp, final)
    except tarfile.TarError as exc:
        raise BuildError(f"{name}: tar could not be parsed near byte {stream.pos}: {exc}") from exc
    finally:
        stream.close()
        if tmp.exists():
            tmp.unlink()
    elapsed = max(time.monotonic() - started, 1e-9)
    read = stream.pos - start
    _log(
        f"{label}: {counts['cif']:,} cif + {counts['pdb']:,} pdb members, {n_skipped} skipped, "
        f"{read:,} bytes in {elapsed:.0f}s ({read / elapsed / 1e6:.1f} MB/s), "
        f"{stream.reconnects} reconnects -> {final}"
    )
    return final


# --------------------------------------------------------------------------
# Part files


class _Part(NamedTuple):
    header: dict[str, Any]
    rows: list[tuple[str, str, int, int]]
    skipped: list[str]
    trailer: dict[str, Any] | None


def _read_part(path: Path) -> _Part:
    header: dict[str, Any] | None = None
    trailer: dict[str, Any] | None = None
    rows: list[tuple[str, str, int, int]] = []
    skipped: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8", newline="\n") as handle:
        for number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\n")
            if number == 1:
                if not line.startswith("#{"):
                    raise ValueError("first line is not the JSON header")
                header = cast("dict[str, Any]", json.loads(line[1:]))
            elif line.startswith("#skipped\t"):
                skipped.append(line.split("\t", 1)[1])
            elif line.startswith("#end\t"):
                trailer = cast("dict[str, Any]", json.loads(line.split("\t", 1)[1]))
            elif line.startswith("#") or not line:
                continue
            else:
                if trailer is not None:
                    raise ValueError(f"line {number}: row after the #end line")
                fields = line.split("\t")
                if len(fields) != 4 or fields[1] not in _FORMATS:
                    raise ValueError(f"line {number}: malformed row {line!r}")
                rows.append((fields[0], fields[1], int(fields[2]), int(fields[3])))
    if header is None:
        raise ValueError("empty part file")
    return _Part(header, rows, skipped, trailer)


_PART_FILE_RE = re.compile(
    r"(?P<tar>.+?\.tar)(?:\.seg(?P<index>[0-9]+)of(?P<count>[0-9]+))?" + re.escape(_PART_SUFFIX)
)
_Segment = tuple[int, int]


def _find_parts(parts_path: Path) -> tuple[dict[str, dict[_Segment | None, Path]], list[str]]:
    """Part files in a directory as ``{tar: {None or (index, count): path}}``."""
    found: dict[str, dict[_Segment | None, Path]] = {}
    unrecognized: list[str] = []
    for path in sorted(parts_path.glob("*" + _PART_SUFFIX)):
        if path.name.startswith("."):
            continue
        match = _PART_FILE_RE.fullmatch(path.name)
        if match is None:
            unrecognized.append(path.name)
            continue
        segment = (int(match["index"]), int(match["count"])) if match["index"] else None
        found.setdefault(match["tar"], {})[segment] = path
    return found, unrecognized


def _read_named(path: Path) -> _Part:
    try:
        return _read_part(path)
    except (OSError, EOFError, ValueError) as exc:
        raise ValueError(f"unreadable part file {path.name}: {exc}") from exc


def _part_problems(part: _Part, name: str, segment: _Segment | None) -> list[str]:
    """What is wrong with one part file on its own. Empty list: nothing."""
    problems: list[str] = []
    if part.header.get("name") != name:
        problems.append(f"header names {part.header.get('name')!r}")
    declared = part.header.get("segment")
    if segment is None:
        if declared is not None:
            problems.append("a whole-tar part whose header describes a segment")
    elif not isinstance(declared, dict) or (
        declared.get("index"), declared.get("count")
    ) != segment:
        problems.append(f"header describes segment {declared!r}")
    counts = {fmt: sum(1 for row in part.rows if row[1] == fmt) for fmt in _FORMATS}
    if part.trailer is None:
        problems.append("no #end line, the part was cut short")
    elif (
        part.trailer.get("n_cif") != counts["cif"]
        or part.trailer.get("n_pdb") != counts["pdb"]
        or part.trailer.get("n_skipped") != len(part.skipped)
    ):
        problems.append(f"rows do not match the #end line {part.trailer}")
    elif segment is not None:
        first, stop = part.trailer.get("first_header", ""), part.trailer.get("stop", "")
        if (
            not (first is None or isinstance(first, int))
            or not (stop is None or isinstance(stop, int))
            or not isinstance(part.trailer.get("reached_end"), bool)
            or (first is None) != (stop is None)
        ):
            problems.append(f"the #end line lacks the segment offsets: {part.trailer}")
    return problems


def _tiling_problems(segments: list[_Part]) -> list[str]:
    """Prove that consecutive segments cover the tar with no gap or overlap.

    Each segment recorded the offset of the header it started parsing at and
    the offset of the header it stopped at. A segment's stop was reached by
    walking the tar member by member from a known header, so it is a real
    header. Requiring the next segment to have started at exactly that offset
    means its resynchronization found the same header, and by induction from
    offset 0 every segment walked the one true chain of members.
    """
    count = len(segments)
    problems: list[str] = []
    trailers = [cast("dict[str, Any]", part.trailer) for part in segments]
    if trailers[0]["first_header"] != 0:
        problems.append(
            f"segment 1/{count} starts at header offset {trailers[0]['first_header']}, expected 0"
        )
    for index, part in enumerate(segments, 1):
        trailer = trailers[index - 1]
        first, stop = trailer["first_header"], trailer["stop"]
        if first is not None and (stop < first or (trailer["reached_end"] and stop == first)):
            problems.append(f"segment {index}/{count}: stop offset {stop} before start {first}")
        if part.rows and (
            first is None or part.rows[0][2] - 512 < first or part.rows[-1][2] + part.rows[-1][3] > stop
        ):
            problems.append(
                f"segment {index}/{count}: rows lie outside its range, bytes {first} to {stop}"
            )
        if index == count:
            if not (trailer["reached_end"] or first is None):
                problems.append(
                    f"segment {count}/{count} stopped at byte {stop} "
                    "without reaching the end of the archive"
                )
            break
        following = trailers[index]["first_header"]
        pair = f"segment {index}/{count}", f"segment {index + 1}/{count}"
        if trailer["reached_end"] or first is None:
            if following is not None:
                where = f"at byte {stop}" if stop is not None else "before its range"
                problems.append(
                    f"{pair[0]} saw the archive end {where} but {pair[1]} "
                    f"starts at a member header at byte {following}"
                )
        elif following is None:
            problems.append(
                f"{pair[0]} stopped at the member header at byte {stop} but {pair[1]} "
                "found no member header before the end of the tar"
            )
        elif following != stop:
            kind = "a gap" if following > stop else "an overlap"
            problems.append(
                f"{pair[0]} stopped at the member header at byte {stop} but {pair[1]} "
                f"starts at the header at byte {following}: {kind} of {abs(following - stop)} bytes"
            )
    return problems


def _assemble(
    name: str, files: dict[_Segment | None, Path]
) -> tuple[_Part | None, list[str], str]:
    """Load a tar's part, or all of its segments as one part.

    Returns the part (None if it cannot be put together), the problems found
    and a short description for the summary table.
    """
    if None in files:
        if len(files) > 1:
            names = ", ".join(sorted(path.name for path in files.values()))
            return None, [f"both a whole-tar part and segment parts are present: {names}"], "mixed"
        part = _read_named(files[None])
        return part, _part_problems(part, name, None), "whole"
    counts = sorted({segment[1] for segment in files if segment is not None})
    if len(counts) != 1:
        return None, [f"segment parts from different layouts are present: of {counts}"], "mixed"
    count = counts[0]
    label = f"{count} seg"
    wanted: list[_Segment] = [(index, count) for index in range(1, count + 1)]
    missing = [str(segment[0]) for segment in wanted if segment not in files]
    if missing or len(files) != count:
        have = ", ".join(sorted(path.name for path in files.values()))
        return (
            None,
            [f"segment {', '.join(missing) or '?'} of {count} missing (present: {have})"],
            label,
        )
    segments = [_read_named(files[segment]) for segment in wanted]
    problems = [
        f"segment {segment[0]}/{count}: {problem}"
        for segment, part in zip(wanted, segments)
        for problem in _part_problems(part, name, segment)
    ]
    if problems:
        return None, problems, label
    first = segments[0].header
    for index, part in enumerate(segments[1:], 2):
        for key in ("url", "size_bytes", "etag", "last_modified"):
            if part.header.get(key) != first.get(key):
                problems.append(
                    f"segment {index}/{count}: {key} {part.header.get(key)!r} differs from "
                    f"segment 1 ({first.get(key)!r}); the segments are not of the same file"
                )
    problems.extend(_tiling_problems(segments))
    rows = [row for part in segments for row in part.rows]
    skipped = [member for part in segments for member in part.skipped]
    trailer = {
        "n_cif": sum(1 for row in rows if row[1] == "cif"),
        "n_pdb": sum(1 for row in rows if row[1] == "pdb"),
        "n_skipped": len(skipped),
    }
    header = {key: value for key, value in first.items() if key != "segment"}
    return _Part(header, rows, skipped, trailer), problems, label



# --------------------------------------------------------------------------
# download_metadata.json


def _fetch_bytes(url: str, *, timeout: float = 60.0) -> bytes:
    def attempt() -> bytes:
        with _open(_request(url), timeout) as resp:
            return resp.read()

    return cast(bytes, _with_retries(url, attempt, max_retries=4, backoff=2.0))


def fetch_download_metadata(*, timeout: float = 60.0) -> list[dict[str, Any]]:
    """EBI's list of bulk tars, restricted to the v6 ones."""
    data = json.loads(_fetch_bytes(METADATA_URL, timeout=timeout))
    if not isinstance(data, list):
        raise MergeError(f"{METADATA_URL}: expected a JSON list")
    return [d for d in data if str(d.get("archive_name", "")).endswith("_v6.tar")]


def _expected_counts(name: str, n: int) -> dict[str, int]:
    """Members per format that a tar must hold, given num_predicted_structures.

    Measured on real v6 tars (see builders/README.md): a proteome or global
    health tar holds one cif and one pdb per predicted structure; each
    Swiss-Prot tar holds one format only.
    """
    if name in KNOWN_COUNT_EXCEPTIONS:
        return {fmt: KNOWN_COUNT_EXCEPTIONS[name].get(fmt, 0) for fmt in _FORMATS}
    swissprot = _SWISSPROT_TAR_RE.fullmatch(name)
    if swissprot:
        return {fmt: (n if fmt == swissprot.group(1) else 0) for fmt in _FORMATS}
    return {fmt: n for fmt in _FORMATS}


# --------------------------------------------------------------------------
# Merge


def _now_utc() -> str:
    return _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _print_table(rows: Sequence[Sequence[str]]) -> None:
    widths = [max(len(row[i]) for row in rows) for i in range(len(rows[0]))]
    for index, row in enumerate(rows):
        cells = [
            cell.ljust(widths[i]) if i == 0 or i == len(row) - 1 else cell.rjust(widths[i])
            for i, cell in enumerate(row)
        ]
        print("  ".join(cells).rstrip(), flush=True)
        if index == 0:
            print("  ".join("-" * w for w in widths), flush=True)


def merge_parts(
    parts_dir: str | os.PathLike[str],
    out: str | os.PathLike[str],
    *,
    metadata: list[dict[str, Any]] | None = None,
) -> Path:
    """Build ``out`` (index-v6.sqlite) and ``out + ".gz"`` from part files.

    ``metadata`` is the list from download_metadata.json and is fetched from
    EBI when omitted. Pass a filtered list to build an index of a few tars.

    Raises :class:`MergeError`, after printing the per-tar table, when:

    * a tar listed in the metadata has no part file, or a part file is not
      listed in the metadata;
    * a part is unreadable, has no ``#end`` line (the build was cut short) or
      its rows disagree with the counts in that line;
    * a tar has zero rows, the same (entry id, format) twice, or rows that
      are not in increasing, non-overlapping tar order;
    * a tar built in segments (see :func:`build_part`) has a segment missing,
      segments of different versions of the file, or segments that do not
      tile it exactly: segment 1 must start at header offset 0, every
      segment's stop offset must equal the next one's first header offset,
      and the last one must have reached the end of the archive;
    * the number of members per format differs from the tar's
      ``num_predicted_structures`` (both formats for proteome and global
      health tars, the one format for each Swiss-Prot tar);
    * the parts were built from more than one base URL.

    Nothing is written to ``out`` unless every check passes.
    """
    parts_path = Path(parts_dir)
    out_path = Path(out)
    if metadata is None:
        metadata = fetch_download_metadata()
    if not metadata:
        raise MergeError("the tar list is empty")

    found, unrecognized = _find_parts(parts_path)
    listed = [str(entry["archive_name"]) for entry in metadata]
    errors: list[str] = [f"{name}: not a part file name" for name in unrecognized]
    if len(set(listed)) != len(listed):
        errors.append("the tar list names a tar more than once")
    for name in sorted(set(found) - set(listed)):
        errors.append(f"{name}: part file present but the tar is not in the tar list")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_db = out_path.with_name(f".{out_path.name}.{os.getpid()}.tmp")
    gz_path = out_path.with_name(out_path.name + ".gz")
    tmp_gz = out_path.with_name(f".{gz_path.name}.{os.getpid()}.tmp")
    for stale in (tmp_db, tmp_gz):
        if stale.exists():
            stale.unlink()

    table: list[list[str]] = [
        ["tar", "type", "tar bytes", "parts", "expected", "cif", "pdb", "skipped", "status"]
    ]
    base_urls: set[str] = set()
    total_rows = 0
    conn = sqlite3.connect(tmp_db, isolation_level=None)
    try:
        conn.execute("PRAGMA journal_mode=OFF")
        conn.execute("PRAGMA synchronous=OFF")
        conn.executescript(_SCHEMA)
        conn.execute("BEGIN")
        for tar_id, entry in enumerate(metadata, 1):
            name = str(entry["archive_name"])
            tar_type = str(entry.get("type") or "")
            n = int(entry["num_predicted_structures"])
            expected = _expected_counts(name, n)
            expected_text = "/".join(str(expected[fmt]) for fmt in _FORMATS)
            files = found.get(name)
            if not files:
                errors.append(f"{name}: no part file in {parts_path}")
                table.append([name, tar_type, "-", "-", expected_text, "-", "-", "-", "MISSING"])
                continue
            try:
                loaded, problems, parts_label = _assemble(name, files)
            except ValueError as exc:
                errors.append(f"{name}: {exc}")
                table.append(
                    [name, tar_type, "-", "-", expected_text, "-", "-", "-", "UNREADABLE"]
                )
                continue
            if loaded is None:
                errors.extend(f"{name}: {problem}" for problem in problems)
                table.append(
                    [name, tar_type, "-", parts_label, expected_text, "-", "-", "-", "FAILED"]
                )
                continue
            part = loaded
            counts = {fmt: sum(1 for row in part.rows if row[1] == fmt) for fmt in _FORMATS}
            if not part.rows:
                problems.append("zero rows")
            if len({(row[0], row[1]) for row in part.rows}) != len(part.rows):
                problems.append("the same entry id and format appears more than once")
            for before, after in zip(part.rows, part.rows[1:]):
                if before[2] + before[3] > after[2] - 512:
                    problems.append(
                        f"rows are not in tar order without overlap: {before[0]} {before[1]} "
                        f"at bytes {before[2]}+{before[3]}, then {after[0]} {after[1]} at {after[2]}"
                    )
                    break
            for fmt in _FORMATS:
                if counts[fmt] != expected[fmt]:
                    problems.append(
                        f"{counts[fmt]} {fmt} members, expected {expected[fmt]} "
                        f"(num_predicted_structures={n}, difference {counts[fmt] - expected[fmt]:+d})"
                    )
            size = part.header.get("size_bytes")
            url = part.header.get("url")
            if not isinstance(size, int) or not isinstance(url, str):
                problems.append("header lacks size_bytes or url")
            notes = "ok"
            if isinstance(size, int) and "size_bytes" in entry and int(entry["size_bytes"]) != size:
                # Reported, not fatal: the index records what the server sent.
                notes = f"ok (tar list says {int(entry['size_bytes'])} bytes)"
            table.append(
                [
                    name,
                    tar_type,
                    f"{size:,}" if isinstance(size, int) else "-",
                    parts_label,
                    expected_text,
                    str(counts["cif"]),
                    str(counts["pdb"]),
                    str(len(part.skipped)),
                    notes if not problems else "FAILED",
                ]
            )
            if problems:
                errors.extend(f"{name}: {problem}" for problem in problems)
                continue
            assert isinstance(url, str) and isinstance(size, int)
            base_urls.add(url.rsplit("/", 1)[0])
            proteome = _PROTEOME_TAR_RE.fullmatch(name)
            conn.execute(
                "INSERT INTO tars VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                (
                    tar_id,
                    name,
                    url,
                    size,
                    part.header.get("last_modified"),
                    part.header.get("etag"),
                    tar_type or None,
                    entry.get("reference_proteome") or (proteome.group(1) if proteome else None),
                    int(proteome.group(2)) if proteome else None,
                    entry.get("species"),
                    entry.get("common_name"),
                    len(part.rows),
                ),
            )
            conn.executemany(
                "INSERT INTO members VALUES (?,?,?,?,?)",
                ((row[0], row[1], tar_id, row[2], row[3]) for row in part.rows),
            )
            total_rows += len(part.rows)

        _print_table(table)
        if len(base_urls) > 1:
            errors.append(f"parts were built from different base URLs: {sorted(base_urls)}")
        if errors:
            raise MergeError(
                f"index not written, {len(errors)} problem(s):\n  " + "\n  ".join(errors)
            )
        meta_rows = {
            "schema_version": SCHEMA_VERSION,
            "afdb_version": AFDB_VERSION,
            "built_at": _now_utc(),
            "base_url": next(iter(base_urls)),
        }
        conn.executemany("INSERT INTO meta VALUES (?,?)", sorted(meta_rows.items()))
        conn.execute("COMMIT")
        conn.execute("VACUUM")
        conn.execute("ANALYZE")
        conn.close()
        with open(tmp_db, "rb") as src, open(tmp_gz, "wb") as raw:
            with gzip.GzipFile(filename="", fileobj=raw, mode="wb", compresslevel=9, mtime=0) as gz:
                shutil.copyfileobj(src, gz, 1 << 20)
            raw.flush()
            os.fsync(raw.fileno())
        os.replace(tmp_db, out_path)
        os.replace(tmp_gz, gz_path)
    finally:
        conn.close()
        for leftover in (tmp_db, tmp_gz):
            if leftover.exists():
                leftover.unlink()
    print(
        f"{out_path}: {len(metadata)} tars, {total_rows:,} rows, "
        f"{out_path.stat().st_size:,} bytes ({gz_path.stat().st_size:,} gzipped)",
        flush=True,
    )
    return out_path


# --------------------------------------------------------------------------
# Spot check against the live server


def _check_model_text(entry_id: str, fmt: str, text: str) -> str | None:
    """Return None when ``text`` is a plausible model for ``entry_id``."""
    match = _ENTRY_RE.fullmatch(entry_id)
    if match is None:
        return "entry id does not look like AF-<accession>-F<n>"
    accession = match.group(1)
    if fmt == "cif":
        first_line = text.split("\n", 1)[0].strip()
        if first_line != f"data_{entry_id}":
            return f"first line is {first_line[:60]!r}, expected 'data_{entry_id}'"
        if "_atom_site." not in text:
            return "no _atom_site category"
        return None
    atom_at = text.find("\nATOM  ")
    if atom_at < 0:
        return "no ATOM records"
    if not text.startswith("HEADER"):
        return "does not start with a HEADER record"
    # AlphaFold PDB files carry the accession (not the entry id) in the TITLE
    # and DBREF records.
    if not re.search(rf"(?<![A-Z0-9]){re.escape(accession)}(?![A-Z0-9])", text[:atom_at]):
        return f"accession {accession} not found in the header records"
    return None


def _range_read(
    url: str, offset: int, length: int, *, timeout: float
) -> tuple[bytes, str | None, str | None]:
    def attempt() -> tuple[bytes, str | None, str | None]:
        request = _request(url, headers={"Range": f"bytes={offset}-{offset + length - 1}"})
        with _open(request, timeout) as resp:
            if resp.status != 206:
                raise VerifyError(f"{url}: HTTP {resp.status} for a range request, expected 206")
            return resp.read(), resp.headers.get("ETag"), resp.headers.get("Last-Modified")

    return cast(
        "tuple[bytes, str | None, str | None]",
        _with_retries(url, attempt, max_retries=3, backoff=2.0),
    )


def verify_index(
    index_path: str | os.PathLike[str],
    n: int = 20,
    *,
    per_tar: int = 0,
    seed: int | None = None,
    timeout: float = 60.0,
) -> dict[str, Any]:
    """Read random members back from the server and check what comes out.

    Picks ``n`` members at random from the whole index, plus ``per_tar`` more
    from every tar. For each one it sends a single Range request for
    ``offset .. offset+length-1``, requires a 206 whose ETag and Last-Modified
    match the index, gunzips the bytes and checks the text: a cif must start
    with ``data_<entry id>`` and contain ``_atom_site``; a pdb must start with
    HEADER, name the accession in its header records and have ATOM records.

    Prints one line per member and a summary, returns the summary as a dict
    and raises :class:`VerifyError` if any member fails.
    """
    rng = random.Random(seed)
    conn = sqlite3.connect(f"file:{Path(index_path).resolve()}?mode=ro", uri=True)
    try:
        tars = {
            int(row[0]): (str(row[1]), str(row[2]), row[3], row[4])
            for row in conn.execute("SELECT tar_id, name, url, etag, last_modified FROM tars")
        }
        total = int(conn.execute("SELECT count(*) FROM members").fetchone()[0])
        if total == 0:
            raise VerifyError(f"{index_path}: the index has no members")
        chosen: set[tuple[str, str, int, int, int]] = set()

        def pick(tar_id: int | None, count: int) -> None:
            where = "" if tar_id is None else "WHERE tar_id = ?"
            args: tuple[int, ...] = () if tar_id is None else (tar_id,)
            (available,) = conn.execute(f"SELECT count(*) FROM members {where}", args).fetchone()
            for position in rng.sample(range(int(available)), min(count, int(available))):
                row = conn.execute(
                    f"SELECT entry_id, fmt, tar_id, offset, length FROM members {where} "
                    "ORDER BY entry_id, fmt, tar_id LIMIT 1 OFFSET ?",
                    (*args, position),
                ).fetchone()
                chosen.add((str(row[0]), str(row[1]), int(row[2]), int(row[3]), int(row[4])))

        pick(None, n)
        if per_tar > 0:
            for tar_id in sorted(tars):
                pick(tar_id, per_tar)
    finally:
        conn.close()

    failures: list[dict[str, Any]] = []
    bytes_read = 0
    by_format = {fmt: 0 for fmt in _FORMATS}
    for entry_id, fmt, tar_id, offset, length in sorted(chosen, key=lambda c: (c[2], c[3])):
        name, url, etag, last_modified = tars[tar_id]
        problem: str | None
        try:
            blob, got_etag, got_modified = _range_read(url, offset, length, timeout=timeout)
            bytes_read += len(blob)
            if len(blob) != length:
                problem = f"got {len(blob)} bytes, expected {length}"
            elif etag is not None and got_etag != etag:
                problem = f"ETag is {got_etag!r}, the index was built against {etag!r}"
            elif last_modified is not None and got_modified != last_modified:
                problem = f"Last-Modified is {got_modified!r}, the index has {last_modified!r}"
            else:
                text = gzip.decompress(blob).decode("utf-8", errors="replace")
                problem = _check_model_text(entry_id, fmt, text)
        except VerifyError as exc:
            problem = str(exc)
        except (OSError, EOFError, http.client.HTTPException) as exc:
            problem = f"{type(exc).__name__}: {exc}"
        by_format[fmt] += 1
        status = "ok" if problem is None else f"FAIL: {problem}"
        print(f"{entry_id:<22} {fmt}  {name}  bytes {offset}+{length}  {status}", flush=True)
        if problem is not None:
            failures.append(
                {"entry_id": entry_id, "fmt": fmt, "tar": name, "offset": offset,
                 "length": length, "problem": problem}
            )
    summary: dict[str, Any] = {
        "checked": len(chosen),
        "passed": len(chosen) - len(failures),
        "failed": failures,
        "by_format": by_format,
        "bytes_read": bytes_read,
        "index_members": total,
    }
    print(
        f"verify_index: {summary['passed']}/{summary['checked']} members read back correctly "
        f"({by_format['cif']} cif, {by_format['pdb']} pdb, {bytes_read:,} bytes)",
        flush=True,
    )
    if failures:
        raise VerifyError(
            f"{len(failures)} of {len(chosen)} spot checks failed: "
            + "; ".join(f"{f['entry_id']} {f['fmt']}: {f['problem']}" for f in failures[:5])
        )
    return summary


# --------------------------------------------------------------------------
# Published summary and the "is a rebuild needed" check


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_summary(
    index_path: str | os.PathLike[str], out: str | os.PathLike[str]
) -> dict[str, Any]:
    """Write the small ``index-v6.json`` that is published next to the index."""
    index = Path(index_path)
    gz_path = index.with_name(index.name + ".gz")
    conn = sqlite3.connect(f"file:{index.resolve()}?mode=ro", uri=True)
    try:
        meta = {str(k): str(v) for k, v in conn.execute("SELECT key, value FROM meta")}
        tars = [
            {"name": row[0], "size_bytes": row[1], "last_modified": row[2], "etag": row[3],
             "n_members": row[4]}
            for row in conn.execute(
                "SELECT name, size_bytes, last_modified, etag, n_members FROM tars ORDER BY tar_id"
            )
        ]
        total = int(conn.execute("SELECT count(*) FROM members").fetchone()[0])
    finally:
        conn.close()
    summary: dict[str, Any] = {
        "schema_version": meta.get("schema_version"),
        "afdb_version": meta.get("afdb_version"),
        "built_at": meta.get("built_at"),
        "base_url": meta.get("base_url"),
        "total_rows": total,
        "index_bytes": index.stat().st_size,
        "index_gz_bytes": gz_path.stat().st_size,
        "index_gz_sha256": _sha256(gz_path),
        "tars": tars,
    }
    target = Path(out)
    tmp = target.with_name(f".{target.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, target)
    return summary


def _head(url: str, *, timeout: float) -> dict[str, Any]:
    def attempt() -> dict[str, Any]:
        with _open(_request(url, method="HEAD"), timeout) as resp:
            length = resp.headers.get("Content-Length")
            return {
                "size_bytes": int(length) if length is not None else None,
                "last_modified": resp.headers.get("Last-Modified"),
                "etag": resp.headers.get("ETag"),
            }

    return cast("dict[str, Any]", _with_retries(url, attempt, max_retries=4, backoff=2.0))


def needs_rebuild(
    published: dict[str, Any] | None,
    *,
    metadata: list[dict[str, Any]] | None = None,
    base_url: str = V6_BASE_URL,
    timeout: float = 60.0,
) -> list[str]:
    """Compare the published summary with the server. Empty list: up to date.

    Sends one HEAD request per tar and returns the reasons a rebuild is
    needed: no published summary, another schema version, a different set of
    tars, or a tar whose size, Last-Modified or ETag changed.
    """
    if not published:
        return ["no published index-v6.json"]
    if str(published.get("schema_version")) != SCHEMA_VERSION:
        return [f"schema_version {published.get('schema_version')!r} != {SCHEMA_VERSION!r}"]
    if metadata is None:
        metadata = fetch_download_metadata(timeout=timeout)
    known = {str(t.get("name")): t for t in published.get("tars", [])}
    listed = [str(entry["archive_name"]) for entry in metadata]
    reasons = [f"{name}: not in the published index" for name in listed if name not in known]
    reasons += [f"{name}: no longer listed by EBI" for name in sorted(set(known) - set(listed))]
    for name in listed:
        if name not in known:
            continue
        now = _head(base_url.rstrip("/") + "/" + name, timeout=timeout)
        for key in ("size_bytes", "last_modified", "etag"):
            if now[key] != known[name].get(key):
                reasons.append(f"{name}: {key} {known[name].get(key)!r} -> {now[key]!r}")
    return reasons


# --------------------------------------------------------------------------
# Command line used by the workflow: python -m scigantic_alphafold._build ...


def _github_output(values: dict[str, str]) -> None:
    target = os.environ.get("GITHUB_OUTPUT")
    if target:
        with open(target, "a", encoding="utf-8") as handle:
            for key, value in values.items():
                handle.write(f"{key}={value}\n")


def plan_jobs(
    metadata: list[dict[str, Any]],
    *,
    base_url: str = V6_BASE_URL,
    reuse_dir: str | os.PathLike[str] | None = None,
    reuse_out: str | os.PathLike[str] | None = None,
    timeout: float = 60.0,
) -> tuple[list[dict[str, str]], list[str]]:
    """The build jobs for a full index: one per tar, or one per segment.

    A tar larger than ``SEGMENT_THRESHOLD_BYTES`` (by the size in the tar
    list) is split into segments of at most ``SEGMENT_BYTES``. Each job is
    ``{"tar": ..., "segment": "" or "3/10", "part": part file name without
    .tsv.gz}``.

    ``reuse_dir`` holds part files of an earlier run (searched recursively).
    A part found there is reused, and its job left out, when it is complete
    and was built from the file the server has now (same size, ETag and
    Last-Modified, one HEAD request per tar). A complete whole-tar part
    stands in for all segments of that tar. Reused parts are copied to
    ``reuse_out`` and their file names returned.
    """
    candidates: dict[str, Path] = {}
    if reuse_dir is not None:
        for path in sorted(Path(reuse_dir).rglob("*" + _PART_SUFFIX)):
            candidates.setdefault(path.name, path)
    jobs: list[dict[str, str]] = []
    reused: list[str] = []
    for entry in metadata:
        name = _tar_name(str(entry["archive_name"]))
        count = segment_count(int(entry["size_bytes"]))
        segments: list[_Segment | None] = (
            [None] if count == 1 else [(index, count) for index in range(1, count + 1)]
        )
        head: dict[str, Any] | None = None

        def reusable(segment: _Segment | None) -> bool:
            nonlocal head
            file_name = part_file_name(name, segment)
            path = candidates.get(file_name)
            if path is None:
                return False
            if head is None:
                head = _head(base_url.rstrip("/") + "/" + name, timeout=timeout)
            try:
                part = _read_named(path)
            except ValueError as exc:
                print(f"not reusing {file_name}: {exc}")
                return False
            problems = _part_problems(part, name, segment)
            problems += [
                f"{key} was {part.header.get(key)!r}, the server now has {head[key]!r}"
                for key in ("size_bytes", "etag", "last_modified")
                if part.header.get(key) != head[key]
            ]
            if problems:
                print(f"not reusing {file_name}: {'; '.join(problems)}")
                return False
            if reuse_out is not None:
                Path(reuse_out).mkdir(parents=True, exist_ok=True)
                shutil.copyfile(path, Path(reuse_out) / file_name)
            reused.append(file_name)
            return True

        if reusable(None):
            continue
        for segment in segments:
            if segment is not None and reusable(segment):
                continue
            file_name = part_file_name(name, segment)
            jobs.append(
                {
                    "tar": name,
                    "segment": "" if segment is None else f"{segment[0]}/{segment[1]}",
                    "part": file_name[: -len(_PART_SUFFIX)],
                }
            )
    return jobs, reused


def _cmd_plan(args: argparse.Namespace) -> int:
    metadata = fetch_download_metadata()
    names = [str(entry["archive_name"]) for entry in metadata]
    published: dict[str, Any] | None = None
    if args.published and Path(args.published).is_file():
        published = json.loads(Path(args.published).read_text(encoding="utf-8"))
    if args.scheduled and not published:
        # Without a published index there is nothing to keep fresh, and a
        # schedule that started a 153 GB build every week until one succeeds
        # would be rude to EBI. The first build is always started by hand.
        reasons = []
        print("no published index-v6.json: the scheduled run does not start the first build")
    elif args.force:
        reasons = ["forced"]
    else:
        reasons = needs_rebuild(published, metadata=metadata, base_url=args.base_url)
        if not reasons:
            print(f"up to date: all {len(names)} tars match the published index")
    for reason in reasons:
        print(f"rebuild: {reason}")
    print(f"{len(names)} tars, {sum(int(e['size_bytes']) for e in metadata):,} bytes")
    jobs: list[dict[str, str]] = []
    reused: list[str] = []
    if reasons:
        jobs, reused = plan_jobs(
            metadata, base_url=args.base_url, reuse_dir=args.reuse_dir, reuse_out=args.reuse_out
        )
        sizes = {str(e["archive_name"]): int(e["size_bytes"]) for e in metadata}
        largest = 0
        for job in jobs:
            size = sizes[job["tar"]]
            if job["segment"]:
                index, count = (int(x) for x in job["segment"].split("/"))
                start, end = segment_bounds(size, index, count)
                size = end - start
            largest = max(largest, size)
            print(f"job: {job['tar']} {job['segment'] or 'whole'} {size:,} bytes")
        for file_name in reused:
            print(f"reused: {file_name}")
        print(
            f"{len(jobs)} build jobs (largest {largest:,} bytes), "
            f"{len(reused)} parts reused from an earlier run"
        )
    _github_output(
        {
            "rebuild": "true" if reasons else "false",
            "jobs": json.dumps(jobs, separators=(",", ":")),
            "n_jobs": str(len(jobs)),
        }
    )
    return 0


def _parse_segment(text: str | None) -> _Segment | None:
    if not text:
        return None
    match = re.fullmatch(r"([0-9]+)/([0-9]+)", text)
    if match is None:
        raise ValueError(f"not a segment: {text!r} (expected e.g. 3/10)")
    return int(match.group(1)), int(match.group(2))


def _filtered_metadata(only: Iterable[str] | None) -> list[dict[str, Any]] | None:
    if not only:
        return None
    wanted = [_tar_name(name) for name in only]
    by_name = {str(e["archive_name"]): e for e in fetch_download_metadata()}
    missing = [name for name in wanted if name not in by_name]
    if missing:
        raise MergeError(f"not in download_metadata.json: {', '.join(missing)}")
    return [by_name[name] for name in wanted]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m scigantic_alphafold._build",
        description="Build and check the AlphaFold DB v6 offset index.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("plan", help="list the tars and decide whether a rebuild is needed")
    p.add_argument("--published", help="path of the currently published index-v6.json, if any")
    p.add_argument("--force", action="store_true", help="rebuild whatever the server says")
    p.add_argument(
        "--scheduled",
        action="store_true",
        help="do nothing when no index has been published yet",
    )
    p.add_argument("--reuse-dir", help="part files of an earlier run that may be reused")
    p.add_argument("--reuse-out", help="where to copy the parts that are reused")
    p.add_argument("--base-url", default=V6_BASE_URL)

    p = sub.add_parser("build-part", help="stream one tar and write its part file")
    p.add_argument("tar")
    p.add_argument("-o", "--out-dir", default="parts")
    p.add_argument("--segment", metavar="K/N", help="build segment K of N only, e.g. 3/10")
    p.add_argument("--base-url", default=V6_BASE_URL)

    p = sub.add_parser("merge", help="merge part files into index-v6.sqlite and .gz")
    p.add_argument("parts_dir")
    p.add_argument("-o", "--out", default="index-v6.sqlite")
    p.add_argument("--only", nargs="+", metavar="TAR", help="build a partial index of these tars")

    p = sub.add_parser("verify", help="read random members back from the server")
    p.add_argument("index")
    p.add_argument("-n", type=int, default=20)
    p.add_argument("--per-tar", type=int, default=0)
    p.add_argument("--seed", type=int)

    p = sub.add_parser("summary", help="write index-v6.json for a built index")
    p.add_argument("index")
    p.add_argument("-o", "--out", default="index-v6.json")

    args = parser.parse_args(argv)
    if args.command == "plan":
        return _cmd_plan(args)
    if args.command == "build-part":
        segment = _parse_segment(args.segment)
        print(build_part(args.tar, args.out_dir, base_url=args.base_url, segment=segment))
    elif args.command == "merge":
        print(merge_parts(args.parts_dir, args.out, metadata=_filtered_metadata(args.only)))
    elif args.command == "verify":
        verify_index(args.index, args.n, per_tar=args.per_tar, seed=args.seed)
    elif args.command == "summary":
        print(json.dumps(write_summary(args.index, args.out), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
