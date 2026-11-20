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
        self, url: str, *, timeout: float = 60.0, max_retries: int = 8, backoff: float = 2.0
    ) -> None:
        self.url = url
        self.timeout = timeout
        self.max_retries = max_retries
        self.backoff = backoff
        self.pos = 0
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
            if self.etag and not self.etag.startswith("W/"):
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
        if self.pos == 0 and self.size < 0:
            if resp.status != 200:
                raise BuildError(f"{self.url}: expected HTTP 200, got {resp.status}")
            length = resp.headers.get("Content-Length")
            if length is None:
                raise BuildError(f"{self.url}: response has no Content-Length")
            self.size = int(length)
            self.etag = etag
            self.last_modified = last_modified
            return
        if etag != self.etag or last_modified != self.last_modified:
            raise TarChangedError(
                f"{self.url} changed during the build: ETag {self.etag!r} -> {etag!r}, "
                f"Last-Modified {self.last_modified!r} -> {last_modified!r}"
            )
        if self.pos == 0:
            # Reconnect before the first byte was read: a plain 200 again.
            if resp.status != 200 or int(resp.headers.get("Content-Length", "-1")) != self.size:
                raise TarChangedError(f"{self.url}: size changed during the build")
            return
        if resp.status != 206:
            raise TarChangedError(
                f"{self.url}: asked for bytes from {self.pos} and got HTTP {resp.status} "
                "instead of 206; the server no longer serves the same file by range"
            )
        content_range = resp.headers.get("Content-Range", "")
        expected = f"bytes {self.pos}-{self.size - 1}/{self.size}"
        if content_range.strip() != expected:
            raise TarChangedError(
                f"{self.url}: Content-Range {content_range!r}, expected {expected!r}"
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


def build_part(
    tar: str,
    out_dir: str | os.PathLike[str],
    *,
    base_url: str = V6_BASE_URL,
    timeout: float = 60.0,
    max_retries: int = 8,
    backoff: float = 2.0,
    progress_seconds: float = 30.0,
) -> Path:
    """Stream one bulk tar and write ``<out_dir>/<tar>.tsv.gz``.

    The tar is never written to disk. For every member the part file records
    the offset of its data bytes (``TarInfo.offset_data``) and its size, which
    is exactly the byte range a reader has to request to get that one
    ``.cif.gz`` or ``.pdb.gz``.

    A dropped connection is resumed with a Range request from the byte that
    had been reached (see :class:`_HTTPTarStream`), with at most
    ``max_retries`` consecutive failures and exponential ``backoff``. The part
    is written under a temporary name and renamed at the end, so a file named
    ``<tar>.tsv.gz`` is always complete.
    """
    name = _tar_name(tar)
    url = base_url.rstrip("/") + "/" + name
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    final = out_path / (name + _PART_SUFFIX)
    tmp = out_path / f".{name}{_PART_SUFFIX}.{os.getpid()}.tmp"

    stream = _HTTPTarStream(url, timeout=timeout, max_retries=max_retries, backoff=backoff)
    counts = {fmt: 0 for fmt in _FORMATS}
    n_skipped = 0
    started = time.monotonic()
    last_report = started
    try:
        header = {
            "name": name,
            "url": url,
            "size_bytes": stream.size,
            "last_modified": stream.last_modified,
            "etag": stream.etag,
        }
        with open(tmp, "wb") as raw:
            with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as out:
                out.write(("#" + json.dumps(header, sort_keys=True) + "\n").encode())
                with tarfile.open(
                    fileobj=cast(IO[bytes], stream), mode="r|", bufsize=_BUFSIZE
                ) as archive:
                    seen: list[tarfile.TarInfo] = getattr(archive, "members")
                    while True:
                        member = archive.next()
                        if member is None:
                            break
                        # tarfile keeps every TarInfo it has seen; with half a
                        # million members that is memory we do not need.
                        seen.clear()
                        match = _MEMBER_RE.fullmatch(member.name) if member.isreg() else None
                        if match is None:
                            n_skipped += 1
                            out.write(f"#skipped\t{_clean(member.name)}\n".encode())
                        else:
                            if member.size <= 0:
                                raise BuildError(f"{name}: member {member.name} is empty")
                            if member.offset_data + member.size > stream.size:
                                raise BuildError(
                                    f"{name}: member {member.name} runs past the end of the tar"
                                )
                            entry_id, fmt = match.group(1), match.group(2)
                            counts[fmt] += 1
                            out.write(
                                f"{entry_id}\t{fmt}\t{member.offset_data}\t{member.size}\n".encode()
                            )
                        now = time.monotonic()
                        if now - last_report >= progress_seconds:
                            last_report = now
                            rate = stream.pos / max(now - started, 1e-9) / 1e6
                            _log(
                                f"{name}: {stream.pos / 1e6:,.0f} of {stream.size / 1e6:,.0f} MB "
                                f"({rate:.1f} MB/s), {sum(counts.values()):,} members, "
                                f"{stream.reconnects} reconnects"
                            )
                    end_of_members = archive.offset
                # tarfile stops quietly at the first block it cannot parse, so
                # prove that nothing but zero padding follows the last member.
                while stream.last_nonzero_end <= end_of_members and stream.read(1 << 20):
                    pass
                if stream.last_nonzero_end > end_of_members:
                    raise BuildError(
                        f"{name}: parsing stopped at byte {end_of_members} but the tar has data "
                        f"up to byte {stream.last_nonzero_end} (corrupt header or concatenated tar)"
                    )
                if stream.pos != stream.size:
                    raise BuildError(f"{name}: read {stream.pos} bytes of {stream.size}")
                trailer = {
                    "n_cif": counts["cif"],
                    "n_pdb": counts["pdb"],
                    "n_skipped": n_skipped,
                    "bytes_read": stream.pos,
                    "end_of_members": end_of_members,
                    "reconnects": stream.reconnects,
                }
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
    _log(
        f"{name}: {counts['cif']:,} cif + {counts['pdb']:,} pdb members, {n_skipped} skipped, "
        f"{stream.size:,} bytes in {elapsed:.0f}s ({stream.size / elapsed / 1e6:.1f} MB/s), "
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
    * a part has zero rows, or the same (entry id, format) twice;
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

    available = {
        p.name[: -len(_PART_SUFFIX)]: p
        for p in sorted(parts_path.glob("*" + _PART_SUFFIX))
        if not p.name.startswith(".")
    }
    listed = [str(entry["archive_name"]) for entry in metadata]
    errors: list[str] = []
    if len(set(listed)) != len(listed):
        errors.append("the tar list names a tar more than once")
    for name in sorted(set(available) - set(listed)):
        errors.append(f"{name}: part file present but the tar is not in the tar list")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_db = out_path.with_name(f".{out_path.name}.{os.getpid()}.tmp")
    gz_path = out_path.with_name(out_path.name + ".gz")
    tmp_gz = out_path.with_name(f".{gz_path.name}.{os.getpid()}.tmp")
    for stale in (tmp_db, tmp_gz):
        if stale.exists():
            stale.unlink()

    table: list[list[str]] = [
        ["tar", "type", "tar bytes", "expected", "cif", "pdb", "skipped", "status"]
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
            path = available.get(name)
            if path is None:
                errors.append(f"{name}: no part file in {parts_path}")
                table.append([name, tar_type, "-", expected_text, "-", "-", "-", "MISSING"])
                continue
            try:
                part = _read_part(path)
            except (OSError, EOFError, ValueError) as exc:
                errors.append(f"{name}: unreadable part file {path.name}: {exc}")
                table.append([name, tar_type, "-", expected_text, "-", "-", "-", "UNREADABLE"])
                continue
            counts = {fmt: sum(1 for row in part.rows if row[1] == fmt) for fmt in _FORMATS}
            problems: list[str] = []
            if part.header.get("name") != name:
                problems.append(f"header names {part.header.get('name')!r}")
            if part.trailer is None:
                problems.append("no #end line, the part was cut short")
            elif (
                part.trailer.get("n_cif") != counts["cif"]
                or part.trailer.get("n_pdb") != counts["pdb"]
                or part.trailer.get("n_skipped") != len(part.skipped)
            ):
                problems.append(f"rows do not match the #end line {part.trailer}")
            if not part.rows:
                problems.append("zero rows")
            if len({(row[0], row[1]) for row in part.rows}) != len(part.rows):
                problems.append("the same entry id and format appears more than once")
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


def _cmd_plan(args: argparse.Namespace) -> int:
    metadata = fetch_download_metadata()
    names = [str(entry["archive_name"]) for entry in metadata]
    if args.force:
        reasons = ["forced"]
    else:
        published: dict[str, Any] | None = None
        if args.published and Path(args.published).is_file():
            published = json.loads(Path(args.published).read_text(encoding="utf-8"))
        reasons = needs_rebuild(published, metadata=metadata, base_url=args.base_url)
    for reason in reasons:
        print(f"rebuild: {reason}")
    if not reasons:
        print(f"up to date: all {len(names)} tars match the published index")
    print(f"{len(names)} tars, {sum(int(e['size_bytes']) for e in metadata):,} bytes")
    _github_output(
        {"rebuild": "true" if reasons else "false", "tars": json.dumps(names, separators=(",", ":"))}
    )
    return 0


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
    p.add_argument("--base-url", default=V6_BASE_URL)

    p = sub.add_parser("build-part", help="stream one tar and write its part file")
    p.add_argument("tar")
    p.add_argument("-o", "--out-dir", default="parts")
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
        print(build_part(args.tar, args.out_dir, base_url=args.base_url))
    elif args.command == "merge":
        print(merge_parts(args.parts_dir, args.out, metadata=_filtered_metadata(args.only)))
    elif args.command == "verify":
        verify_index(args.index, args.n, per_tar=args.per_tar, seed=args.seed)
    elif args.command == "summary":
        print(json.dumps(write_summary(args.index, args.out), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
