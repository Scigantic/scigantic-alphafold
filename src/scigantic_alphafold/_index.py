"""Locate, cache and query the offset index.

The index is one SQLite file that maps (entry id, format) to (tar, byte offset, length).
Nothing here runs at import time: the index is resolved, downloaded if needed and opened
on the first call that needs it.

Where the index comes from
--------------------------
1. ``SCIGANTIC_ALPHAFOLD_INDEX``, when set: a local path (absolute, relative to the
   current directory, or starting with ``~``) or an http(s) URL, to either a ``.sqlite``
   file or a gzipped ``.sqlite.gz`` (told apart by content, not by name).
2. Otherwise the default release URL, :data:`DEFAULT_INDEX_URL`. The release also carries
   a small ``index-v6.json`` with the SHA-256 of the gzipped index. A download from the
   default URL is compared with it and refused on a mismatch; when that JSON is absent
   (HTTP 404) or has no checksum field the download is used unverified.

A plain local ``.sqlite`` is opened in place. Everything else is materialized once under
the cache directory (``SCIGANTIC_ALPHAFOLD_CACHE``, default ``~/.cache/scigantic-alphafold``):
written to a unique temporary name in that directory, gunzipped, checked, then renamed
over the final name with ``os.replace``. Two processes that race each do their own
download and the last rename wins; a reader never sees a half-written file. Temporary
files left by a killed process are removed once they are a day old.

What the index is trusted for
-----------------------------
The index decides which URL is requested and what a downloaded file is called, so its
content is checked before use. ``meta.afdb_version`` must be digits. Every tar URL must
be http(s) and lie under ``meta.base_url``, and ``base_url`` must be an https URL on
``ftp.ebi.ac.uk``. Entry ids, formats, offsets and lengths must be well formed and inside
the tar. To use a private mirror or a test server, set
``SCIGANTIC_ALPHAFOLD_ALLOW_ANY_HOST=1``: that lifts the https and EBI host requirements,
nothing else.

How ids are resolved
--------------------
Accepted forms, in any letter case and with surrounding whitespace: a UniProt accession
(``P69905``), an isoform accession (``P12345-2``), an entry id (``AF-P69905-F1``,
``AF-P12345-2-F1``), the same without the ``AF-`` prefix (``P69905-F1``), and a member
file name (``AF-P69905-F1-model_v6.cif.gz``, ``AF-P69905-F1-model_v6.pdb``).

A file name carries a version and usually a format, and both are honored:

* The version must be the one the index describes. ``...-model_v4.cif`` against a v6
  index raises :class:`InvalidIdError` instead of quietly returning the v6 file.
* When ``fmt`` is not passed, the extension decides (``.pdb.gz`` reads the PDB file).
  With no extension either, the format is ``cif``.
* When ``fmt`` is passed and disagrees with the extension, :class:`InvalidIdError`.

* An id that names a fragment (``-F<n>``) resolves to exactly that entry.
* An accession without a fragment resolves to ``AF-<accession>-F1`` when that is the only
  fragment of that accession. Isoform entries (``AF-<accession>-2-F1``) never compete
  with the canonical entry: they are returned only when the isoform accession is asked for.
* When the accession has several fragments (F1..Fn), or has no entry of its own but
  isoform entries exist, :class:`AmbiguousEntryError` lists the candidates. Nothing is
  picked silently, because F1 of a fragmented prediction is only part of the protein.
* :func:`locate` never raises for these cases: it returns every row for the accession,
  isoforms and fragments included.

The same entry in several tars
------------------------------
A Swiss-Prot protein of a proteome organism is stored both in the Swiss-Prot tar and in
the proteome tar. The copies are the same file. The rule is: the tar whose name sorts
first in plain byte order wins. Proteome tars are named ``UP...`` and the Swiss-Prot tars
``swissprot_...``, and uppercase sorts before lowercase, so the proteome copy is used.
The rule depends only on tar names, so it gives the same answer for every index build.
:func:`locate` returns all copies, the chosen one first.
"""

from __future__ import annotations

import gzip
import hashlib
import http.client
import json
import os
import re
import shutil
import sqlite3
import tempfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ._errors import (
    AlphaFoldError,
    AmbiguousEntryError,
    IndexSchemaError,
    IndexUnavailableError,
    InvalidIdError,
    NotInBulkSetError,
)
from ._version import USER_AGENT

DEFAULT_INDEX_URL = (
    "https://github.com/Scigantic/scigantic-alphafold/releases/download/index-v6/index-v6.sqlite.gz"
)
ENV_INDEX = "SCIGANTIC_ALPHAFOLD_INDEX"
ENV_CACHE = "SCIGANTIC_ALPHAFOLD_CACHE"
ENV_ALLOW_ANY_HOST = "SCIGANTIC_ALPHAFOLD_ALLOW_ANY_HOST"
TRUSTED_HOST = "ftp.ebi.ac.uk"
SUMMARY_NAME = "index-v6.json"
SUMMARY_SHA256_KEY = "index_gz_sha256"
SUPPORTED_SCHEMA_VERSIONS = frozenset({"1"})
FORMATS = ("cif", "pdb")

_GZIP_MAGIC = b"\x1f\x8b"
_SQLITE_MAGIC = b"SQLite format 3\x00"
_CHUNK = 1 << 20
_INDEX_TIMEOUT = 120.0
_SUMMARY_TIMEOUT = 30.0
_PAGE = 5000
_STALE_PART_SECONDS = 24 * 3600.0
_REPLACE_ATTEMPTS = 5
_REPLACE_PAUSE = 0.2
_sleep = time.sleep

_ID_RE = re.compile(
    r"^(?:AF-)?(?P<acc>[A-Z0-9]+(?:-[0-9]+)?)"
    r"(?:-F(?P<frag>[0-9]+))?"
    r"(?:-MODEL_V(?P<version>[0-9]+)(?:\.(?P<ext>CIF|PDB)(?:\.GZ)?)?)?$"
)
_MEMBER_ENTRY_RE = re.compile(r"^AF-[A-Z0-9]+(?:-[0-9]+)?-F[0-9]+$")
_VERSION_RE = re.compile(r"^[0-9]+$")
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_ENTRY_RE = re.compile(r"^AF-(?P<acc>[A-Z0-9]+)(?:-(?P<iso>[0-9]+))?-F(?P<frag>[0-9]+)$")

_BULK_SET = (
    "The AlphaFold DB v6 bulk set is Swiss-Prot plus 46 proteome tars (16 model organisms "
    "and 30 global health organisms), not the whole database. This library reads only those "
    "tars; EBI's per-file endpoint is not used. A cached index never expires and EBI can add "
    "tars, so if yours is old, refresh_index() is worth trying."
)


@dataclass(frozen=True)
class Member:
    """Where one stored file lives: ``length`` gzip bytes at ``offset`` inside the tar at ``url``."""

    entry_id: str
    fmt: str
    tar: str
    url: str
    offset: int
    length: int


@dataclass(frozen=True)
class Tar:
    """One bulk tar, as recorded in the index when it was built."""

    name: str
    url: str
    size_bytes: int
    last_modified: str | None
    etag: str | None
    type: str | None
    reference_proteome: str | None
    tax_id: int | None
    species: str | None
    common_name: str | None
    n_members: int


@dataclass(frozen=True)
class ParsedId:
    """An id taken apart. ``version`` and ``fmt`` are set only when a file name was given."""

    accession: str
    fragment: int | None
    version: str | None
    fmt: str | None


def parse_id_full(text: str) -> ParsedId:
    """Take an id apart. The accession keeps an isoform suffix (``P12345-2``)."""
    if not isinstance(text, str):
        raise InvalidIdError(f"expected a UniProt accession or an entry id, got {type(text).__name__}")
    match = _ID_RE.match(text.strip().upper())
    if match is None:
        raise InvalidIdError(
            f"{text!r} is not a UniProt accession (P69905), an isoform accession (P12345-2) "
            "or an AlphaFold DB entry id (AF-P69905-F1)"
        )
    frag = match.group("frag")
    ext = match.group("ext")
    return ParsedId(
        match.group("acc"),
        int(frag) if frag is not None else None,
        match.group("version"),
        ext.lower() if ext is not None else None,
    )


def parse_id(text: str) -> tuple[str, int | None]:
    """Split an id into (accession, fragment number or None)."""
    parsed = parse_id_full(text)
    return parsed.accession, parsed.fragment


def check_fmt(fmt: str) -> str:
    value = fmt.strip().lower() if isinstance(fmt, str) else fmt
    if value not in FORMATS:
        raise ValueError(f"fmt must be 'cif' or 'pdb', got {fmt!r}")
    return value


def _is_url(source: str) -> bool:
    return source.startswith(("http://", "https://"))


def _entry_key(entry_id: str) -> tuple[int, int, str]:
    match = _ENTRY_RE.match(entry_id)
    if match is None:
        return (1 << 30, 0, entry_id)
    return (int(match.group("iso") or 0), int(match.group("frag")), entry_id)


def _abspath(text: str | os.PathLike[str]) -> Path:
    """Expand ``~`` and make the path absolute, so that it can be a file URI."""
    return Path(os.path.abspath(os.path.expanduser(os.fspath(text))))


def allow_any_host() -> bool:
    return os.environ.get(ENV_ALLOW_ANY_HOST, "").strip() == "1"


_PARTIAL_TAR = "UP000000429_85962_HELPY_v6.tar"


def _not_downloadable(url: str) -> str:
    if url == DEFAULT_INDEX_URL:
        return (
            f"The scigantic-alphafold index could not be downloaded ({url} returned 404); it may not "
            "be published yet or may be mid-update. Try again in a minute. To work without it, build "
            f"a partial index from one small tar (165 MB is streamed): `scigantic-alphafold build-part "
            f"{_PARTIAL_TAR} -o parts`, then `scigantic-alphafold merge parts -o index-v6.sqlite "
            f"--only {_PARTIAL_TAR}`, then set {ENV_INDEX}=index-v6.sqlite."
        )
    return f"No index at {url} (HTTP 404). Check {ENV_INDEX}, or point it at a local index-v6.sqlite file."


def _http_get(url: str, dest: Path | None, timeout: float) -> tuple[int, bytes]:
    """GET ``url`` into ``dest`` (or into memory). Returns (status, body); 404 is returned, not raised."""
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept-Encoding": "identity"})
    again = f"Try again, or point {ENV_INDEX} at a local index-v6.sqlite file."
    try:
        with urllib.request.urlopen(request, timeout=timeout) as response:
            if dest is None:
                return 200, bytes(response.read())
            with open(dest, "wb") as out:
                shutil.copyfileobj(response, out, _CHUNK)
            return 200, b""
    except urllib.error.HTTPError as exc:
        code = exc.code
        exc.close()
        if code == 404:
            return 404, b""
        raise IndexUnavailableError(f"Could not download {url}: HTTP {code}. {again}") from None
    except (OSError, http.client.HTTPException, ValueError) as exc:
        raise IndexUnavailableError(
            f"Could not download {url}: {str(exc) or type(exc).__name__}. {again}"
        ) from exc


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(_CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def _check_index(path: Path) -> None:
    """Open ``path`` read-only just long enough to check that this release can use it."""
    try:
        with open(path, "rb") as handle:
            magic = handle.read(len(_SQLITE_MAGIC))
    except (OSError, ValueError) as exc:
        raise IndexUnavailableError(f"Cannot read the index at {path}: {exc}") from exc
    if magic != _SQLITE_MAGIC:
        raise IndexUnavailableError(
            f"{path} is not a SQLite index (it does not start with the SQLite header). "
            "Call refresh_index() to download it again."
        )
    try:
        conn = sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)
        try:
            meta = {str(k): str(v) for k, v in conn.execute("SELECT key, value FROM meta")}
        finally:
            conn.close()
    except (sqlite3.Error, ValueError) as exc:
        raise IndexUnavailableError(
            f"{path} is not a scigantic-alphafold index ({exc}). Call refresh_index() to download it again."
        ) from exc
    version = meta.get("schema_version", "missing")
    if version not in SUPPORTED_SCHEMA_VERSIONS:
        raise IndexSchemaError(
            f"The index at {path} has schema_version {version}, and this release of scigantic-alphafold "
            f"understands only {', '.join(sorted(SUPPORTED_SCHEMA_VERSIONS))}. "
            "Upgrade the package (pip install -U scigantic-alphafold) or use an index built for this release."
        )
    afdb = meta.get("afdb_version")
    if afdb is None or not _VERSION_RE.match(afdb):
        raise IndexSchemaError(
            f"The index at {path} has afdb_version {afdb!r}; it must be a number such as '6'. "
            "The index is damaged or was not built by scigantic-alphafold."
        )


class Index:
    """One index, opened lazily and safe to share between threads.

    A single read-only SQLite connection is shared behind a lock. Lookups take
    microseconds, so the lock is never held across network I/O.
    """

    def __init__(self, source: str | None = None, cache_dir: str | os.PathLike[str] | None = None) -> None:
        self.source = source if source else DEFAULT_INDEX_URL
        # Both paths are made absolute here: a relative path cannot be a SQLite file URI,
        # and an unexpanded "~" would create a directory literally named "~".
        self._source_path: Path | None = None if _is_url(self.source) else _abspath(self.source)
        self.cache_dir = _abspath(cache_dir) if cache_dir else Path.home() / ".cache" / "scigantic-alphafold"
        self._lock = threading.RLock()
        self._conn: sqlite3.Connection | None = None
        self._path: Path | None = None
        self._tars: dict[str, Tar] | None = None
        self._meta: dict[str, str] | None = None

    # -- locating and materializing the file ---------------------------------------------

    def _local_source(self) -> Path:
        path = self._source_path
        assert path is not None
        if not path.is_file():
            shown = self.source if str(path) == self.source else f"{self.source} ({path})"
            raise IndexUnavailableError(
                f"{ENV_INDEX} points at {shown}, which is not a file. "
                "It must be a local index-v6.sqlite or index-v6.sqlite.gz, or an http(s) URL to one."
            )
        return path

    def _cache_name(self) -> str:
        if self.source == DEFAULT_INDEX_URL:
            return "index-v6.sqlite"
        if _is_url(self.source):
            return f"index-{hashlib.sha256(self.source.encode()).hexdigest()[:16]}.sqlite"
        local = self._local_source().resolve()
        stat = local.stat()
        key = f"{local}|{stat.st_size}|{stat.st_mtime_ns}"
        return f"index-local-{hashlib.sha256(key.encode()).hexdigest()[:16]}.sqlite"

    def _published_sha256(self) -> str | None:
        """The checksum the release publishes for the gzipped index, or None when there is none."""
        url = self.source.rsplit("/", 1)[0] + "/" + SUMMARY_NAME
        status, body = _http_get(url, None, _SUMMARY_TIMEOUT)
        if status == 404:
            return None
        try:
            summary = json.loads(body.decode("utf-8"))
        except ValueError as exc:
            raise IndexUnavailableError(
                f"{url} is not valid JSON ({exc}), so the index download cannot be verified. Try again in a minute."
            ) from None
        value = summary.get(SUMMARY_SHA256_KEY) if isinstance(summary, dict) else None
        if value is None:
            return None
        if not isinstance(value, str) or not _SHA256_RE.match(value.strip().lower()):
            raise IndexUnavailableError(
                f"{url} has a {SUMMARY_SHA256_KEY} that is not a SHA-256 ({value!r}), so the index "
                "download cannot be verified. Try again in a minute."
            )
        return value.strip().lower()

    def _sweep_stale_parts(self) -> None:
        """Remove temporary files a killed process left behind, once they are a day old."""
        try:
            cutoff = time.time() - _STALE_PART_SECONDS
            for entry in self.cache_dir.iterdir():
                if entry.name.startswith(".index-") and entry.name.endswith(".part"):
                    try:
                        if entry.stat().st_mtime < cutoff:
                            entry.unlink()
                    except OSError:
                        pass
        except OSError:
            pass

    def _replace(self, ready: Path, final: Path, force: bool) -> None:
        """Rename ``ready`` over ``final``, with nothing of ours holding ``final`` open.

        On Windows a file that any process has open cannot be replaced. Our own
        connection is closed first. If another process holds the file: a first open
        (not a refresh) simply uses the good copy that process put there; a refresh
        retries briefly and then says why it failed.
        """
        self.close()
        failure: OSError | None = None
        for attempt in range(_REPLACE_ATTEMPTS):
            if attempt:
                _sleep(_REPLACE_PAUSE)
            try:
                os.replace(ready, final)
                return
            except PermissionError as exc:
                failure = exc
                if not force and final.is_file():
                    try:
                        _check_index(final)
                        return
                    except AlphaFoldError:
                        pass
        raise IndexUnavailableError(
            f"Could not replace the cached index at {final}: {failure}. Another process probably has it "
            "open (on Windows an open file cannot be replaced). Close that process and try again."
        )

    def _materialize(self, final: Path, force: bool) -> None:
        """Fetch or copy the source, gunzip it if needed, check it, then rename it into place."""
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        temps: list[Path] = []

        def temp() -> Path:
            fd, name = tempfile.mkstemp(dir=self.cache_dir, prefix=".index-", suffix=".part")
            os.close(fd)
            temps.append(Path(name))
            return temps[-1]

        try:
            if _is_url(self.source):
                raw = temp()
                status, _ = _http_get(self.source, raw, _INDEX_TIMEOUT)
                if status == 404:
                    raise IndexUnavailableError(_not_downloadable(self.source))
                if self.source == DEFAULT_INDEX_URL:
                    expected = self._published_sha256()
                    actual = _file_sha256(raw)
                    if expected is not None and actual != expected:
                        raise IndexUnavailableError(
                            f"The index downloaded from {self.source} does not match the checksum published "
                            f"next to it in {SUMMARY_NAME} (expected sha256 {expected}, got {actual}). It was "
                            "not installed. The release may be mid-update: try again in a minute."
                        )
            else:
                raw = self._local_source()
            with open(raw, "rb") as handle:
                gzipped = handle.read(2) == _GZIP_MAGIC
            ready = raw
            if gzipped:
                ready = temp()
                try:
                    with gzip.open(raw, "rb") as source, open(ready, "wb") as out:
                        shutil.copyfileobj(source, out, _CHUNK)
                except (OSError, EOFError) as exc:
                    raise IndexUnavailableError(
                        f"The index from {self.source} is not a complete gzip file ({exc}). Try again."
                    ) from exc
            _check_index(ready)
            self._replace(ready, final, force)
        finally:
            for leftover in temps:
                try:
                    leftover.unlink()
                except OSError:
                    pass

    def _resolve_path(self, force: bool) -> Path:
        if not _is_url(self.source):
            local = self._local_source()
            with open(local, "rb") as handle:
                if handle.read(2) != _GZIP_MAGIC:
                    return local
        final = self.cache_dir / self._cache_name()
        self._sweep_stale_parts()
        if force or not final.is_file():
            self._materialize(final, force)
        return final

    # -- connection --------------------------------------------------------------------------

    def _open(self, force: bool = False) -> sqlite3.Connection:
        with self._lock:
            if self._conn is not None and not force:
                return self._conn
            try:
                path = self._resolve_path(force)
                _check_index(path)
                conn = sqlite3.connect(path.as_uri() + "?mode=ro", uri=True, check_same_thread=False)
            except AlphaFoldError:
                raise
            except (OSError, ValueError, sqlite3.Error) as exc:
                # A cache directory that cannot be created, a path with a null byte, a full disk:
                # whatever it is, the caller gets one of the library's own errors.
                raise IndexUnavailableError(
                    f"Cannot open the index from {self.source} (cache directory {self.cache_dir}): "
                    f"{str(exc) or type(exc).__name__}"
                ) from exc
            self.close()
            self._conn = conn
            self._path = path
            return conn

    def _query(self, sql: str, params: tuple[Any, ...] = ()) -> list[tuple[Any, ...]]:
        with self._lock:
            conn = self._open()
            try:
                return conn.execute(sql, params).fetchall()
            except sqlite3.Error as exc:
                raise IndexUnavailableError(
                    f"The index at {self._path} could not be read ({exc}). Call refresh_index()."
                ) from exc

    def close(self) -> None:
        with self._lock:
            if self._conn is not None:
                self._conn.close()
            self._conn = None
            self._tars = None
            self._meta = None

    @property
    def path(self) -> Path:
        with self._lock:
            self._open()
            assert self._path is not None
            return self._path

    def refresh(self) -> Path:
        """Fetch the index again and replace the cached copy. Returns the path now in use.

        The new copy is downloaded and checked first. Only then is this object's
        connection closed and the file replaced, so a failed download leaves the working
        index in use, and if the replace itself fails the old file is opened again on
        the next call.
        """
        with self._lock:
            self._open(force=True)
            assert self._path is not None
            return self._path

    # -- queries -----------------------------------------------------------------------------

    def meta(self) -> dict[str, str]:
        with self._lock:
            if self._meta is None:
                rows = self._query("SELECT key, value FROM meta")
                self._meta = {str(key): str(value) for key, value in rows}
            return dict(self._meta)

    def _tar_map(self) -> dict[str, Tar]:
        with self._lock:
            if self._tars is None:
                rows = self._query(
                    "SELECT name, url, size_bytes, last_modified, etag, type, reference_proteome, "
                    "tax_id, species, common_name, n_members FROM tars"
                )
                self._tars = {row[0]: Tar(*row) for row in sorted(rows)}
            return self._tars

    def tars(self) -> list[Tar]:
        return list(self._tar_map().values())

    def tar(self, name: str) -> Tar:
        try:
            return self._tar_map()[name]
        except KeyError:
            raise NotInBulkSetError(
                f"No tar named {name!r} in the index. tars() lists the {len(self._tar_map())} that are."
            ) from None

    def _wanted(self, id: str, fmt: str | None) -> tuple[ParsedId, str | None]:
        """Parse ``id`` and settle the format: the explicit ``fmt``, else the id's extension, else None."""
        parsed = parse_id_full(id)
        wanted = check_fmt(fmt) if fmt is not None else None
        if parsed.version is not None:
            indexed = self.meta().get("afdb_version", "")
            if not indexed.isdigit() or int(parsed.version) != int(indexed):
                raise InvalidIdError(
                    f"{id.strip()!r} names AlphaFold DB v{int(parsed.version)}, but the index describes "
                    f"v{indexed} and only that version can be read."
                )
        if parsed.fmt is not None:
            if wanted is not None and wanted != parsed.fmt:
                raise InvalidIdError(
                    f"{id.strip()!r} is a {parsed.fmt} file name but fmt={wanted!r} was asked for. "
                    "Pass the matching fmt, or leave fmt out and the extension decides."
                )
            wanted = parsed.fmt
        return parsed, wanted

    def locate(self, id: str, fmt: str | None = None) -> list[Member]:
        parsed, wanted = self._wanted(id, fmt)
        return self._rows(parsed, wanted)

    def _rows(self, parsed: ParsedId, wanted: str | None) -> list[Member]:
        accession, fragment = parsed.accession, parsed.fragment
        sql = (
            "SELECT m.entry_id, m.fmt, t.name, t.url, m.offset, m.length "
            "FROM members AS m JOIN tars AS t ON t.tar_id = m.tar_id WHERE "
        )
        params: tuple[Any, ...]
        if fragment is None:
            # Primary-key range scan: '-' is followed by '.' in byte order, so this is
            # exactly the rows whose entry_id starts with 'AF-<accession>-'.
            sql += "m.entry_id >= ? AND m.entry_id < ?"
            params = (f"AF-{accession}-", f"AF-{accession}.")
        else:
            sql += "m.entry_id = ?"
            params = (f"AF-{accession}-F{fragment}",)
        if wanted is not None:
            sql += " AND m.fmt = ?"
            params += (wanted,)
        members = [Member(*row) for row in self._query(sql, params)]
        for member in members:
            # The index is data from the network. These two fields end up in file names.
            if not isinstance(member.entry_id, str) or not _MEMBER_ENTRY_RE.match(member.entry_id):
                raise IndexSchemaError(f"The index holds a malformed entry id: {member.entry_id!r}")
            if member.fmt not in FORMATS:
                raise IndexSchemaError(f"The index holds an unknown format for {member.entry_id}: {member.fmt!r}")
        members.sort(key=lambda m: (_entry_key(str(m.entry_id)), str(m.fmt), str(m.tar).encode()))
        return members

    def resolve(self, id: str, fmt: str | None = None) -> Member:
        """The one member a fetch of ``id`` reads, or a clear error.

        ``fmt=None`` means the extension of a file-name id, and ``cif`` otherwise.
        """
        parsed, settled = self._wanted(id, fmt)
        accession, fragment = parsed.accession, parsed.fragment
        wanted = settled or "cif"
        members = self._rows(parsed, wanted)
        if not members:
            other = self._rows(parsed, None)
            if other:
                formats = ", ".join(sorted({str(m.fmt) for m in other}))
                raise NotInBulkSetError(
                    f"{id.strip()} is in the AlphaFold DB v6 bulk tars only as {formats}, not as {wanted}."
                )
            raise NotInBulkSetError(f"{id.strip()} is not in any AlphaFold DB v6 bulk tar. {_BULK_SET}")
        if fragment is None:
            own = re.compile(rf"^AF-{re.escape(accession)}-F[0-9]+$")
            candidates = list(dict.fromkeys(m.entry_id for m in members))
            canonical = [entry for entry in candidates if own.match(entry)]
            if len(canonical) != 1:
                if canonical:
                    why = f"{accession} is split into {len(canonical)} fragments"
                    candidates = canonical
                else:
                    why = f"{accession} has no entry of its own, only isoform entries"
                raise AmbiguousEntryError(
                    f"{why}; ask for one by its entry id: {', '.join(candidates)}", candidates
                )
            members = [m for m in members if m.entry_id == canonical[0]]
        self.check_member(members[0])
        return members[0]

    def check_member(self, member: Member) -> None:
        """Refuse a row that would send a request somewhere unexpected or outside its tar.

        The index is data from the network and it decides which URL and byte range are
        requested. (Entry ids and formats are checked as rows are read.)
        """
        size = self.tar(member.tar).size_bytes
        ok = isinstance(member.offset, int) and isinstance(member.length, int) and isinstance(size, int)
        if not ok or member.offset < 0 or member.length <= 0 or member.offset + member.length > size:
            raise IndexSchemaError(
                f"The index places {member.entry_id} at offset {member.offset!r}, length {member.length!r} "
                f"of {member.tar}, which is {size!r} bytes long. That range is not inside the tar."
            )
        self.check_tar_url(member.url)

    def check_tar_url(self, url: str) -> None:
        """Only fetch from where the index says its tars live, and by default only from EBI."""
        base = (self.meta().get("base_url") or "").strip()
        relaxed = allow_any_host()
        opt_out = f"For a private mirror or a test server, set {ENV_ALLOW_ANY_HOST}=1."
        try:
            base_parts = urllib.parse.urlsplit(base)
            base_host = base_parts.hostname
            scheme = urllib.parse.urlsplit(url).scheme if isinstance(url, str) else ""
        except ValueError:
            raise IndexSchemaError(f"The index holds a URL that cannot be parsed: {base!r} or {url!r}") from None
        if scheme not in ("http", "https"):
            raise IndexSchemaError(f"The index holds a tar URL that is not http(s): {url!r}")
        if not relaxed:
            if base_parts.scheme != "https" or base_host != TRUSTED_HOST:
                raise IndexSchemaError(
                    f"The index's base_url is {base!r}, which is not an https URL on {TRUSTED_HOST}, "
                    f"so nothing is fetched from it. {opt_out}"
                )
            if scheme != "https":
                raise IndexSchemaError(f"The index holds a tar URL that is not https: {url!r}. {opt_out}")
        if base or not relaxed:
            if not url.startswith(base.rstrip("/") + "/"):
                raise IndexSchemaError(
                    f"The index holds a tar URL outside its own base_url: {url!r} is not under {base!r}. "
                    "The index is damaged or was not built by scigantic-alphafold."
                )

    def has(self, id: str) -> bool:
        return bool(self.locate(id))

    def entries(self, *, proteome: str | None = None, tar: str | Tar | None = None, fmt: str = "cif") -> Iterator[str]:
        wanted = check_fmt(fmt)
        if proteome is not None and tar is not None:
            raise ValueError("pass proteome or tar, not both")
        names: list[str] | None = None
        if tar is not None:
            names = [self.tar(tar.name if isinstance(tar, Tar) else tar).name]
        elif proteome is not None:
            key = proteome.strip().upper()
            names = [t.name for t in self.tars() if (t.reference_proteome or "").upper() == key]
            if not names:
                raise NotInBulkSetError(
                    f"No bulk tar for proteome {proteome!r}. {_BULK_SET} tars() lists them."
                )
        return self._iter_entries(names, wanted)

    def _iter_entries(self, names: list[str] | None, fmt: str) -> Iterator[str]:
        sql = "SELECT DISTINCT entry_id FROM members WHERE fmt = ? AND entry_id > ?"
        extra: tuple[Any, ...] = ()
        if names is not None:
            marks = ", ".join("?" for _ in names)
            sql += f" AND tar_id IN (SELECT tar_id FROM tars WHERE name IN ({marks}))"
            extra = tuple(names)
        sql += " ORDER BY entry_id LIMIT ?"
        last = ""
        while True:
            rows = self._query(sql, (fmt, last, *extra, _PAGE))
            for (entry_id,) in rows:
                yield str(entry_id)
            if len(rows) < _PAGE:
                return
            last = str(rows[-1][0])

    def info(self) -> dict[str, Any]:
        meta = self.meta()
        path = self.path
        counts = self._query("SELECT COUNT(*), COALESCE(SUM(n_members), 0) FROM tars")[0]
        return {
            "path": str(path),
            "source": self.source,
            "size_bytes": path.stat().st_size,
            "schema_version": meta.get("schema_version"),
            "afdb_version": meta.get("afdb_version"),
            "built_at": meta.get("built_at"),
            "base_url": meta.get("base_url"),
            "tars": int(counts[0]),
            "members": int(counts[1]),
        }


_shared_lock = threading.Lock()
_shared: Index | None = None
_shared_key: tuple[str, str] | None = None


def get_index() -> Index:
    """The process-wide index for the current environment settings."""
    global _shared, _shared_key
    key = (os.environ.get(ENV_INDEX, "").strip(), os.environ.get(ENV_CACHE, "").strip())
    with _shared_lock:
        if _shared is None or _shared_key != key:
            if _shared is not None:
                _shared.close()
            _shared = Index(key[0] or None, key[1] or None)
            _shared_key = key
        return _shared


def reset() -> None:
    """Close and forget the process-wide index (the next call opens it again)."""
    global _shared, _shared_key
    with _shared_lock:
        if _shared is not None:
            _shared.close()
        _shared = None
        _shared_key = None


def locate(id: str, fmt: str | None = None) -> list[Member]:
    """Every stored copy of an accession or entry: all fragments, isoforms, formats and tars.

    Returns an empty list when nothing matches. The first row for a given entry and
    format is the copy that ``fetch`` reads. A file-name id (``AF-P69905-F1-model_v6.pdb.gz``)
    limits the rows to its own format unless ``fmt`` says the same thing; a conflicting
    ``fmt`` or another version number raises :class:`InvalidIdError`.
    """
    return get_index().locate(id, fmt)


def has(id: str) -> bool:
    """True when the bulk tars hold at least one file for this accession or entry id."""
    return get_index().has(id)


def entries(*, proteome: str | None = None, tar: str | Tar | None = None, fmt: str = "cif") -> Iterator[str]:
    """Iterate entry ids in sorted order: one proteome (``UP000005640``), one tar, or everything."""
    return get_index().entries(proteome=proteome, tar=tar, fmt=fmt)


def tars() -> list[Tar]:
    """The bulk tars the index covers, sorted by name."""
    return get_index().tars()


def index_info() -> dict[str, Any]:
    """Where the index is, what it was built from and how much it holds."""
    return get_index().info()


def refresh_index() -> Path:
    """Download the index again, replace the cached copy and return its path."""
    return get_index().refresh()
