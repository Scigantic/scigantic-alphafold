"""Locate, cache and query the offset index.

The index is one SQLite file that maps (entry id, format) to (tar, byte offset, length).
Nothing here runs at import time: the index is resolved, downloaded if needed and opened
on the first call that needs it.

Where the index comes from
--------------------------
1. ``SCIGANTIC_ALPHAFOLD_INDEX``, when set: a local path or an http(s) URL, to either a
   ``.sqlite`` file or a gzipped ``.sqlite.gz`` (told apart by content, not by name).
2. Otherwise the default release URL, :data:`DEFAULT_INDEX_URL`.

A plain local ``.sqlite`` is opened in place. Everything else is materialized once under
the cache directory (``SCIGANTIC_ALPHAFOLD_CACHE``, default ``~/.cache/scigantic-alphafold``):
written to a unique temporary name in that directory, gunzipped, checked, then renamed
over the final name with ``os.replace``. Two processes that race each do their own
download and the last rename wins; a reader never sees a half-written file.

How ids are resolved
--------------------
Accepted forms, in any letter case and with surrounding whitespace: a UniProt accession
(``P69905``), an isoform accession (``P12345-2``), an entry id (``AF-P69905-F1``,
``AF-P12345-2-F1``), the same without the ``AF-`` prefix (``P69905-F1``), and a member
file name (``AF-P69905-F1-model_v6.cif.gz``).

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
import os
import re
import shutil
import sqlite3
import tempfile
import threading
import urllib.error
import urllib.request
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ._errors import (
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
SUPPORTED_SCHEMA_VERSIONS = frozenset({"1"})
FORMATS = ("cif", "pdb")

_GZIP_MAGIC = b"\x1f\x8b"
_SQLITE_MAGIC = b"SQLite format 3\x00"
_CHUNK = 1 << 20
_INDEX_TIMEOUT = 120.0
_PAGE = 5000

_ID_RE = re.compile(
    r"^(?:AF-)?(?P<acc>[A-Z0-9]+(?:-[0-9]+)?)"
    r"(?:-F(?P<frag>[0-9]+))?"
    r"(?:-MODEL_V[0-9]+(?:\.(?:CIF|PDB))?(?:\.GZ)?)?$"
)
_ENTRY_RE = re.compile(r"^AF-(?P<acc>[A-Z0-9]+)(?:-(?P<iso>[0-9]+))?-F(?P<frag>[0-9]+)$")

_BULK_SET = (
    "The AlphaFold DB v6 bulk set is Swiss-Prot plus 46 proteome tars (16 model organisms "
    "and 30 global health organisms), not the whole database. This library reads only those "
    "tars; EBI's per-file endpoint is not used."
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


def parse_id(text: str) -> tuple[str, int | None]:
    """Split an id into (accession, fragment number or None). The accession keeps an isoform suffix."""
    if not isinstance(text, str):
        raise InvalidIdError(f"expected a UniProt accession or an entry id, got {type(text).__name__}")
    match = _ID_RE.match(text.strip().upper())
    if match is None:
        raise InvalidIdError(
            f"{text!r} is not a UniProt accession (P69905), an isoform accession (P12345-2) "
            "or an AlphaFold DB entry id (AF-P69905-F1)"
        )
    frag = match.group("frag")
    return match.group("acc"), (int(frag) if frag is not None else None)


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


def _not_published(url: str) -> str:
    if url == DEFAULT_INDEX_URL:
        return (
            f"The scigantic-alphafold index has not been published yet ({url} returned 404). "
            f"Build one with `scigantic-alphafold build-part <tar>` and `scigantic-alphafold merge "
            f"<parts dir> -o index-v6.sqlite`, then point {ENV_INDEX} at the local file, for example "
            f"{ENV_INDEX}=/path/to/index-v6.sqlite."
        )
    return f"No index at {url} (HTTP 404). Check {ENV_INDEX}, or point it at a local index-v6.sqlite file."


def _read_schema_version(path: Path) -> str:
    """Open ``path`` read-only just long enough to read and check ``meta.schema_version``."""
    try:
        with open(path, "rb") as handle:
            magic = handle.read(len(_SQLITE_MAGIC))
    except OSError as exc:
        raise IndexUnavailableError(f"Cannot read the index at {path}: {exc}") from exc
    if magic != _SQLITE_MAGIC:
        raise IndexUnavailableError(
            f"{path} is not a SQLite index (it does not start with the SQLite header). "
            "Call refresh_index() to download it again."
        )
    try:
        conn = sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)
        try:
            row = conn.execute("SELECT value FROM meta WHERE key = 'schema_version'").fetchone()
        finally:
            conn.close()
    except sqlite3.Error as exc:
        raise IndexUnavailableError(
            f"{path} is not a scigantic-alphafold index ({exc}). Call refresh_index() to download it again."
        ) from exc
    version = str(row[0]) if row is not None else "missing"
    if version not in SUPPORTED_SCHEMA_VERSIONS:
        raise IndexSchemaError(
            f"The index at {path} has schema_version {version}, and this release of scigantic-alphafold "
            f"understands only {', '.join(sorted(SUPPORTED_SCHEMA_VERSIONS))}. "
            "Upgrade the package (pip install -U scigantic-alphafold) or use an index built for this release."
        )
    return version


class Index:
    """One index, opened lazily and safe to share between threads.

    A single read-only SQLite connection is shared behind a lock. Lookups take
    microseconds, so the lock is never held across network I/O.
    """

    def __init__(self, source: str | None = None, cache_dir: str | os.PathLike[str] | None = None) -> None:
        self.source = source if source else DEFAULT_INDEX_URL
        self.cache_dir = Path(cache_dir) if cache_dir else Path.home() / ".cache" / "scigantic-alphafold"
        self._lock = threading.RLock()
        self._conn: sqlite3.Connection | None = None
        self._path: Path | None = None
        self._tars: dict[str, Tar] | None = None
        self._meta: dict[str, str] | None = None

    # -- locating and materializing the file ---------------------------------------------

    def _local_source(self) -> Path:
        path = Path(self.source).expanduser()
        if not path.is_file():
            raise IndexUnavailableError(
                f"{ENV_INDEX} points at {self.source}, which is not a file. "
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

    def _download(self, dest: Path) -> None:
        request = urllib.request.Request(
            self.source, headers={"User-Agent": USER_AGENT, "Accept-Encoding": "identity"}
        )
        try:
            with urllib.request.urlopen(request, timeout=_INDEX_TIMEOUT) as response, open(dest, "wb") as out:
                shutil.copyfileobj(response, out, _CHUNK)
        except urllib.error.HTTPError as exc:
            code = exc.code
            exc.close()
            if code == 404:
                raise IndexUnavailableError(_not_published(self.source)) from None
            raise IndexUnavailableError(
                f"Could not download the index from {self.source}: HTTP {code}. "
                f"Try again, or point {ENV_INDEX} at a local index-v6.sqlite file."
            ) from None
        except OSError as exc:
            raise IndexUnavailableError(
                f"Could not download the index from {self.source}: {exc}. "
                f"Try again, or point {ENV_INDEX} at a local index-v6.sqlite file."
            ) from exc

    def _materialize(self, final: Path) -> None:
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
                self._download(raw)
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
            _read_schema_version(ready)
            os.replace(ready, final)
        finally:
            for leftover in temps:
                try:
                    leftover.unlink()
                except FileNotFoundError:
                    pass

    def _resolve_path(self, force: bool) -> Path:
        if not _is_url(self.source):
            local = self._local_source()
            with open(local, "rb") as handle:
                if handle.read(2) != _GZIP_MAGIC:
                    return local
        final = self.cache_dir / self._cache_name()
        if force or not final.is_file():
            self._materialize(final)
        return final

    # -- connection --------------------------------------------------------------------------

    def _open(self, force: bool = False) -> sqlite3.Connection:
        with self._lock:
            if self._conn is not None and not force:
                return self._conn
            path = self._resolve_path(force)
            _read_schema_version(path)
            conn = sqlite3.connect(path.as_uri() + "?mode=ro", uri=True, check_same_thread=False)
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
        """Fetch the index again and replace the cached copy. Returns the path now in use."""
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

    def locate(self, id: str, fmt: str | None = None) -> list[Member]:
        accession, fragment = parse_id(id)
        wanted = check_fmt(fmt) if fmt is not None else None
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
        members.sort(key=lambda m: (_entry_key(m.entry_id), m.fmt, m.tar.encode()))
        return members

    def resolve(self, id: str, fmt: str = "cif") -> Member:
        """The one member a fetch of ``id`` reads, or a clear error."""
        accession, fragment = parse_id(id)
        wanted = check_fmt(fmt)
        members = self.locate(id, wanted)
        if not members:
            other = self.locate(id)
            if other:
                formats = ", ".join(sorted({m.fmt for m in other}))
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
        return members[0]

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
    format is the copy that ``fetch`` reads.
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
