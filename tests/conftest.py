"""Offline fixtures: synthetic bulk tars, a matching index, and a local HTTP server.

The server implements what the reader relies on (Range, Content-Range, Last-Modified,
ETag) and can be told to misbehave per path: ignore Range, fail a number of times, or
serve a huge body that must never be read.
"""

from __future__ import annotations

import gzip
import io
import os
import sqlite3
import tarfile
import threading
from collections.abc import Iterator
from dataclasses import dataclass, field
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

import pytest

from scigantic_alphafold import _fetch, _index

PROTEOME_TAR = "UP000000001_1_TEST_v6.tar"
SWISSPROT_TAR = "swissprot_cif_v6.tar"
LAST_MODIFIED = "Sun, 24 Aug 2025 10:52:01 GMT"

SCHEMA = """
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

# (entry id, format) per tar. Covers: a plain entry in both formats, a fragmented
# accession, an accession with an isoform, an accession with isoforms only, the same
# entry in two tars, and an entry stored as cif only.
PROTEOME_MEMBERS = [
    ("AF-P69905-F1", "cif"),
    ("AF-P69905-F1", "pdb"),
    ("AF-Q00001-F1", "cif"),
    ("AF-Q00001-F2", "cif"),
    ("AF-Q00001-F10", "cif"),
    ("AF-P12345-F1", "cif"),
    ("AF-P12345-2-F1", "cif"),
    ("AF-A0A000ISO1-2-F1", "cif"),
    ("AF-A0A000ISO1-3-F1", "cif"),
]
SWISSPROT_MEMBERS = [
    ("AF-P11111-F1", "cif"),
    ("AF-P69905-F1", "cif"),
]


def structure_text(entry_id: str, fmt: str, origin: str) -> str:
    if fmt == "cif":
        return f"data_{entry_id}\n#\n_entry.id {entry_id}\n# stored in {origin}\n" + "ATOM 1 N N . MET A 1 1 ? 0.0 0.0 0.0\n" * 40
    return f"HEADER    {entry_id}\nREMARK   stored in {origin}\n" + "ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00 90.00           N\n" * 40 + "END\n"


def build_tar(members: list[tuple[str, str]], origin: str) -> tuple[bytes, dict[tuple[str, str], tuple[int, int]]]:
    """A plain tar of individually gzipped members, plus each member's data offset and length."""
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w", format=tarfile.USTAR_FORMAT) as archive:
        for entry_id, fmt in members:
            payload = gzip.compress(structure_text(entry_id, fmt, origin).encode(), mtime=0)
            info = tarfile.TarInfo(f"{entry_id}-model_v6.{fmt}.gz")
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))
    data = buffer.getvalue()
    where: dict[tuple[str, str], tuple[int, int]] = {}
    with tarfile.open(fileobj=io.BytesIO(data)) as archive:
        for info in archive:
            entry_id, _, rest = info.name.partition("-model_v6.")
            where[(entry_id, rest.split(".")[0])] = (info.offset_data, info.size)
    return data, where


class Handler(BaseHTTPRequestHandler):
    protocol_version = "HTTP/1.1"
    server: "RangeServer"

    def log_message(self, format: str, *args: Any) -> None:
        pass

    def _reply(self, status: int, headers: dict[str, str], body: bytes = b"") -> None:
        self.send_response(status)
        for key, value in headers.items():
            self.send_header(key, value)
        self.send_header("Connection", "close")
        self.end_headers()
        self.close_connection = True
        if body:
            self.wfile.write(body)

    def do_GET(self) -> None:
        srv = self.server
        path = self.path
        range_header = self.headers.get("Range")
        with srv.lock:
            srv.requests.append({"path": path, "range": range_header, "user_agent": self.headers.get("User-Agent")})
            queued = srv.fail_first.get(path) or []
            forced = queued.pop(0) if queued else None
        if srv.delay:
            threading.Event().wait(srv.delay)
        if forced is not None:
            self._reply(forced, {"Content-Length": "0"})
            return
        if path in srv.huge:
            self._send_huge(srv.huge[path])
            return
        body = srv.files.get(path)
        if body is None:
            self._reply(404, {"Content-Length": "0"})
            return
        common = {"Accept-Ranges": "bytes", "Last-Modified": srv.last_modified.get(path, LAST_MODIFIED), "ETag": srv.etags.get(path, f'"{len(body):x}-test"')}
        if range_header and path not in srv.ignore_range:
            first_text, _, last_text = range_header.removeprefix("bytes=").partition("-")
            first, last = int(first_text), min(int(last_text), len(body) - 1)
            if first >= len(body):
                self._reply(416, {"Content-Range": f"bytes */{len(body)}", "Content-Length": "0"})
                return
            chunk = body[first : last + 1]
            if path in srv.truncate:
                # Promise the full range, deliver less, then hang up.
                self._reply(206, {**common, "Content-Range": f"bytes {first}-{last}/{len(body)}", "Content-Length": str(len(chunk))}, chunk[: len(chunk) // 2])
                with srv.lock:
                    srv.truncate.discard(path)
                return
            self._reply(206, {**common, "Content-Range": f"bytes {first}-{last}/{len(body)}", "Content-Length": str(len(chunk))}, chunk)
            return
        self._reply(200, {**common, "Content-Length": str(len(body))}, body)

    def _send_huge(self, size: int) -> None:
        """Answer 200 with ``size`` bytes, ignoring Range, and count how much the client took."""
        self._reply(200, {"Accept-Ranges": "bytes", "Last-Modified": LAST_MODIFIED, "Content-Length": str(size)})
        block = b"\0" * 65536
        sent = 0
        try:
            while sent < size:
                self.wfile.write(block)
                sent += len(block)
        except OSError:
            pass
        with self.server.lock:
            self.server.huge_sent.append(sent)
        self.server.huge_done.set()


class RangeServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(self) -> None:
        super().__init__(("127.0.0.1", 0), Handler)
        self.lock = threading.Lock()
        self.files: dict[str, bytes] = {}
        self.etags: dict[str, str] = {}
        self.last_modified: dict[str, str] = {}
        self.ignore_range: set[str] = set()
        self.truncate: set[str] = set()
        self.fail_first: dict[str, list[int]] = {}
        self.huge: dict[str, int] = {}
        self.huge_sent: list[int] = []
        self.huge_done = threading.Event()
        self.requests: list[dict[str, Any]] = []
        self.delay = 0.0

    @property
    def base(self) -> str:
        return f"http://127.0.0.1:{self.server_address[1]}"

    def count(self, path: str) -> int:
        with self.lock:
            return sum(1 for request in self.requests if request["path"] == path)


@dataclass
class Afdb:
    """Everything a test needs to know about the synthetic database."""

    server: RangeServer
    root: Path
    index_path: Path
    cache_dir: Path
    tars: dict[str, bytes]
    where: dict[str, dict[tuple[str, str], tuple[int, int]]]
    monkeypatch: pytest.MonkeyPatch
    sleeps: list[float] = field(default_factory=list)

    def tar_path(self, name: str) -> str:
        return f"/v6/{name}"

    def tar_url(self, name: str) -> str:
        return self.server.base + self.tar_path(name)

    def text(self, entry_id: str, fmt: str = "cif", tar: str = PROTEOME_TAR) -> str:
        return structure_text(entry_id, fmt, tar)

    def use(self, source: str | os.PathLike[str] | None) -> None:
        """Point the library at another index source (None means the default URL)."""
        if source is None:
            self.monkeypatch.delenv(_index.ENV_INDEX, raising=False)
        else:
            self.monkeypatch.setenv(_index.ENV_INDEX, str(source))
        _index.reset()

    def write_index(self, path: Path, *, schema_version: str = "1", etags: dict[str, str] | None = None) -> Path:
        if path.exists():
            path.unlink()
        conn = sqlite3.connect(path)
        conn.executescript(SCHEMA)
        conn.executemany(
            "INSERT INTO meta VALUES (?, ?)",
            [("schema_version", schema_version), ("afdb_version", "6"), ("built_at", "2026-11-20T00:00:00Z"), ("base_url", self.server.base + "/v6/")],
        )
        described = {
            # Inserted in this order on purpose: tar_id order must not decide which copy wins.
            SWISSPROT_TAR: ("swissprot", None, None, None, None),
            PROTEOME_TAR: ("proteome", "UP000000001", 1, "Testus exemplaris", "test organism"),
        }
        for tar_id, (name, extra) in enumerate(described.items(), start=1):
            data = self.tars[name]
            etag = (etags or {}).get(name, f'"{len(data):x}-test"')
            conn.execute(
                "INSERT INTO tars VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (tar_id, name, self.tar_url(name), len(data), LAST_MODIFIED, etag, *extra, len(self.where[name])),
            )
            conn.executemany(
                "INSERT INTO members VALUES (?, ?, ?, ?, ?)",
                [(entry_id, fmt, tar_id, offset, length) for (entry_id, fmt), (offset, length) in self.where[name].items()],
            )
        conn.commit()
        conn.close()
        return path


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    if os.environ.get("SCIGANTIC_ALPHAFOLD_LIVE") == "1":
        return
    skip = pytest.mark.skip(reason="live test: set SCIGANTIC_ALPHAFOLD_LIVE=1 to run it")
    for item in items:
        if "live" in item.keywords:
            item.add_marker(skip)


@pytest.fixture
def afdb(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Iterator[Afdb]:
    server = RangeServer()
    thread = threading.Thread(target=server.serve_forever, kwargs={"poll_interval": 0.01}, daemon=True)
    thread.start()

    tars: dict[str, bytes] = {}
    where: dict[str, dict[tuple[str, str], tuple[int, int]]] = {}
    for name, members in ((PROTEOME_TAR, PROTEOME_MEMBERS), (SWISSPROT_TAR, SWISSPROT_MEMBERS)):
        tars[name], where[name] = build_tar(members, name)
        server.files[f"/v6/{name}"] = tars[name]

    db = Afdb(server, tmp_path, tmp_path / "index-v6.sqlite", tmp_path / "cache", tars, where, monkeypatch)
    db.write_index(db.index_path)
    plain = db.index_path.read_bytes()
    (tmp_path / "index-v6.sqlite.gz").write_bytes(gzip.compress(plain))
    server.files["/index/index-v6.sqlite"] = plain
    server.files["/index/index-v6.sqlite.gz"] = gzip.compress(plain)

    monkeypatch.setenv(_index.ENV_CACHE, str(db.cache_dir))
    monkeypatch.setenv(_index.ENV_INDEX, str(db.index_path))
    monkeypatch.setattr(_fetch, "_sleep", db.sleeps.append)
    _index.reset()
    try:
        yield db
    finally:
        _index.reset()
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)
