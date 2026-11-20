"""Tests for the index builder.

Everything here is offline except the tests marked ``live``: a synthetic tar
of gzipped members is served from a local HTTP server that supports Range and
can be told to drop connections, answer 503 or change its ETag.
"""

from __future__ import annotations

import base64
import gzip
import http.server
import io
import json
import os
import random
import re
import socket
import sqlite3
import sys
import tarfile
import threading
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Iterator

import pytest

try:
    from scigantic_alphafold import _build
except ImportError:  # the package is not installed: run from the source tree
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from scigantic_alphafold import _build

live = pytest.mark.skipif(
    os.environ.get("SCIGANTIC_ALPHAFOLD_LIVE") != "1",
    reason="touches the real EBI server; set SCIGANTIC_ALPHAFOLD_LIVE=1 to run",
)

LAST_MODIFIED = "Sun, 24 Aug 2025 10:52:01 GMT"


# --------------------------------------------------------------------------
# Synthetic data


def _filler(rng: random.Random, n: int) -> str:
    return base64.b64encode(rng.randbytes(n)).decode()


def cif_text(entry_id: str, rng: random.Random) -> str:
    return (
        f"data_{entry_id}\n#\n_entry.id {entry_id}\n#\nloop_\n_atom_site.group_PDB\n"
        f"_atom_site.id\nATOM 1\n# {_filler(rng, rng.randrange(4000, 30000))}\n"
    )


def pdb_text(entry_id: str, rng: random.Random) -> str:
    accession = entry_id.split("-")[1]
    return (
        "HEADER                                            01-AUG-25\n"
        f"TITLE     ALPHAFOLD MONOMER V2.0 PREDICTION FOR A PROTEIN ({accession})\n"
        f"DBREF  XXXX A    1   217  UNP    {accession}   {accession}_TEST     1    217\n"
        f"REMARK   1 {_filler(rng, rng.randrange(3000, 20000))}\n"
        "ATOM      1  N   MET A   1     -12.230  20.373  10.905  1.00 88.69           N\n"
        "END\n"
    )


class Tar:
    def __init__(self, data: bytes, originals: dict[tuple[str, str], bytes], skipped: list[str]):
        self.data = data
        self.originals = originals  # (entry_id, fmt) -> the .gz bytes stored in the tar
        self.skipped = skipped


def make_tar(
    entry_ids: list[str],
    *,
    formats: tuple[str, ...] = ("cif", "pdb"),
    extras: bool = True,
    seed: int = 0,
    tar_format: int = tarfile.GNU_FORMAT,
) -> Tar:
    rng = random.Random(seed)
    originals: dict[tuple[str, str], bytes] = {}
    skipped: list[str] = []
    buffer = io.BytesIO()

    def add(archive: tarfile.TarFile, name: str, payload: bytes) -> None:
        info = tarfile.TarInfo(name)
        info.size = len(payload)
        info.mtime = 1755966064
        archive.addfile(info, io.BytesIO(payload))

    with tarfile.open(fileobj=buffer, mode="w", format=tar_format) as archive:
        for position, entry_id in enumerate(entry_ids):
            for fmt in formats:
                text = cif_text(entry_id, rng) if fmt == "cif" else pdb_text(entry_id, rng)
                blob = gzip.compress(text.encode(), mtime=0)
                originals[(entry_id, fmt)] = blob
                add(archive, f"{entry_id}-model_v6.{fmt}.gz", blob)
            if extras and position == 1:
                # Things that must be listed as skipped, never indexed.
                names = [
                    "README.txt",
                    "AF-P12345-F1-model_v4.cif.gz",
                    "AF-P12345-F1-predicted_aligned_error_v6.json.gz",
                    "nested/AF-P12345-F1-model_v6.cif.gz",
                ]
                if tar_format != tarfile.USTAR_FORMAT:
                    names.append("long-" + "x" * 150 + ".txt")
                for name in names:
                    add(archive, name, b"not a model " * 40)
                    skipped.append(name)
                directory = tarfile.TarInfo("extra")
                directory.type = tarfile.DIRTYPE
                archive.addfile(directory)
                skipped.append("extra")
    return Tar(buffer.getvalue(), originals, skipped)


def entry_ids(prefix: str, n: int) -> list[str]:
    ids = [f"AF-{prefix}{i:04d}-F1" for i in range(n)]
    ids[-1] = f"AF-{prefix}0000-2-F1"  # an isoform-style id
    return ids


# --------------------------------------------------------------------------
# Local HTTP server with Range, drops, 503s and a mutable ETag


class Site:
    def __init__(self) -> None:
        self.files: dict[str, bytes] = {}
        self.etags: dict[str, str] = {}
        self.drops: list[int] = []  # for each coming tar GET: bytes to send before hanging up
        self.fail_5xx = 0  # number of coming GETs to answer with 503
        self.change_etag_after: int | None = None  # change every ETag after this many GETs
        self.requests: list[tuple[str, str, str | None, str | None]] = []
        self.base_url = ""

    def add(self, name: str, data: bytes) -> None:
        self.files[name] = data
        self.etags[name] = f'"{len(data):x}-{abs(hash(name)) % 10**8:x}"'

    def gets(self, name: str) -> list[tuple[str, str, str | None, str | None]]:
        return [r for r in self.requests if r[0] == "GET" and r[1] == name]


def _handler(site: Site) -> type[http.server.BaseHTTPRequestHandler]:
    class Handler(http.server.BaseHTTPRequestHandler):
        def log_message(self, format: str, *args: Any) -> None:
            pass

        def do_HEAD(self) -> None:
            self._serve(send_body=False)

        def do_GET(self) -> None:
            self._serve(send_body=True)

        def _serve(self, *, send_body: bool) -> None:
            name = self.path.rsplit("/", 1)[-1]
            rng = self.headers.get("Range")
            site.requests.append((self.command, name, rng, self.headers.get("If-Range")))
            assert self.headers.get("User-Agent", "").startswith("scigantic-alphafold/")
            if name not in site.files:
                self.send_error(404)
                return
            if send_body and site.fail_5xx > 0:
                site.fail_5xx -= 1
                self.send_error(503)
                return
            if site.change_etag_after is not None and send_body:
                if len([r for r in site.requests if r[0] == "GET"]) > site.change_etag_after:
                    site.etags = {k: '"changed"' for k in site.etags}
            body = site.files[name]
            total = len(body)
            start, end = 0, total - 1
            status = 200
            if rng:
                match = re.fullmatch(r"bytes=(\d+)-(\d*)", rng)
                assert match, rng
                start = int(match.group(1))
                end = int(match.group(2)) if match.group(2) else total - 1
                if start >= total:
                    self.send_error(416)
                    return
                end = min(end, total - 1)
                status = 206
            payload = body[start : end + 1]
            self.send_response(status)
            self.send_header("Content-Type", "application/x-tar")
            self.send_header("Content-Length", str(len(payload)))
            self.send_header("Accept-Ranges", "bytes")
            self.send_header("ETag", site.etags[name])
            self.send_header("Last-Modified", LAST_MODIFIED)
            if status == 206:
                self.send_header("Content-Range", f"bytes {start}-{end}/{total}")
            self.end_headers()
            if not send_body:
                return
            if site.drops and name.endswith(".tar"):
                keep = site.drops.pop(0)
                self.wfile.write(payload[:keep])
                self.wfile.flush()
                try:
                    self.connection.shutdown(socket.SHUT_RDWR)
                except OSError:
                    pass
                self.close_connection = True
                return
            try:
                self.wfile.write(payload)
            except (BrokenPipeError, ConnectionResetError):
                pass

    return Handler


@pytest.fixture()
def site() -> Iterator[Site]:
    state = Site()
    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), _handler(state))
    server.daemon_threads = True
    thread = threading.Thread(target=server.serve_forever, args=(0.01,), daemon=True)
    thread.start()
    state.base_url = f"http://127.0.0.1:{server.server_address[1]}/v6"
    try:
        yield state
    finally:
        server.shutdown()
        server.server_close()
        thread.join(timeout=5)


@pytest.fixture(autouse=True)
def no_sleep(monkeypatch: pytest.MonkeyPatch) -> list[float]:
    delays: list[float] = []
    monkeypatch.setattr(_build, "_sleep", delays.append)
    return delays


def read_part(path: Path) -> tuple[dict[str, Any], list[tuple[str, str, int, int]], list[str], Any]:
    lines = gzip.decompress(path.read_bytes()).decode().splitlines()
    assert lines[0].startswith("#{")
    header = json.loads(lines[0][1:])
    rows: list[tuple[str, str, int, int]] = []
    skipped: list[str] = []
    trailer = None
    for line in lines[1:]:
        if line.startswith("#skipped\t"):
            skipped.append(line.split("\t", 1)[1])
        elif line.startswith("#end\t"):
            trailer = json.loads(line.split("\t", 1)[1])
        else:
            entry_id, fmt, offset, length = line.split("\t")
            rows.append((entry_id, fmt, int(offset), int(length)))
    return header, rows, skipped, trailer


def assert_rows_match(tar: Tar, rows: list[tuple[str, str, int, int]]) -> None:
    """Every row's byte range is exactly the stored .gz and gunzips to a model."""
    assert {(r[0], r[1]) for r in rows} == set(tar.originals)
    assert len(rows) == len(tar.originals)
    for entry_id, fmt, offset, length in rows:
        blob = tar.data[offset : offset + length]
        assert blob == tar.originals[(entry_id, fmt)]
        text = gzip.decompress(blob).decode()
        assert _build._check_model_text(entry_id, fmt, text) is None


NAME = "UP000000001_111_TESTA_v6.tar"


# --------------------------------------------------------------------------
# build_part


@pytest.mark.parametrize(
    "tar_format", [tarfile.GNU_FORMAT, tarfile.PAX_FORMAT, tarfile.USTAR_FORMAT]
)
def test_build_part_is_exact(site: Site, tmp_path: Path, tar_format: int) -> None:
    tar = make_tar(entry_ids("P", 30), tar_format=tar_format)
    assert len(tar.data) > 3 * _build._BUFSIZE  # several stream reads
    site.add(NAME, tar.data)

    part = _build.build_part(NAME, tmp_path / "parts", base_url=site.base_url)

    assert part == tmp_path / "parts" / (NAME + ".tsv.gz")
    header, rows, skipped, trailer = read_part(part)
    assert header == {
        "name": NAME,
        "url": f"{site.base_url}/{NAME}",
        "size_bytes": len(tar.data),
        "last_modified": LAST_MODIFIED,
        "etag": site.etags[NAME],
    }
    assert_rows_match(tar, rows)
    assert skipped == tar.skipped
    assert trailer["n_cif"] == 30 and trailer["n_pdb"] == 30
    assert trailer["n_skipped"] == len(tar.skipped)
    assert trailer["bytes_read"] == len(tar.data)
    assert trailer["reconnects"] == 0
    # Rows are in tar order and offsets sit on 512-byte block boundaries.
    assert [r[2] for r in rows] == sorted(r[2] for r in rows)
    assert all(r[2] % 512 == 0 for r in rows)
    assert ("AF-P0000-2-F1", "cif") in {(r[0], r[1]) for r in rows}
    # One plain GET, no Range, nothing left behind but the part.
    assert site.gets(NAME) == [("GET", NAME, None, None)]
    assert [p.name for p in (tmp_path / "parts").iterdir()] == [NAME + ".tsv.gz"]


def test_build_part_resumes_after_dropped_connections(site: Site, tmp_path: Path) -> None:
    tar = make_tar(entry_ids("Q", 30))
    site.add(NAME, tar.data)
    clean = read_part(_build.build_part(NAME, tmp_path / "clean", base_url=site.base_url))
    site.requests.clear()

    # Drop after 100000 bytes, then after 1 more byte, then mid-way again,
    # then before sending anything at all.
    site.drops = [100_000, 1, 300_001, 0]
    part = _build.build_part(NAME, tmp_path / "dropped", base_url=site.base_url)

    header, rows, skipped, trailer = read_part(part)
    assert (header, rows, skipped) == clean[:3]
    assert_rows_match(tar, rows)
    assert trailer["reconnects"] == 4
    assert trailer["bytes_read"] == len(tar.data)
    etag = site.etags[NAME]
    assert site.gets(NAME) == [
        ("GET", NAME, None, None),
        ("GET", NAME, "bytes=100000-", etag),
        ("GET", NAME, "bytes=100001-", etag),
        ("GET", NAME, "bytes=400002-", etag),
        ("GET", NAME, "bytes=400002-", etag),
    ]


def test_build_part_retries_5xx_with_backoff(
    site: Site, tmp_path: Path, no_sleep: list[float]
) -> None:
    tar = make_tar(entry_ids("R", 6))
    site.add(NAME, tar.data)
    site.fail_5xx = 3
    part = _build.build_part(NAME, tmp_path, base_url=site.base_url, backoff=2.0)
    assert_rows_match(tar, read_part(part)[1])
    assert no_sleep == [2.0, 4.0, 8.0]
    assert len(site.gets(NAME)) == 4


def test_build_part_gives_up_after_bounded_retries(site: Site, tmp_path: Path) -> None:
    site.add(NAME, make_tar(entry_ids("R", 6)).data)
    site.fail_5xx = 100
    with pytest.raises(_build.BuildError, match="gave up"):
        _build.build_part(NAME, tmp_path, base_url=site.base_url, max_retries=2)
    assert len(site.gets(NAME)) == 3  # the first attempt and two retries
    assert list(tmp_path.iterdir()) == []


def test_build_part_gives_up_when_every_resume_is_dropped(site: Site, tmp_path: Path) -> None:
    site.add(NAME, make_tar(entry_ids("R", 6)).data)
    site.drops = [1000] + [0] * 50
    with pytest.raises(_build.BuildError, match="gave up at byte 1000"):
        _build.build_part(NAME, tmp_path, base_url=site.base_url, max_retries=3)
    # The drop itself is failure 1; three resumes that deliver nothing follow.
    assert [r[2] for r in site.gets(NAME)] == [None] + ["bytes=1000-"] * 3
    assert list(tmp_path.iterdir()) == []


def test_build_part_does_not_retry_4xx(site: Site, tmp_path: Path) -> None:
    with pytest.raises(urllib.error.HTTPError) as caught:
        _build.build_part(NAME, tmp_path, base_url=site.base_url)
    assert caught.value.code == 404
    assert len(site.gets(NAME)) == 1
    assert list(tmp_path.iterdir()) == []


def test_build_part_stops_if_the_tar_changes_mid_build(site: Site, tmp_path: Path) -> None:
    site.add(NAME, make_tar(entry_ids("S", 30)).data)
    site.drops = [100_000]
    site.change_etag_after = 1
    with pytest.raises(_build.TarChangedError, match="changed during the build"):
        _build.build_part(NAME, tmp_path, base_url=site.base_url)
    assert list(tmp_path.iterdir()) == []


def test_build_part_rejects_a_truncated_tar(site: Site, tmp_path: Path) -> None:
    tar = make_tar(entry_ids("T", 10))
    site.add(NAME, tar.data[: len(tar.data) // 2 + 100])
    with pytest.raises(_build.BuildError):
        _build.build_part(NAME, tmp_path, base_url=site.base_url)
    assert list(tmp_path.iterdir()) == []


def test_build_part_rejects_a_corrupt_header_mid_archive(site: Site, tmp_path: Path) -> None:
    tar = make_tar(entry_ids("T", 10), extras=False)
    clean = make_tar(entry_ids("T", 10), extras=False)
    site.add(NAME, clean.data)
    rows = read_part(_build.build_part(NAME, tmp_path / "a", base_url=site.base_url))[1]
    header_at = rows[7][2] - 512
    broken = bytearray(tar.data)
    broken[header_at : header_at + 512] = b"\xff" * 512
    site.add(NAME, bytes(broken))
    with pytest.raises(_build.BuildError, match="parsing stopped|could not be parsed"):
        _build.build_part(NAME, tmp_path / "b", base_url=site.base_url)
    assert not (tmp_path / "b" / (NAME + ".tsv.gz")).exists()


def test_build_part_rejects_data_after_the_end_marker(site: Site, tmp_path: Path) -> None:
    first = make_tar(entry_ids("T", 4), extras=False, seed=1)
    second = make_tar(entry_ids("U", 4), extras=False, seed=2)
    site.add(NAME, first.data + second.data)
    with pytest.raises(_build.BuildError, match="parsing stopped at byte"):
        _build.build_part(NAME, tmp_path, base_url=site.base_url)
    assert list(tmp_path.iterdir()) == []


def test_build_part_replaces_an_older_part_only_on_success(site: Site, tmp_path: Path) -> None:
    tar = make_tar(entry_ids("V", 5))
    site.add(NAME, tar.data)
    part = _build.build_part(NAME, tmp_path, base_url=site.base_url)
    before = part.read_bytes()
    site.fail_5xx = 100
    with pytest.raises(_build.BuildError):
        _build.build_part(NAME, tmp_path, base_url=site.base_url, max_retries=1)
    assert part.read_bytes() == before
    assert [p.name for p in tmp_path.iterdir()] == [part.name]


@pytest.mark.parametrize("bad", ["", "../x.tar", "a/b.tar", "notatar.zip", "http://h/x.tar"])
def test_build_part_rejects_bad_tar_names(tmp_path: Path, bad: str) -> None:
    with pytest.raises(ValueError):
        _build.build_part(bad, tmp_path, base_url="http://127.0.0.1:1")


# --------------------------------------------------------------------------
# merge_parts


class Corpus:
    def __init__(self, site: Site, parts: Path) -> None:
        self.site = site
        self.parts = parts
        self.tars: dict[str, Tar] = {}
        self.metadata: list[dict[str, Any]] = []

    def add(
        self,
        name: str,
        tar: Tar,
        *,
        n: int,
        type_: str,
        build: bool = True,
        **fields: Any,
    ) -> None:
        self.tars[name] = tar
        self.site.add(name, tar.data)
        self.metadata.append(
            {"archive_name": name, "num_predicted_structures": n, "size_bytes": len(tar.data),
             "type": type_, **fields}
        )
        if build:
            _build.build_part(name, self.parts, base_url=self.site.base_url)


A = "UP000000001_111_TESTA_v6.tar"
B = "UP000000002_2222_9TEST1_v6.tar"
SC = "swissprot_cif_v6.tar"
SP = "swissprot_pdb_v6.tar"


@pytest.fixture()
def corpus(site: Site, tmp_path: Path) -> Corpus:
    c = Corpus(site, tmp_path / "parts")
    c.add(A, make_tar(entry_ids("A", 12), seed=1), n=12, type_="proteome",
          species="Testus alpha", common_name="Alpha", reference_proteome="UP000000001")
    c.add(B, make_tar(entry_ids("B", 5), seed=2, extras=False), n=5, type_="global_health",
          species="Testus beta", common_name=None, reference_proteome="UP000000002")
    # Swiss-Prot repeats entries that are also in proteome tars.
    shared = entry_ids("A", 12)[:3] + [f"AF-S{i:04d}-F1" for i in range(6)]
    c.add(SC, make_tar(shared, formats=("cif",), seed=3, extras=False), n=9, type_="swissprot",
          species=None, common_name=None, reference_proteome=None, label="Swiss-Prot (CIF files)")
    c.add(SP, make_tar(shared, formats=("pdb",), seed=4, extras=False), n=9, type_="swissprot",
          species=None, common_name=None, reference_proteome=None, label="Swiss-Prot (PDB files)")
    return c


CONTRACT_SCHEMA = {
    "meta": "CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT NOT NULL)",
    "tars": (
        "CREATE TABLE tars( tar_id INTEGER PRIMARY KEY, name TEXT NOT NULL UNIQUE, "
        "url TEXT NOT NULL, size_bytes INTEGER NOT NULL, last_modified TEXT, etag TEXT, "
        "type TEXT, reference_proteome TEXT, tax_id INTEGER, species TEXT, common_name TEXT, "
        "n_members INTEGER NOT NULL )"
    ),
    "members": (
        "CREATE TABLE members( entry_id TEXT NOT NULL, fmt TEXT NOT NULL, "
        "tar_id INTEGER NOT NULL, offset INTEGER NOT NULL, length INTEGER NOT NULL, "
        "PRIMARY KEY(entry_id, fmt, tar_id) ) WITHOUT ROWID"
    ),
}


def test_merge_builds_exactly_the_contract_index(
    corpus: Corpus, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    out = tmp_path / "out" / "index-v6.sqlite"
    assert _build.merge_parts(corpus.parts, out, metadata=corpus.metadata) == out
    printed = capsys.readouterr().out

    conn = sqlite3.connect(out)
    schema = {
        name: " ".join(sql.split())
        for name, sql in conn.execute("SELECT name, sql FROM sqlite_master WHERE type='table'")
        if not name.startswith("sqlite_")
    }
    assert schema == CONTRACT_SCHEMA
    assert conn.execute("PRAGMA integrity_check").fetchone() == ("ok",)

    meta = dict(conn.execute("SELECT key, value FROM meta"))
    assert set(meta) == {"schema_version", "afdb_version", "built_at", "base_url"}
    assert meta["schema_version"] == "1" and meta["afdb_version"] == "6"
    assert meta["base_url"] == corpus.site.base_url
    assert re.fullmatch(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ", meta["built_at"])

    base, etags = corpus.site.base_url, corpus.site.etags
    assert conn.execute("SELECT * FROM tars ORDER BY tar_id").fetchall() == [
        (1, A, f"{base}/{A}", len(corpus.tars[A].data), LAST_MODIFIED, etags[A],
         "proteome", "UP000000001", 111, "Testus alpha", "Alpha", 24),
        (2, B, f"{base}/{B}", len(corpus.tars[B].data), LAST_MODIFIED, etags[B],
         "global_health", "UP000000002", 2222, "Testus beta", None, 10),
        (3, SC, f"{base}/{SC}", len(corpus.tars[SC].data), LAST_MODIFIED, etags[SC],
         "swissprot", None, None, None, None, 9),
        (4, SP, f"{base}/{SP}", len(corpus.tars[SP].data), LAST_MODIFIED, etags[SP],
         "swissprot", None, None, None, None, 9),
    ]

    names = {1: A, 2: B, 3: SC, 4: SP}
    members = conn.execute("SELECT entry_id, fmt, tar_id, offset, length FROM members").fetchall()
    assert len(members) == 24 + 10 + 9 + 9
    expected = {(e, f, t) for t, name in names.items() for (e, f) in corpus.tars[name].originals}
    assert {(e, f, t) for e, f, t, _, _ in members} == expected
    for entry_id, fmt, tar_id, offset, length in members:
        tar = corpus.tars[names[tar_id]]
        blob = tar.data[offset : offset + length]
        assert blob == tar.originals[(entry_id, fmt)]
        assert gzip.decompress(blob).decode().startswith(("data_" + entry_id, "HEADER"))

    # The accession lookup from the contract: a primary-key range scan that
    # also finds the isoform and every tar holding the entry.
    found = conn.execute(
        "SELECT entry_id, fmt, tar_id FROM members WHERE entry_id >= ? AND entry_id < ? "
        "ORDER BY entry_id, fmt, tar_id",
        ("AF-A0000-", "AF-A0000."),
    ).fetchall()
    assert found == [
        ("AF-A0000-2-F1", "cif", 1), ("AF-A0000-2-F1", "pdb", 1),
        ("AF-A0000-F1", "cif", 1), ("AF-A0000-F1", "cif", 3),
        ("AF-A0000-F1", "pdb", 1), ("AF-A0000-F1", "pdb", 4),
    ]
    conn.close()

    gz = out.with_name("index-v6.sqlite.gz")
    assert gzip.decompress(gz.read_bytes()) == out.read_bytes()
    assert sorted(p.name for p in out.parent.iterdir()) == ["index-v6.sqlite", "index-v6.sqlite.gz"]
    for name in names.values():
        assert name in printed
    assert "12/12" in printed and "9/0" in printed and "0/9" in printed


def test_merge_fetches_the_tar_list_when_not_given(
    corpus: Corpus, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    listing = corpus.metadata + [
        {"archive_name": "UP000000001_111_TESTA_v4.tar", "num_predicted_structures": 1,
         "size_bytes": 1, "type": "proteome"}
    ]
    corpus.site.files["download_metadata.json"] = json.dumps(listing).encode()
    corpus.site.etags["download_metadata.json"] = '"m"'
    monkeypatch.setattr(
        _build, "METADATA_URL", corpus.site.base_url.rsplit("/", 1)[0] + "/download_metadata.json"
    )
    out = _build.merge_parts(corpus.parts, tmp_path / "index-v6.sqlite")
    conn = sqlite3.connect(out)
    assert conn.execute("SELECT count(*) FROM tars").fetchone() == (4,)
    conn.close()


def _expect_merge_failure(corpus: Corpus, tmp_path: Path, pattern: str) -> str:
    out = tmp_path / "out" / "index-v6.sqlite"
    out.parent.mkdir(exist_ok=True)
    out.write_bytes(b"previous index")
    with pytest.raises(_build.MergeError, match=pattern) as caught:
        _build.merge_parts(corpus.parts, out, metadata=corpus.metadata)
    # Nothing was written: the previous index is untouched, no temp files.
    assert out.read_bytes() == b"previous index"
    assert [p.name for p in out.parent.iterdir()] == ["index-v6.sqlite"]
    return str(caught.value)


def _rewrite_part(path: Path, edit: Any) -> None:
    lines = gzip.decompress(path.read_bytes()).decode().splitlines()
    path.write_bytes(gzip.compress(("\n".join(edit(lines)) + "\n").encode()))


def test_merge_fails_when_a_listed_tar_has_no_part(corpus: Corpus, tmp_path: Path) -> None:
    (corpus.parts / (B + ".tsv.gz")).unlink()
    message = _expect_merge_failure(corpus, tmp_path, "no part file")
    assert B in message and A not in message


def test_merge_fails_on_a_part_that_is_not_listed(corpus: Corpus, tmp_path: Path) -> None:
    corpus.metadata = [m for m in corpus.metadata if m["archive_name"] != B]
    message = _expect_merge_failure(corpus, tmp_path, "not in the tar list")
    assert B in message


def test_merge_fails_on_count_mismatch_and_reports_the_difference(
    corpus: Corpus, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    corpus.metadata[0]["num_predicted_structures"] = 13
    message = _expect_merge_failure(corpus, tmp_path, "expected 13")
    assert f"{A}: 12 cif members, expected 13" in message
    assert f"{A}: 12 pdb members, expected 13" in message
    assert "difference -1" in message
    assert "FAILED" in capsys.readouterr().out


def test_merge_fails_when_one_format_is_short(corpus: Corpus, tmp_path: Path) -> None:
    def drop_one_pdb(lines: list[str]) -> list[str]:
        victim = next(i for i, line in enumerate(lines) if "\tpdb\t" in line)
        trailer = json.loads(lines[-1].split("\t", 1)[1])
        trailer["n_pdb"] -= 1
        return lines[:victim] + lines[victim + 1 : -1] + ["#end\t" + json.dumps(trailer)]

    _rewrite_part(corpus.parts / (A + ".tsv.gz"), drop_one_pdb)
    message = _expect_merge_failure(corpus, tmp_path, "11 pdb members, expected 12")
    assert "cif members" not in message


def test_merge_fails_when_a_swissprot_tar_holds_the_wrong_format(
    corpus: Corpus, tmp_path: Path
) -> None:
    def flip(lines: list[str]) -> list[str]:
        body = [line.replace("\tcif\t", "\tpdb\t") for line in lines[1:-1]]
        return [lines[0]] + body + ['#end\t{"n_cif": 0, "n_pdb": 9, "n_skipped": 0}']

    _rewrite_part(corpus.parts / (SC + ".tsv.gz"), flip)
    message = _expect_merge_failure(corpus, tmp_path, "0 cif members, expected 9")
    assert "9 pdb members, expected 0" in message


def test_merge_fails_on_a_part_with_zero_rows(corpus: Corpus, tmp_path: Path) -> None:
    empty = "UP000000003_333_EMPTY_v6.tar"
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w") as archive:
        info = tarfile.TarInfo("README.txt")
        info.size = 5
        archive.addfile(info, io.BytesIO(b"hello"))
    corpus.add(empty, Tar(buffer.getvalue(), {}, ["README.txt"]), n=0, type_="proteome")
    assert read_part(corpus.parts / (empty + ".tsv.gz"))[2] == ["README.txt"]
    message = _expect_merge_failure(corpus, tmp_path, "zero rows")
    assert empty in message


def test_merge_fails_on_a_part_that_was_cut_short(corpus: Corpus, tmp_path: Path) -> None:
    _rewrite_part(corpus.parts / (A + ".tsv.gz"), lambda lines: lines[:-1])
    _expect_merge_failure(corpus, tmp_path, "no #end line")


def test_merge_fails_on_a_truncated_gzip(corpus: Corpus, tmp_path: Path) -> None:
    path = corpus.parts / (A + ".tsv.gz")
    path.write_bytes(path.read_bytes()[:-20])
    _expect_merge_failure(corpus, tmp_path, "unreadable part file")


def test_merge_fails_on_duplicate_rows(corpus: Corpus, tmp_path: Path) -> None:
    def duplicate(lines: list[str]) -> list[str]:
        body = lines[1:-1]
        first_cif = next(i for i, line in enumerate(body) if "\tcif\t" in line)
        last_cif = max(i for i, line in enumerate(body) if "\tcif\t" in line)
        body[last_cif] = body[first_cif]
        return [lines[0]] + body + [lines[-1]]

    _rewrite_part(corpus.parts / (A + ".tsv.gz"), duplicate)
    _expect_merge_failure(corpus, tmp_path, "more than once")


def test_merge_fails_on_a_part_for_another_tar(corpus: Corpus, tmp_path: Path) -> None:
    (corpus.parts / (A + ".tsv.gz")).write_bytes((corpus.parts / (B + ".tsv.gz")).read_bytes())
    _expect_merge_failure(corpus, tmp_path, "header names")


def test_merge_reports_every_problem_at_once(corpus: Corpus, tmp_path: Path) -> None:
    (corpus.parts / (B + ".tsv.gz")).unlink()
    corpus.metadata[0]["num_predicted_structures"] = 99
    message = _expect_merge_failure(corpus, tmp_path, "3 problem")
    assert "no part file" in message and "expected 99" in message


def test_known_count_exception_is_explicit(
    corpus: Corpus, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    corpus.metadata[0]["num_predicted_structures"] = 13
    monkeypatch.setitem(_build.KNOWN_COUNT_EXCEPTIONS, A, {"cif": 12, "pdb": 12})
    _build.merge_parts(corpus.parts, tmp_path / "index-v6.sqlite", metadata=corpus.metadata)


# --------------------------------------------------------------------------
# verify_index, summary, rebuild check, command line


@pytest.fixture()
def index(corpus: Corpus, tmp_path: Path) -> Path:
    return _build.merge_parts(corpus.parts, tmp_path / "index-v6.sqlite", metadata=corpus.metadata)


def test_verify_index_reads_members_back(corpus: Corpus, index: Path) -> None:
    corpus.site.requests.clear()
    summary = _build.verify_index(index, 10, per_tar=2, seed=7)
    assert summary["passed"] == summary["checked"] >= 10
    assert summary["failed"] == []
    assert summary["by_format"]["cif"] > 0 and summary["by_format"]["pdb"] > 0
    ranged = [r for r in corpus.site.requests if r[2]]
    assert len(ranged) == summary["checked"]  # one range request per member
    assert {r[1] for r in ranged} == {A, B, SC, SP}  # per_tar reaches every tar
    assert _build.verify_index(index, 1000)["checked"] == 52  # capped at the index size


def test_verify_index_fails_on_wrong_offsets(index: Path) -> None:
    conn = sqlite3.connect(index)
    conn.execute("UPDATE members SET offset = offset + 1")
    conn.commit()
    conn.close()
    with pytest.raises(_build.VerifyError, match="spot checks failed"):
        _build.verify_index(index, 5, seed=1)


def test_verify_index_fails_when_rows_point_at_the_wrong_entry(index: Path) -> None:
    conn = sqlite3.connect(index)
    rows = conn.execute(
        "SELECT entry_id, offset, length FROM members WHERE tar_id = 1 AND fmt = 'cif' "
        "ORDER BY entry_id LIMIT 2"
    ).fetchall()
    conn.execute(
        "UPDATE members SET offset = ?, length = ? WHERE tar_id = 1 AND fmt = 'cif' "
        "AND entry_id = ?", (rows[1][1], rows[1][2], rows[0][0]),
    )
    conn.execute("DELETE FROM members WHERE NOT (tar_id = 1 AND fmt = 'cif' AND entry_id = ?)",
                 (rows[0][0],))
    conn.commit()
    conn.close()
    with pytest.raises(_build.VerifyError, match="first line is"):
        _build.verify_index(index, 5)


def test_verify_index_fails_when_the_tar_changed(corpus: Corpus, index: Path) -> None:
    corpus.site.etags[A] = '"different"'
    with pytest.raises(_build.VerifyError, match="ETag"):
        _build.verify_index(index, 0, per_tar=1)


def test_check_model_text() -> None:
    rng = random.Random(0)
    check = _build._check_model_text
    assert check("AF-P69905-F1", "cif", cif_text("AF-P69905-F1", rng)) is None
    assert check("AF-P69905-F1", "pdb", pdb_text("AF-P69905-F1", rng)) is None
    assert check("AF-P69905-F1", "cif", cif_text("AF-P69906-F1", rng)) is not None
    assert check("AF-P69905-F1", "cif", "data_AF-P69905-F1\n#\n") is not None
    assert check("AF-P69905-F1", "pdb", pdb_text("AF-P69906-F1", rng)) is not None
    assert check("AF-P6990-F1", "pdb", pdb_text("AF-P69905-F1", rng)) is not None
    assert check("AF-P69905-F1", "pdb", "HEADER\nTITLE (P69905)\nEND\n") is not None
    assert check("P69905", "cif", cif_text("AF-P69905-F1", rng)) is not None


def test_summary_and_rebuild_check(corpus: Corpus, index: Path, tmp_path: Path) -> None:
    import hashlib

    summary = _build.write_summary(index, tmp_path / "index-v6.json")
    assert json.loads((tmp_path / "index-v6.json").read_text()) == summary
    gz = index.with_name(index.name + ".gz")
    assert summary["index_gz_sha256"] == hashlib.sha256(gz.read_bytes()).hexdigest()
    assert summary["index_gz_bytes"] == gz.stat().st_size
    assert summary["index_bytes"] == index.stat().st_size
    assert summary["total_rows"] == 52
    assert summary["schema_version"] == "1" and summary["built_at"]
    assert summary["tars"][0] == {
        "name": A, "size_bytes": len(corpus.tars[A].data), "last_modified": LAST_MODIFIED,
        "etag": corpus.site.etags[A], "n_members": 24,
    }

    def reasons(published: Any) -> list[str]:
        return _build.needs_rebuild(
            published, metadata=corpus.metadata, base_url=corpus.site.base_url
        )

    corpus.site.requests.clear()
    assert reasons(summary) == []
    assert [r[0] for r in corpus.site.requests] == ["HEAD"] * 4  # HEAD only, no downloads
    assert reasons(None) == ["no published index-v6.json"]
    assert "schema_version" in reasons({**summary, "schema_version": "0"})[0]
    assert reasons({**summary, "tars": summary["tars"][1:]}) == [f"{A}: not in the published index"]
    corpus.metadata.pop()
    assert reasons(summary) == [f"{SP}: no longer listed by EBI"]
    corpus.site.etags[B] = '"new"'
    assert any(r.startswith(f"{B}: etag") for r in reasons(summary))


def test_command_line(
    site: Site, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    tar = make_tar(entry_ids("C", 4), extras=False)
    site.add(A, tar.data)
    listing = [{"archive_name": A, "num_predicted_structures": 4, "size_bytes": len(tar.data),
                "type": "proteome", "reference_proteome": "UP000000001"}]
    site.files["download_metadata.json"] = json.dumps(listing).encode()
    site.etags["download_metadata.json"] = '"m"'
    monkeypatch.setattr(
        _build, "METADATA_URL", site.base_url.rsplit("/", 1)[0] + "/download_metadata.json"
    )
    github_output = tmp_path / "github_output"
    monkeypatch.setenv("GITHUB_OUTPUT", str(github_output))
    parts, index, summary = tmp_path / "parts", tmp_path / "index-v6.sqlite", tmp_path / "s.json"

    def outputs() -> dict[str, str]:
        pairs = (line.split("=", 1) for line in github_output.read_text().splitlines())
        result = {key: value for key, value in pairs}
        github_output.unlink()
        return result

    assert _build.main(["plan", "--base-url", site.base_url]) == 0
    assert outputs() == {"rebuild": "true", "tars": json.dumps([A], separators=(",", ":"))}
    assert _build.main(["build-part", A, "-o", str(parts), "--base-url", site.base_url]) == 0
    assert _build.main(["merge", str(parts), "-o", str(index), "--only", A]) == 0
    assert _build.main(["verify", str(index), "-n", "3", "--per-tar", "1"]) == 0
    assert _build.main(["summary", str(index), "-o", str(summary)]) == 0
    capsys.readouterr()
    args = ["plan", "--base-url", site.base_url, "--published", str(summary)]
    assert _build.main(args) == 0
    assert outputs()["rebuild"] == "false"
    assert "up to date" in capsys.readouterr().out
    assert _build.main(args + ["--force"]) == 0
    assert outputs()["rebuild"] == "true"
    site.etags[A] = '"new"'
    assert _build.main(args) == 0
    assert outputs()["rebuild"] == "true"


# --------------------------------------------------------------------------
# Live tests against EBI (SCIGANTIC_ALPHAFOLD_LIVE=1)


def _live_head_names(tar: str, n_bytes: int = 200_000) -> list[str]:
    request = urllib.request.Request(
        f"{_build.V6_BASE_URL}/{tar}",
        headers={"User-Agent": _build.USER_AGENT, "Range": f"bytes=0-{n_bytes - 1}"},
    )
    with urllib.request.urlopen(request, timeout=60) as resp:
        assert resp.status == 206
        head = resp.read()
    names: list[str] = []
    try:
        with tarfile.open(fileobj=io.BytesIO(head), mode="r|") as archive:
            for member in archive:
                names.append(member.name)
    except tarfile.ReadError:
        pass  # the head ends in the middle of a member
    return names


@live
@pytest.mark.live
@pytest.mark.parametrize("fmt", ["cif", "pdb"])
def test_live_swissprot_member_naming(fmt: str) -> None:
    names = _live_head_names(f"swissprot_{fmt}_v6.tar")
    assert len(names) >= 3
    for name in names:
        match = _build._MEMBER_RE.fullmatch(name)
        assert match is not None, name
        assert match.group(2) == fmt


@live
@pytest.mark.live
def test_live_build_merge_verify_smallest_tar(tmp_path: Path) -> None:
    """Streams the smallest v6 tar (HELPY, 173 MB) from EBI."""
    tar = "UP000000429_85962_HELPY_v6.tar"
    entry = next(e for e in _build.fetch_download_metadata() if e["archive_name"] == tar)
    part = _build.build_part(tar, tmp_path / "parts")
    _, rows, skipped, trailer = read_part(part)
    n = entry["num_predicted_structures"]
    assert trailer["n_cif"] == trailer["n_pdb"] == n
    assert len(rows) == 2 * n and skipped == []
    index = _build.merge_parts(tmp_path / "parts", tmp_path / "index-v6.sqlite", metadata=[entry])
    assert _build.verify_index(index, 10)["passed"] == 10
