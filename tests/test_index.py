from __future__ import annotations

import dataclasses
import gzip
import hashlib
import http.client
import json
import os
import subprocess
import sys
import threading
import time
import urllib.request
from pathlib import Path

import pytest

import scigantic_alphafold as af
from scigantic_alphafold import _index

from conftest import PROTEOME_TAR, SWISSPROT_TAR, Afdb


# -- import hygiene ---------------------------------------------------------------------------


def test_import_touches_neither_network_nor_cache(tmp_path: Path) -> None:
    cache = tmp_path / "cache"
    code = (
        "import socket, sqlite3\n"
        "def boom(*a, **k): raise AssertionError('import opened a socket or a database')\n"
        "socket.socket.connect = boom\n"
        "sqlite3.connect = boom\n"
        "import scigantic_alphafold as af\n"
        "print(af.__version__)\n"
    )
    env = {**os.environ, "SCIGANTIC_ALPHAFOLD_CACHE": str(cache), "HOME": str(tmp_path), "USERPROFILE": str(tmp_path)}
    env.pop("SCIGANTIC_ALPHAFOLD_INDEX", None)
    done = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env)
    assert done.returncode == 0, done.stderr
    assert done.stdout.strip() == af.__version__
    assert not cache.exists()
    assert not (tmp_path / ".cache").exists()


def test_public_api_surface() -> None:
    contract = {
        "fetch", "fetch_bytes", "download", "fetch_many", "locate", "has", "entries", "tars",
        "index_info", "refresh_index", "metadata", "proteomes",
    }  # fmt: skip
    extras = {"__version__", "Member", "Tar"}
    errors = {name for name in af.__all__ if name.endswith("Error")}
    assert set(af.__all__) == contract | extras | errors
    assert {"NotInBulkSetError", "StaleIndexError", "AmbiguousEntryError", "FetchManyError"} <= errors
    for name in af.__all__:
        assert hasattr(af, name)
    assert issubclass(af.NotInBulkSetError, LookupError)
    assert issubclass(af.StaleIndexError, RuntimeError)
    assert issubclass(af.AmbiguousEntryError, LookupError)
    assert issubclass(af.InvalidIdError, ValueError)


def test_member_and_tar_are_frozen() -> None:
    member = af.Member("AF-P69905-F1", "cif", "t.tar", "http://x/t.tar", 512, 10)
    assert [f.name for f in dataclasses.fields(af.Member)] == ["entry_id", "fmt", "tar", "url", "offset", "length"]
    with pytest.raises(dataclasses.FrozenInstanceError):
        member.offset = 1  # type: ignore[misc]
    assert dataclasses.is_dataclass(af.Tar) and af.Tar.__dataclass_params__.frozen  # type: ignore[attr-defined]


# -- id forms ---------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "text",
    [
        "P69905",
        "p69905",
        "  P69905\n",
        "AF-P69905-F1",
        "af-p69905-f1",
        "\tAF-P69905-F1  ",
        "P69905-F1",
        "AF-P69905-F1-model_v6.cif.gz",
        "AF-P69905-F1-model_v6.cif",
        "AF-P69905-F1-model_v6",
    ],
)
def test_every_accepted_form_of_a_plain_accession(afdb: Afdb, text: str) -> None:
    member = _index.get_index().resolve(text, "cif")
    assert member.entry_id == "AF-P69905-F1"
    assert af.has(text)


@pytest.mark.parametrize("text", ["P12345-2", "p12345-2", " AF-P12345-2-F1 ", "af-p12345-2-f1", "P12345-2-F1"])
def test_isoform_accession_forms(afdb: Afdb, text: str) -> None:
    assert _index.get_index().resolve(text, "cif").entry_id == "AF-P12345-2-F1"


@pytest.mark.parametrize("text", ["", "   ", "P69905 F1", "AF-", "P69905/..", "AF-P69905-F", "P69905-F1-F2", "P6'9905"])
def test_malformed_ids_are_rejected(afdb: Afdb, text: str) -> None:
    with pytest.raises(af.InvalidIdError):
        af.locate(text)
    with pytest.raises(ValueError):
        af.has(text)


def test_non_string_id_is_rejected(afdb: Afdb) -> None:
    with pytest.raises(af.InvalidIdError):
        af.locate(69905)  # type: ignore[arg-type]


def test_bad_format_is_rejected(afdb: Afdb) -> None:
    with pytest.raises(ValueError, match="'cif' or 'pdb'"):
        af.locate("P69905", fmt="mmtf")


# -- which entry a bare accession means ---------------------------------------------------------


def test_bare_accession_means_the_canonical_entry_not_an_isoform(afdb: Afdb) -> None:
    assert _index.get_index().resolve("P12345", "cif").entry_id == "AF-P12345-F1"


def test_locate_returns_everything_for_an_accession(afdb: Afdb) -> None:
    assert [m.entry_id for m in af.locate("P12345")] == ["AF-P12345-F1", "AF-P12345-2-F1"]
    assert [m.entry_id for m in af.locate("P12345-2")] == ["AF-P12345-2-F1"]
    # Fragments in numeric order, not text order (F10 after F2).
    assert [m.entry_id for m in af.locate("Q00001")] == ["AF-Q00001-F1", "AF-Q00001-F2", "AF-Q00001-F10"]
    assert [(m.fmt, m.tar) for m in af.locate("P69905")] == [
        ("cif", PROTEOME_TAR),
        ("cif", SWISSPROT_TAR),
        ("pdb", PROTEOME_TAR),
    ]
    assert [m.fmt for m in af.locate("P69905", fmt="pdb")] == ["pdb"]


def test_locate_does_not_match_a_longer_accession(afdb: Afdb) -> None:
    # P1234 is a prefix of P12345 as text, but not the same accession.
    assert af.locate("P1234") == []
    assert af.locate("P12345-20") == []
    assert not af.has("P1234")


def test_fragmented_accession_is_ambiguous(afdb: Afdb) -> None:
    with pytest.raises(af.AmbiguousEntryError) as caught:
        _index.get_index().resolve("Q00001", "cif")
    assert caught.value.candidates == ("AF-Q00001-F1", "AF-Q00001-F2", "AF-Q00001-F10")
    message = str(caught.value)
    assert "3 fragments" in message
    for candidate in caught.value.candidates:
        assert candidate in message
    # Asking for a fragment by id is not ambiguous.
    assert _index.get_index().resolve("AF-Q00001-F2", "cif").entry_id == "AF-Q00001-F2"
    assert _index.get_index().resolve("q00001-f10", "cif").entry_id == "AF-Q00001-F10"


def test_isoforms_only_is_ambiguous_not_a_silent_pick(afdb: Afdb) -> None:
    with pytest.raises(af.AmbiguousEntryError) as caught:
        _index.get_index().resolve("A0A000ISO1", "cif")
    assert caught.value.candidates == ("AF-A0A000ISO1-2-F1", "AF-A0A000ISO1-3-F1")
    assert "only isoform entries" in str(caught.value)
    assert af.has("A0A000ISO1")
    assert _index.get_index().resolve("A0A000ISO1-3", "cif").entry_id == "AF-A0A000ISO1-3-F1"


def test_missing_fragment_of_a_known_accession(afdb: Afdb) -> None:
    with pytest.raises(af.NotInBulkSetError):
        _index.get_index().resolve("AF-P69905-F2", "cif")
    assert not af.has("AF-P69905-F2")


# -- the same entry in two tars -------------------------------------------------------------------


def test_same_entry_in_two_tars_picks_by_tar_name(afdb: Afdb) -> None:
    copies = af.locate("P69905", fmt="cif")
    assert [m.tar for m in copies] == [PROTEOME_TAR, SWISSPROT_TAR]
    chosen = _index.get_index().resolve("P69905", "cif")
    assert chosen == copies[0]
    assert chosen.tar == PROTEOME_TAR
    assert (chosen.offset, chosen.length) == afdb.where[PROTEOME_TAR][("AF-P69905-F1", "cif")]
    # The Swiss-Prot tar has the lower tar_id in this index, so the rule is not insertion order.


# -- not in the bulk set ----------------------------------------------------------------------------


def test_not_in_bulk_set_message(afdb: Afdb) -> None:
    with pytest.raises(af.NotInBulkSetError) as caught:
        _index.get_index().resolve("A0A999ZZZ9", "cif")
    message = str(caught.value)
    assert "A0A999ZZZ9" in message
    assert "Swiss-Prot plus 46 proteome" in message
    assert "per-file endpoint is not used" in message
    assert isinstance(caught.value, LookupError)
    assert af.locate("A0A999ZZZ9") == []
    assert not af.has("A0A999ZZZ9")


def test_wrong_format_message_names_the_formats_that_exist(afdb: Afdb) -> None:
    with pytest.raises(af.NotInBulkSetError, match="only as cif, not as pdb"):
        _index.get_index().resolve("P11111", "pdb")


# -- tars, entries, info ----------------------------------------------------------------------------


def test_tars(afdb: Afdb) -> None:
    listed = af.tars()
    assert [t.name for t in listed] == [PROTEOME_TAR, SWISSPROT_TAR]
    proteome = listed[0]
    assert proteome.url == afdb.tar_url(PROTEOME_TAR)
    assert proteome.size_bytes == len(afdb.tars[PROTEOME_TAR])
    assert (proteome.type, proteome.reference_proteome, proteome.tax_id) == ("proteome", "UP000000001", 1)
    assert proteome.n_members == 9
    assert listed[1].reference_proteome is None and listed[1].tax_id is None


def test_entries(afdb: Afdb, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(_index, "_PAGE", 2)  # force several pages
    everything = list(af.entries())
    assert everything == sorted(everything)
    assert len(everything) == len(set(everything)) == 9  # P69905 is in two tars, listed once
    assert list(af.entries(fmt="pdb")) == ["AF-P69905-F1"]
    assert list(af.entries(tar=SWISSPROT_TAR)) == ["AF-P11111-F1", "AF-P69905-F1"]
    assert list(af.entries(tar=af.tars()[1])) == ["AF-P11111-F1", "AF-P69905-F1"]
    assert len(list(af.entries(proteome="UP000000001"))) == 8
    assert list(af.entries(proteome="up000000001", fmt="pdb")) == ["AF-P69905-F1"]


def test_entries_rejects_unknown_filters_before_iterating(afdb: Afdb) -> None:
    with pytest.raises(af.NotInBulkSetError, match="UP999999999"):
        af.entries(proteome="UP999999999")
    with pytest.raises(af.NotInBulkSetError, match="nope.tar"):
        af.entries(tar="nope.tar")
    with pytest.raises(ValueError, match="not both"):
        af.entries(proteome="UP000000001", tar=PROTEOME_TAR)
    with pytest.raises(ValueError):
        af.entries(fmt="xyz")


def test_index_info(afdb: Afdb) -> None:
    info = af.index_info()
    assert info["path"] == str(afdb.index_path)
    assert info["schema_version"] == "1"
    assert info["afdb_version"] == "6"
    assert info["built_at"] == "2026-11-20T00:00:00Z"
    assert info["tars"] == 2
    assert info["members"] == 11
    assert info["size_bytes"] == afdb.index_path.stat().st_size


# -- where the index comes from ----------------------------------------------------------------------


def test_local_sqlite_is_opened_in_place_and_read_only(afdb: Afdb) -> None:
    assert af.index_info()["path"] == str(afdb.index_path)
    assert not afdb.cache_dir.exists()
    conn = _index.get_index()._open()
    with pytest.raises(Exception, match="readonly"):
        conn.execute("DELETE FROM members")
    assert af.refresh_index() == afdb.index_path


def test_local_gz_is_unpacked_into_the_cache(afdb: Afdb) -> None:
    afdb.use(afdb.root / "index-v6.sqlite.gz")
    assert af.has("P69905")
    path = Path(af.index_info()["path"])
    assert path.parent == afdb.cache_dir
    assert path.read_bytes() == afdb.index_path.read_bytes()
    assert [p.name for p in afdb.cache_dir.iterdir()] == [path.name]
    assert af.refresh_index() == path


def test_local_gz_detected_by_content_not_by_name(afdb: Afdb) -> None:
    odd = afdb.root / "snapshot.bin"
    odd.write_bytes(gzip.compress(afdb.index_path.read_bytes()))
    afdb.use(odd)
    assert af.has("P69905")


@pytest.mark.parametrize("name", ["index-v6.sqlite.gz", "index-v6.sqlite"])
def test_index_from_a_url(afdb: Afdb, name: str) -> None:
    afdb.use(f"{afdb.server.base}/index/{name}")
    assert afdb.server.count(f"/index/{name}") == 0
    assert af.has("P69905")
    path = Path(af.index_info()["path"])
    assert path.parent == afdb.cache_dir
    assert path.read_bytes() == afdb.index_path.read_bytes()
    assert [p.name for p in afdb.cache_dir.iterdir()] == [path.name]  # no temp files left
    request = afdb.server.requests[0]
    assert request["user_agent"].startswith(f"scigantic-alphafold/{af.__version__} (+https://github.com/")
    # Cached: a new process (simulated by a reset) does not download again.
    _index.reset()
    assert af.has("P69905")
    assert afdb.server.count(f"/index/{name}") == 1


def test_default_url_caches_under_a_fixed_name(afdb: Afdb, monkeypatch: pytest.MonkeyPatch) -> None:
    url = f"{afdb.server.base}/index/index-v6.sqlite.gz"
    monkeypatch.setattr(_index, "DEFAULT_INDEX_URL", url)
    afdb.use(None)
    assert af.index_info()["path"] == str(afdb.cache_dir / "index-v6.sqlite")
    assert af.index_info()["source"] == url


def test_refresh_downloads_again_and_replaces(afdb: Afdb) -> None:
    path = "/index/index-v6.sqlite.gz"
    afdb.use(afdb.server.base + path)
    assert not af.has("Z99999")
    newer = afdb.write_index(afdb.root / "newer.sqlite")
    import sqlite3

    conn = sqlite3.connect(newer)
    conn.execute("INSERT INTO members VALUES ('AF-Z99999-F1', 'cif', 1, 512, 10)")
    conn.commit()
    conn.close()
    afdb.server.files[path] = gzip.compress(newer.read_bytes())
    assert not af.has("Z99999")  # still the cached copy
    refreshed = af.refresh_index()
    assert refreshed.parent == afdb.cache_dir
    assert afdb.server.count(path) == 2
    assert af.has("Z99999")
    assert [p.name for p in afdb.cache_dir.iterdir()] == [refreshed.name]


def test_default_index_404_says_not_published(afdb: Afdb, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(_index, "DEFAULT_INDEX_URL", f"{afdb.server.base}/releases/index-v6.sqlite.gz")
    afdb.use(None)
    with pytest.raises(af.IndexUnavailableError) as caught:
        af.has("P69905")
    message = str(caught.value)
    assert "could not be downloaded" in message and "returned 404" in message
    assert "may not be published yet or may be mid-update" in message
    assert "has not been published yet" not in message  # a 404 does not prove that
    # The advice is a partial build that works with one tar: merge needs --only for that.
    tar = "UP000000429_85962_HELPY_v6.tar"
    assert f"scigantic-alphafold build-part {tar} -o parts" in message
    assert f"scigantic-alphafold merge parts -o index-v6.sqlite --only {tar}" in message
    assert "SCIGANTIC_ALPHAFOLD_INDEX=index-v6.sqlite" in message
    assert not any(afdb.cache_dir.iterdir())  # nothing half-written is left behind


def test_custom_url_404_and_500(afdb: Afdb) -> None:
    afdb.use(f"{afdb.server.base}/index/missing.sqlite.gz")
    with pytest.raises(af.IndexUnavailableError, match="HTTP 404"):
        af.has("P69905")
    afdb.server.fail_first["/index/index-v6.sqlite.gz"] = [500]
    afdb.use(f"{afdb.server.base}/index/index-v6.sqlite.gz")
    with pytest.raises(af.IndexUnavailableError, match="HTTP 500"):
        af.has("P69905")
    assert af.has("P69905")  # the next call tries again and succeeds


def test_missing_local_path(afdb: Afdb) -> None:
    afdb.use(afdb.root / "no-such-index.sqlite")
    with pytest.raises(af.IndexUnavailableError, match="not a file"):
        af.has("P69905")


def test_unknown_schema_version(afdb: Afdb) -> None:
    afdb.use(afdb.write_index(afdb.root / "future.sqlite", schema_version="2"))
    with pytest.raises(af.IndexSchemaError) as caught:
        af.has("P69905")
    message = str(caught.value)
    assert "schema_version 2" in message and "understands only 1" in message and "pip install -U" in message


def test_unknown_schema_version_from_a_url_is_not_cached(afdb: Afdb) -> None:
    future = afdb.write_index(afdb.root / "future.sqlite", schema_version="7")
    afdb.server.files["/index/future.sqlite.gz"] = gzip.compress(future.read_bytes())
    afdb.use(f"{afdb.server.base}/index/future.sqlite.gz")
    with pytest.raises(af.IndexSchemaError):
        af.has("P69905")
    assert not any(afdb.cache_dir.iterdir())


def test_a_file_that_is_not_an_index(afdb: Afdb) -> None:
    junk = afdb.root / "junk.sqlite"
    junk.write_bytes(b"<html>not an index</html>")
    afdb.use(junk)
    with pytest.raises(af.IndexUnavailableError, match="not a SQLite index"):
        af.has("P69905")
    afdb.server.files["/index/cut.sqlite.gz"] = afdb.server.files["/index/index-v6.sqlite.gz"][:200]
    afdb.use(f"{afdb.server.base}/index/cut.sqlite.gz")
    with pytest.raises(af.IndexUnavailableError, match="not a complete gzip"):
        af.has("P69905")
    assert not any(afdb.cache_dir.iterdir())


# -- concurrency --------------------------------------------------------------------------------------


def test_racing_downloads_leave_one_good_cache_file(afdb: Afdb) -> None:
    """Eight independent Index objects stand in for eight processes with a cold cache."""
    source = f"{afdb.server.base}/index/index-v6.sqlite.gz"
    afdb.server.delay = 0.05
    start = threading.Barrier(8)
    outcomes: list[object] = []

    def worker() -> None:
        index = _index.Index(source, afdb.cache_dir)
        start.wait()
        try:
            outcomes.append(index.resolve("P69905", "cif").entry_id)
        except Exception as exc:  # pragma: no cover - reported by the assert below
            outcomes.append(exc)
        finally:
            index.close()

    threads = [threading.Thread(target=worker) for _ in range(8)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(timeout=30)
    assert outcomes == ["AF-P69905-F1"] * 8
    files = list(afdb.cache_dir.iterdir())
    assert len(files) == 1
    assert files[0].read_bytes() == afdb.index_path.read_bytes()


def test_racing_processes_leave_one_good_cache_file(afdb: Afdb) -> None:
    source = f"{afdb.server.base}/index/index-v6.sqlite.gz"
    afdb.server.delay = 0.05
    code = "import scigantic_alphafold as af; print(af.locate('P69905', fmt='pdb')[0].entry_id)"
    env = {**os.environ, _index.ENV_INDEX: source, _index.ENV_CACHE: str(afdb.cache_dir)}
    running = [subprocess.Popen([sys.executable, "-c", code], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env) for _ in range(4)]
    for process in running:
        out, err = process.communicate(timeout=60)
        assert process.returncode == 0, err
        assert out.strip() == "AF-P69905-F1"
    files = list(afdb.cache_dir.iterdir())
    assert len(files) == 1
    assert files[0].read_bytes() == afdb.index_path.read_bytes()


@pytest.mark.skipif(sys.platform == "win32", reason="Windows cannot replace a file another connection has open")
def test_a_reader_survives_a_refresh_by_someone_else(afdb: Afdb) -> None:
    source = f"{afdb.server.base}/index/index-v6.sqlite.gz"
    reader = _index.Index(source, afdb.cache_dir)
    other = _index.Index(source, afdb.cache_dir)
    try:
        assert reader.has("P69905")
        other.refresh()
        assert reader.has("P69905")  # its open file was replaced, not truncated
        assert reader.resolve("P12345", "cif").entry_id == "AF-P12345-F1"
    finally:
        reader.close()
        other.close()


def test_one_index_shared_by_many_threads(afdb: Afdb) -> None:
    ids = ["P69905", "P12345", "P12345-2", "AF-Q00001-F2", "P11111"] * 40
    expect = ["AF-P69905-F1", "AF-P12345-F1", "AF-P12345-2-F1", "AF-Q00001-F2", "AF-P11111-F1"] * 40
    got: dict[int, str] = {}
    failures: list[Exception] = []

    def worker(position: int) -> None:
        try:
            got[position] = _index.get_index().resolve(ids[position], "cif").entry_id
        except Exception as exc:  # pragma: no cover
            failures.append(exc)

    threads = [threading.Thread(target=worker, args=(i,)) for i in range(len(ids))]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join(timeout=30)
    assert not failures
    assert [got[i] for i in range(len(ids))] == expect


# -- review fixes: paths ----------------------------------------------------------------------------------


def test_relative_index_path(afdb: Afdb, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(afdb.root)
    afdb.use("index-v6.sqlite")
    assert af.has("P69905")
    assert af.index_info()["path"] == str(afdb.index_path)
    assert af.fetch("P69905") == afdb.text("AF-P69905-F1")
    afdb.use(os.path.join(".", "index-v6.sqlite.gz"))
    assert af.has("P69905")


def test_relative_cache_dir(afdb: Afdb, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(afdb.root)
    monkeypatch.setenv(_index.ENV_CACHE, os.path.join("rel", "cache"))
    afdb.use(f"{afdb.server.base}/index/index-v6.sqlite.gz")
    assert af.has("P69905")
    path = Path(af.index_info()["path"])
    assert path.is_absolute() and path.parent == afdb.root / "rel" / "cache"
    assert af.refresh_index() == path


def test_tilde_is_expanded_in_both_paths(afdb: Afdb, monkeypatch: pytest.MonkeyPatch) -> None:
    home = afdb.root / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))
    monkeypatch.chdir(afdb.root)
    (home / "idx.sqlite.gz").write_bytes((afdb.root / "index-v6.sqlite.gz").read_bytes())
    monkeypatch.setenv(_index.ENV_CACHE, os.path.join("~", "af-cache"))
    afdb.use(os.path.join("~", "idx.sqlite.gz"))
    assert af.has("P69905")
    assert Path(af.index_info()["path"]).parent == home / "af-cache"
    assert not (afdb.root / "~").exists()  # no directory literally named "~"


@pytest.mark.parametrize("source", ["url", "gz"])
def test_unusable_cache_dir_is_a_library_error(afdb: Afdb, monkeypatch: pytest.MonkeyPatch, source: str) -> None:
    blocker = afdb.root / "a-file"
    blocker.write_text("not a directory")
    monkeypatch.setenv(_index.ENV_CACHE, str(blocker / "cache"))
    afdb.use(f"{afdb.server.base}/index/index-v6.sqlite.gz" if source == "url" else afdb.root / "index-v6.sqlite.gz")
    with pytest.raises(af.IndexUnavailableError, match="Cannot open the index"):
        af.has("P69905")
    with pytest.raises(af.AlphaFoldError):
        af.fetch("P69905")


@pytest.mark.parametrize("bad", ["bad\0path.sqlite", "", "   ", "."])
def test_odd_index_paths_never_escape_as_other_errors(afdb: Afdb, monkeypatch: pytest.MonkeyPatch, bad: str) -> None:
    index = _index.Index(bad or None, afdb.cache_dir)
    if not bad.strip():
        monkeypatch.setattr(_index, "DEFAULT_INDEX_URL", f"{afdb.server.base}/releases/index-v6.sqlite.gz")
        index = _index.Index(None, afdb.cache_dir)
    try:
        with pytest.raises(af.AlphaFoldError):
            index.has("P69905")
    finally:
        index.close()


# -- review fixes: refresh on a platform that cannot replace an open file ------------------------------------


def test_refresh_closes_its_own_connection_before_replacing(afdb: Afdb, monkeypatch: pytest.MonkeyPatch) -> None:
    """Simulated Windows: os.replace fails while this process has the cached index open."""
    afdb.use(f"{afdb.server.base}/index/index-v6.sqlite.gz")
    assert af.has("P69905")
    index = _index.get_index()
    assert index._conn is not None
    real_replace = os.replace
    held_open: list[bool] = []

    def windows_replace(src: str, dst: str) -> None:
        held_open.append(index._conn is not None)
        if index._conn is not None:
            raise PermissionError(13, "The process cannot access the file because it is being used by another process")
        real_replace(src, dst)

    monkeypatch.setattr(os, "replace", windows_replace)
    path = af.refresh_index()
    assert held_open == [False]
    assert af.has("P69905") and Path(af.index_info()["path"]) == path


def test_refresh_blocked_by_another_process_keeps_the_old_index(afdb: Afdb, monkeypatch: pytest.MonkeyPatch) -> None:
    afdb.use(f"{afdb.server.base}/index/index-v6.sqlite.gz")
    assert af.has("P69905")
    attempts: list[str] = []

    def always_locked(src: str, dst: str) -> None:
        attempts.append(dst)
        raise PermissionError(13, "being used by another process")

    with monkeypatch.context() as patched:
        patched.setattr(os, "replace", always_locked)
        with pytest.raises(af.IndexUnavailableError, match="Another process probably has it open"):
            af.refresh_index()
    assert len(attempts) == _index._REPLACE_ATTEMPTS
    assert af.has("P69905")  # the old file is opened again
    assert [p.name for p in afdb.cache_dir.iterdir()] == ["index-" + hashlib.sha256(f"{afdb.server.base}/index/index-v6.sqlite.gz".encode()).hexdigest()[:16] + ".sqlite"]


def test_first_open_uses_the_copy_another_process_holds(afdb: Afdb, monkeypatch: pytest.MonkeyPatch) -> None:
    """Two cold starts race on Windows: the loser cannot replace the winner's open file, and need not."""
    source = f"{afdb.server.base}/index/index-v6.sqlite.gz"
    winner = _index.Index(source, afdb.cache_dir)
    loser = _index.Index(source, afdb.cache_dir)
    try:
        assert winner.has("P69905")  # the winner finished and holds the file open
        real_is_file = Path.is_file
        looked: list[Path] = []

        def cold(self: Path) -> bool:
            # The loser looked before the winner had finished, so it downloads too.
            if self.parent == afdb.cache_dir and self.suffix == ".sqlite" and not looked:
                looked.append(self)
                return False
            return real_is_file(self)

        replaces: list[str] = []

        def locked(src: str, dst: str) -> None:
            replaces.append(dst)
            raise PermissionError(13, "being used by another process")

        with monkeypatch.context() as patched:
            patched.setattr(Path, "is_file", cold)
            patched.setattr(os, "replace", locked)
            assert loser.resolve("P69905", "cif").entry_id == "AF-P69905-F1"
        assert len(looked) == 1 and len(replaces) == 1
        assert len(list(afdb.cache_dir.iterdir())) == 1  # the loser's temp files are gone
    finally:
        winner.close()
        loser.close()


# -- review fixes: file-name ids -----------------------------------------------------------------------------


def test_file_name_id_decides_the_format_when_fmt_is_left_out(afdb: Afdb) -> None:
    index = _index.get_index()
    assert index.resolve("AF-P69905-F1-model_v6.pdb.gz").fmt == "pdb"
    assert index.resolve("af-p69905-f1-model_v6.pdb").fmt == "pdb"
    assert index.resolve("AF-P69905-F1-model_v6.cif.gz").fmt == "cif"
    assert index.resolve("AF-P69905-F1-model_v6").fmt == "cif"
    assert index.resolve("P69905").fmt == "cif"
    assert index.resolve("AF-P69905-F1-model_v6.pdb.gz", "pdb").fmt == "pdb"
    assert [m.fmt for m in af.locate("AF-P69905-F1-model_v6.pdb.gz")] == ["pdb"]
    assert af.has("AF-P69905-F1-model_v6.pdb.gz")
    assert not af.has("AF-P11111-F1-model_v6.pdb.gz")  # stored as cif only


def test_file_name_id_conflicting_with_fmt_is_rejected(afdb: Afdb) -> None:
    index = _index.get_index()
    with pytest.raises(af.InvalidIdError, match="is a pdb file name but fmt='cif'"):
        index.resolve("AF-P69905-F1-model_v6.pdb.gz", "cif")
    with pytest.raises(af.InvalidIdError, match="is a cif file name but fmt='pdb'"):
        af.locate("AF-P69905-F1-model_v6.cif", fmt="pdb")


@pytest.mark.parametrize("text", ["AF-P69905-F1-model_v4.cif", "AF-P69905-F1-model_v4", "P69905-F1-model_v7.pdb.gz"])
def test_file_name_id_of_another_version_is_rejected(afdb: Afdb, text: str) -> None:
    with pytest.raises(af.InvalidIdError, match="the index describes v6"):
        _index.get_index().resolve(text)
    with pytest.raises(af.InvalidIdError):
        af.locate(text)


def test_gz_without_a_format_is_not_an_id(afdb: Afdb) -> None:
    with pytest.raises(af.InvalidIdError):
        af.locate("AF-P69905-F1-model_v6.gz")


# -- review fixes: the index is checked, not trusted ------------------------------------------------------------


@pytest.mark.parametrize("value", ["6/../../../escaped", "v6", "", "6 "])
def test_afdb_version_must_be_digits(afdb: Afdb, value: str) -> None:
    afdb.use(afdb.edit_index("badversion.sqlite", f"UPDATE meta SET value = '{value}' WHERE key = 'afdb_version'"))
    with pytest.raises(af.IndexSchemaError) as caught:
        af.has("P69905")
    assert repr(value) in str(caught.value) and "afdb_version" in str(caught.value)


def test_missing_afdb_version(afdb: Afdb) -> None:
    afdb.use(afdb.edit_index("noversion.sqlite", "DELETE FROM meta WHERE key = 'afdb_version'"))
    with pytest.raises(af.IndexSchemaError, match="afdb_version None"):
        af.has("P69905")


def test_tar_url_outside_base_url_is_never_requested(afdb: Afdb) -> None:
    hostile = f"{afdb.server.base}/internal/admin"
    afdb.server.files["/internal/admin"] = afdb.tars[PROTEOME_TAR]
    afdb.use(afdb.edit_index("hostile.sqlite", f"UPDATE tars SET url = '{hostile}'"))
    with pytest.raises(af.IndexSchemaError) as caught:
        af.fetch("P69905")
    assert hostile in str(caught.value) and "outside its own base_url" in str(caught.value)
    assert afdb.server.requests == []


def test_local_server_is_refused_without_the_opt_out(afdb: Afdb, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(_index.ENV_ALLOW_ANY_HOST)
    with pytest.raises(af.IndexSchemaError) as caught:
        af.fetch("P69905")
    message = str(caught.value)
    assert afdb.server.base in message and "ftp.ebi.ac.uk" in message
    assert "SCIGANTIC_ALPHAFOLD_ALLOW_ANY_HOST=1" in message
    assert afdb.server.requests == []
    # Looking things up needs no host at all.
    assert af.has("P69905") and af.locate("P69905")
    monkeypatch.setenv(_index.ENV_ALLOW_ANY_HOST, "0")
    with pytest.raises(af.IndexSchemaError):
        af.fetch("P69905")


@pytest.mark.parametrize(
    ("base", "url", "expect"),
    [
        ("https://ftp.ebi.ac.uk/pub/databases/alphafold/v6", "https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/x.tar", None),
        ("https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/", "https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/x.tar", None),
        ("http://ftp.ebi.ac.uk/pub/v6", "http://ftp.ebi.ac.uk/pub/v6/x.tar", "not an https URL"),
        ("https://ftp.ebi.ac.uk.evil.example/v6", "https://ftp.ebi.ac.uk.evil.example/v6/x.tar", "not an https URL on ftp.ebi.ac.uk"),
        ("https://ftp.ebi.ac.uk@evil.example/v6", "https://ftp.ebi.ac.uk@evil.example/v6/x.tar", "not an https URL on ftp.ebi.ac.uk"),
        ("https://ftp.ebi.ac.uk/pub/v6", "https://ftp.ebi.ac.uk/pub/v6-other/x.tar", "outside its own base_url"),
        ("https://ftp.ebi.ac.uk/pub/v6", "https://evil.example/x.tar", "outside its own base_url"),
        ("", "https://ftp.ebi.ac.uk/pub/v6/x.tar", "not an https URL"),
    ],
)
def test_tar_url_policy_without_the_opt_out(afdb: Afdb, monkeypatch: pytest.MonkeyPatch, base: str, url: str, expect: str | None) -> None:
    monkeypatch.delenv(_index.ENV_ALLOW_ANY_HOST)
    afdb.use(afdb.edit_index("policy.sqlite", f"UPDATE meta SET value = '{base}' WHERE key = 'base_url'"))
    index = _index.get_index()
    if expect is None:
        index.check_tar_url(url)
    else:
        with pytest.raises(af.IndexSchemaError, match=expect):
            index.check_tar_url(url)


def test_opt_out_never_allows_other_schemes(afdb: Afdb) -> None:
    afdb.use(afdb.edit_index("file.sqlite", "UPDATE tars SET url = 'file:///etc/passwd'", "UPDATE meta SET value = 'file:///etc' WHERE key = 'base_url'"))
    with pytest.raises(af.IndexSchemaError, match="not http"):
        af.fetch("P69905")


def test_malformed_rows_are_refused(afdb: Afdb, tmp_path: Path) -> None:
    offset, length = afdb.where[PROTEOME_TAR][("AF-P69905-F1", "cif")]
    afdb.use(afdb.edit_index("rows.sqlite", f"INSERT INTO members VALUES ('AF-Z00001-../../../escaped', 'cif', 2, {offset}, {length})"))
    with pytest.raises(af.IndexSchemaError, match="malformed entry id"):
        af.download("Z00001", tmp_path / "out")
    assert not (tmp_path / "escaped").exists() and afdb.server.requests == []
    size = len(afdb.tars[PROTEOME_TAR])
    for name, statement in (
        ("past.sqlite", f"UPDATE members SET length = {size} WHERE entry_id = 'AF-P69905-F1'"),
        ("negative.sqlite", "UPDATE members SET offset = -5 WHERE entry_id = 'AF-P69905-F1'"),
        ("empty.sqlite", "UPDATE members SET length = 0 WHERE entry_id = 'AF-P69905-F1'"),
    ):
        afdb.use(afdb.edit_index(name, statement))
        with pytest.raises(af.IndexSchemaError, match="not inside the tar"):
            af.fetch("P69905")
    assert afdb.server.requests == []


# -- review fixes: the published checksum ------------------------------------------------------------------------


def use_default_url(afdb: Afdb, monkeypatch: pytest.MonkeyPatch, summary: object = None) -> str:
    monkeypatch.setattr(_index, "DEFAULT_INDEX_URL", f"{afdb.server.base}/index/index-v6.sqlite.gz")
    if summary is not None:
        body = summary if isinstance(summary, bytes) else json.dumps(summary).encode()
        afdb.server.files["/index/index-v6.json"] = body
    afdb.use(None)
    return hashlib.sha256(afdb.server.files["/index/index-v6.sqlite.gz"]).hexdigest()


def test_checksum_match(afdb: Afdb, monkeypatch: pytest.MonkeyPatch) -> None:
    digest = hashlib.sha256(afdb.server.files["/index/index-v6.sqlite.gz"]).hexdigest()
    use_default_url(afdb, monkeypatch, {"index_gz_sha256": digest.upper(), "schema_version": "1"})
    assert af.has("P69905")
    assert [r["path"] for r in afdb.server.requests] == ["/index/index-v6.sqlite.gz", "/index/index-v6.json"]
    assert afdb.server.requests[1]["user_agent"].startswith("scigantic-alphafold/")


def test_checksum_mismatch_fails_closed(afdb: Afdb, monkeypatch: pytest.MonkeyPatch) -> None:
    actual = use_default_url(afdb, monkeypatch, {"index_gz_sha256": "0" * 64})
    with pytest.raises(af.IndexUnavailableError) as caught:
        af.has("P69905")
    message = str(caught.value)
    assert "does not match the checksum" in message and "0" * 64 in message and actual in message
    assert "mid-update" in message
    assert not any(afdb.cache_dir.iterdir())  # nothing was installed


def test_checksum_mismatch_on_refresh_keeps_the_working_index(afdb: Afdb, monkeypatch: pytest.MonkeyPatch) -> None:
    use_default_url(afdb, monkeypatch)
    assert af.has("P69905")
    afdb.server.files["/index/index-v6.json"] = json.dumps({"index_gz_sha256": "f" * 64}).encode()
    with pytest.raises(af.IndexUnavailableError, match="does not match the checksum"):
        af.refresh_index()
    assert af.has("P69905")
    assert [p.name for p in afdb.cache_dir.iterdir()] == ["index-v6.sqlite"]


@pytest.mark.parametrize("summary", [None, {"schema_version": "1"}, {"index_gz_sha256": None}])
def test_no_published_checksum_proceeds(afdb: Afdb, monkeypatch: pytest.MonkeyPatch, summary: object) -> None:
    use_default_url(afdb, monkeypatch, summary)
    assert af.has("P69905")
    assert afdb.server.count("/index/index-v6.json") == 1


@pytest.mark.parametrize("summary", [b"<html>oops</html>", {"index_gz_sha256": "abc"}, {"index_gz_sha256": 5}])
def test_unusable_summary_fails_closed(afdb: Afdb, monkeypatch: pytest.MonkeyPatch, summary: object) -> None:
    use_default_url(afdb, monkeypatch, summary)
    with pytest.raises(af.IndexUnavailableError, match="cannot be verified"):
        af.has("P69905")
    assert not any(afdb.cache_dir.iterdir())


def test_summary_server_error_fails_closed(afdb: Afdb, monkeypatch: pytest.MonkeyPatch) -> None:
    use_default_url(afdb, monkeypatch, {"index_gz_sha256": "0" * 64})
    afdb.server.fail_first["/index/index-v6.json"] = [503]
    with pytest.raises(af.IndexUnavailableError, match="HTTP 503"):
        af.has("P69905")


def test_other_sources_are_not_checked_against_the_summary(afdb: Afdb) -> None:
    afdb.server.files["/index/index-v6.json"] = json.dumps({"index_gz_sha256": "0" * 64}).encode()
    afdb.use(f"{afdb.server.base}/index/index-v6.sqlite.gz")
    assert af.has("P69905")
    assert afdb.server.count("/index/index-v6.json") == 0


def test_real_summary_key_matches_the_builder() -> None:
    """The key this reader looks for is the one the builder writes."""
    import inspect

    from scigantic_alphafold import _build

    assert f'"{_index.SUMMARY_SHA256_KEY}"' in inspect.getsource(_build.write_summary)


# -- review fixes: smaller items ---------------------------------------------------------------------------------


def test_a_broken_chunked_download_is_a_library_error(afdb: Afdb, monkeypatch: pytest.MonkeyPatch) -> None:
    afdb.use(f"{afdb.server.base}/index/index-v6.sqlite.gz")

    def broken(request: urllib.request.Request, timeout: float = 0.0) -> object:
        raise http.client.IncompleteRead(b"partial", 100)

    with monkeypatch.context() as patched:
        patched.setattr(urllib.request, "urlopen", broken)
        with pytest.raises(af.IndexUnavailableError, match="IncompleteRead"):
            af.has("P69905")
    assert not any(afdb.cache_dir.iterdir())


def test_stale_part_files_are_swept_after_a_day(afdb: Afdb) -> None:
    afdb.cache_dir.mkdir()
    old = afdb.cache_dir / ".index-killed.part"
    fresh = afdb.cache_dir / ".index-inflight.part"
    unrelated = afdb.cache_dir / "notes.part"
    for path in (old, fresh, unrelated):
        path.write_bytes(b"x")
    two_days_ago = time.time() - 2 * 24 * 3600
    os.utime(old, (two_days_ago, two_days_ago))
    os.utime(unrelated, (two_days_ago, two_days_ago))
    afdb.use(f"{afdb.server.base}/index/index-v6.sqlite.gz")
    assert af.has("P69905")
    names = {p.name for p in afdb.cache_dir.iterdir()}
    assert ".index-killed.part" not in names
    assert {".index-inflight.part", "notes.part"} <= names


def test_not_in_bulk_set_suggests_a_refresh(afdb: Afdb) -> None:
    with pytest.raises(af.NotInBulkSetError, match=r"refresh_index\(\)"):
        _index.get_index().resolve("A0A999ZZZ9", "cif")


def test_check_index_waits_out_a_rename_in_progress(afdb: Afdb, monkeypatch: pytest.MonkeyPatch) -> None:
    """Windows refuses to open a file for a moment while another process renames over it."""
    calls = {"n": 0}

    def flaky_open(path: object, mode: str = "r") -> object:
        calls["n"] += 1
        if calls["n"] <= 2:
            raise PermissionError(13, "Permission denied")
        return open(path, mode)  # type: ignore[call-overload]

    monkeypatch.setattr(_index, "open", flaky_open, raising=False)
    monkeypatch.setattr(_index, "_sleep", lambda seconds: None)
    _index._check_index(afdb.index_path)
    assert calls["n"] == 3


def test_check_index_gives_up_on_a_file_that_stays_locked(afdb: Afdb, monkeypatch: pytest.MonkeyPatch) -> None:
    def locked_open(path: object, mode: str = "r") -> object:
        raise PermissionError(13, "Permission denied")

    monkeypatch.setattr(_index, "open", locked_open, raising=False)
    monkeypatch.setattr(_index, "_sleep", lambda seconds: None)
    with pytest.raises(af.IndexUnavailableError, match="Cannot read the index"):
        _index._check_index(afdb.index_path)
