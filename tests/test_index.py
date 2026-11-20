from __future__ import annotations

import dataclasses
import gzip
import subprocess
import sys
import threading
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
    env = {"PATH": "/usr/bin:/bin", "SCIGANTIC_ALPHAFOLD_CACHE": str(cache), "HOME": str(tmp_path)}
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
    assert "has not been published yet" in message
    assert "SCIGANTIC_ALPHAFOLD_INDEX=/path/to/index-v6.sqlite" in message
    assert "build-part" in message and "merge" in message
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
    env = {"PATH": "/usr/bin:/bin", _index.ENV_INDEX: source, _index.ENV_CACHE: str(afdb.cache_dir)}
    running = [subprocess.Popen([sys.executable, "-c", code], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env) for _ in range(4)]
    for process in running:
        out, err = process.communicate(timeout=60)
        assert process.returncode == 0, err
        assert out.strip() == "AF-P69905-F1"
    files = list(afdb.cache_dir.iterdir())
    assert len(files) == 1
    assert files[0].read_bytes() == afdb.index_path.read_bytes()


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
