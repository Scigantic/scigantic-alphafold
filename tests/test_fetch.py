from __future__ import annotations

import gzip
import io
import os
import sqlite3
import sys
import tarfile
import time
import urllib.request
from pathlib import Path

import pytest

import scigantic_alphafold as af
from scigantic_alphafold import _fetch, _index
from scigantic_alphafold._version import USER_AGENT

from conftest import PROTEOME_TAR, SWISSPROT_TAR, Afdb


def tar_requests(afdb: Afdb, name: str = PROTEOME_TAR) -> list[dict[str, object]]:
    return [r for r in afdb.server.requests if r["path"] == afdb.tar_path(name)]


# -- the happy path -------------------------------------------------------------------------------


def test_fetch_is_one_range_request_with_the_honest_user_agent(afdb: Afdb) -> None:
    text = af.fetch("P69905")
    assert text == afdb.text("AF-P69905-F1")
    assert text.startswith("data_AF-P69905-F1\n")
    offset, length = afdb.where[PROTEOME_TAR][("AF-P69905-F1", "cif")]
    assert afdb.server.requests == [
        {"path": afdb.tar_path(PROTEOME_TAR), "range": f"bytes={offset}-{offset + length - 1}", "user_agent": USER_AGENT}
    ]
    assert USER_AGENT == f"scigantic-alphafold/{af.__version__} (+https://github.com/Scigantic/scigantic-alphafold)"


def test_fetch_pdb(afdb: Afdb) -> None:
    text = af.fetch("p69905", fmt="pdb")
    assert text == afdb.text("AF-P69905-F1", "pdb")
    assert "\nATOM " in text


@pytest.mark.parametrize("text", ["P69905", " af-p69905-f1 ", "P69905-F1", "AF-P69905-F1-model_v6.cif.gz"])
def test_fetch_accepts_every_id_form(afdb: Afdb, text: str) -> None:
    assert af.fetch(text) == afdb.text("AF-P69905-F1")


def test_fetch_bytes_is_the_stored_gzip(afdb: Afdb) -> None:
    raw = af.fetch_bytes("P69905")
    offset, length = afdb.where[PROTEOME_TAR][("AF-P69905-F1", "cif")]
    assert raw == afdb.tars[PROTEOME_TAR][offset : offset + length]
    assert raw[:2] == b"\x1f\x8b"
    assert gzip.decompress(raw).decode() == afdb.text("AF-P69905-F1")
    # The same bytes tarfile itself would extract.
    with tarfile.open(fileobj=io.BytesIO(afdb.tars[PROTEOME_TAR])) as archive:
        extracted = archive.extractfile("AF-P69905-F1-model_v6.cif.gz")
        assert extracted is not None and extracted.read() == raw


def test_fetch_resolution_rules(afdb: Afdb) -> None:
    assert af.fetch("P12345") == afdb.text("AF-P12345-F1")
    assert af.fetch("P12345-2") == afdb.text("AF-P12345-2-F1")
    assert af.fetch("AF-Q00001-F2") == afdb.text("AF-Q00001-F2")
    with pytest.raises(af.AmbiguousEntryError, match="AF-Q00001-F1, AF-Q00001-F2, AF-Q00001-F10"):
        af.fetch("Q00001")
    with pytest.raises(af.AmbiguousEntryError, match="AF-A0A000ISO1-2-F1, AF-A0A000ISO1-3-F1"):
        af.fetch("A0A000ISO1")
    assert len(tar_requests(afdb)) == 3  # the ambiguous lookups never reached the network


def test_entry_in_two_tars_is_read_from_the_proteome_tar(afdb: Afdb) -> None:
    assert af.fetch("P69905") == afdb.text("AF-P69905-F1", tar=PROTEOME_TAR)
    assert len(tar_requests(afdb, SWISSPROT_TAR)) == 0
    # An entry that only the Swiss-Prot tar has is read from there.
    assert af.fetch("P11111") == afdb.text("AF-P11111-F1", tar=SWISSPROT_TAR)
    assert len(tar_requests(afdb, SWISSPROT_TAR)) == 1


def test_not_in_bulk_set_makes_no_request(afdb: Afdb) -> None:
    with pytest.raises(af.NotInBulkSetError, match="not in any AlphaFold DB v6 bulk tar") as caught:
        af.fetch("A0A999ZZZ9")
    assert "Swiss-Prot plus 46 proteome" in str(caught.value)
    assert "per-file endpoint is not used" in str(caught.value)
    with pytest.raises(af.NotInBulkSetError, match="only as cif, not as pdb"):
        af.fetch("P11111", fmt="pdb")
    assert afdb.server.requests == []


# -- stale index ------------------------------------------------------------------------------------


def test_changed_etag_is_stale(afdb: Afdb) -> None:
    afdb.server.etags[afdb.tar_path(PROTEOME_TAR)] = '"deadbeef-new"'
    with pytest.raises(af.StaleIndexError) as caught:
        af.fetch("P69905")
    message = str(caught.value)
    assert PROTEOME_TAR in message and '"deadbeef-new"' in message and "refresh_index()" in message
    assert isinstance(caught.value, RuntimeError)
    assert len(tar_requests(afdb)) == 1  # stale is not retried
    # The other tar is untouched and still readable.
    assert af.fetch("P11111") == afdb.text("AF-P11111-F1", tar=SWISSPROT_TAR)


def test_weak_and_strong_forms_of_the_same_etag_match(afdb: Afdb) -> None:
    path = afdb.tar_path(PROTEOME_TAR)
    afdb.server.etags[path] = f'W/"{len(afdb.tars[PROTEOME_TAR]):x}-test"'
    assert af.fetch("P69905") == afdb.text("AF-P69905-F1")


def test_changed_last_modified_is_stale(afdb: Afdb) -> None:
    afdb.server.last_modified[afdb.tar_path(PROTEOME_TAR)] = "Mon, 01 Jan 2029 00:00:00 GMT"
    with pytest.raises(af.StaleIndexError, match="Last-Modified is now Mon, 01 Jan 2029"):
        af.fetch("P69905")


def test_changed_size_is_stale(afdb: Afdb) -> None:
    path = afdb.tar_path(PROTEOME_TAR)
    grown = afdb.tars[PROTEOME_TAR] + b"\0" * 1024
    afdb.server.files[path] = grown
    afdb.server.etags[path] = f'"{len(afdb.tars[PROTEOME_TAR]):x}-test"'  # validators unchanged, size changed
    with pytest.raises(af.StaleIndexError, match="size is now"):
        af.fetch("P69905")


def test_corrupted_bytes_are_stale(afdb: Afdb) -> None:
    path = afdb.tar_path(PROTEOME_TAR)
    offset, length = afdb.where[PROTEOME_TAR][("AF-P69905-F1", "cif")]
    data = bytearray(afdb.tars[PROTEOME_TAR])
    data[offset + length // 2] ^= 0xFF
    afdb.server.files[path] = bytes(data)  # same size, same validators, different bytes
    with pytest.raises(af.StaleIndexError) as caught:
        af.fetch("P69905")
    assert "AF-P69905-F1" in str(caught.value) and "refresh_index()" in str(caught.value)
    with pytest.raises(af.StaleIndexError):
        af.fetch_bytes("P69905")
    assert len(tar_requests(afdb)) == 2


def test_shifted_offsets_are_stale(afdb: Afdb) -> None:
    """An index whose offsets are off (a tar repacked in place) lands on tar headers, not gzip."""
    shifted = afdb.root / "shifted.sqlite"
    afdb.write_index(shifted)
    conn = sqlite3.connect(shifted)
    conn.execute("UPDATE members SET offset = offset - 512")
    conn.commit()
    conn.close()
    afdb.use(shifted)
    with pytest.raises(af.StaleIndexError, match="not the gzip file of AF-P69905-F1"):
        af.fetch("P69905")


def test_range_past_the_end_is_stale(afdb: Afdb) -> None:
    afdb.server.files[afdb.tar_path(PROTEOME_TAR)] = afdb.tars[PROTEOME_TAR][:100]
    with pytest.raises(af.StaleIndexError, match="HTTP 416"):
        af.fetch("P69905")
    assert len(tar_requests(afdb)) == 1


def test_index_without_validators_still_fetches(afdb: Afdb) -> None:
    bare = afdb.root / "bare.sqlite"
    afdb.write_index(bare)
    conn = sqlite3.connect(bare)
    conn.execute("UPDATE tars SET etag = NULL, last_modified = NULL")
    conn.commit()
    conn.close()
    afdb.use(bare)
    afdb.server.etags[afdb.tar_path(PROTEOME_TAR)] = '"anything"'
    assert af.fetch("P69905") == afdb.text("AF-P69905-F1")


# -- a server that ignores Range ----------------------------------------------------------------------


def test_a_200_with_the_whole_tar_is_aborted_without_reading_it(afdb: Afdb) -> None:
    size = 2 * 1024**3
    afdb.server.huge[afdb.tar_path(PROTEOME_TAR)] = size
    started = time.monotonic()
    with pytest.raises(af.FetchError) as caught:
        af.fetch("P69905")
    elapsed = time.monotonic() - started
    assert caught.value.status == 200
    assert "HTTP 200 instead of 206" in str(caught.value)
    assert "closed without reading" in str(caught.value)
    assert not isinstance(caught.value, af.StaleIndexError)
    assert afdb.server.huge_done.wait(timeout=10)
    sent = afdb.server.huge_sent[0]
    # Only what fit in the socket buffers before the close was ever sent.
    assert sent < 32 * 1024**2, f"the server sent {sent} bytes"
    assert elapsed < 5
    assert len(tar_requests(afdb)) == 1  # and it is not retried


def test_a_small_200_is_rejected_too(afdb: Afdb) -> None:
    afdb.server.ignore_range.add(afdb.tar_path(PROTEOME_TAR))
    with pytest.raises(af.FetchError, match="HTTP 200 instead of 206"):
        af.fetch_bytes("P69905")
    assert len(tar_requests(afdb)) == 1


# -- retries ------------------------------------------------------------------------------------------


@pytest.mark.parametrize("failures", [[500], [503, 502], [500, 503, 504]])
def test_5xx_is_retried_with_backoff_until_it_succeeds(afdb: Afdb, failures: list[int]) -> None:
    afdb.server.fail_first[afdb.tar_path(PROTEOME_TAR)] = list(failures)
    assert af.fetch("P69905") == afdb.text("AF-P69905-F1")
    assert len(tar_requests(afdb)) == len(failures) + 1
    assert len(afdb.sleeps) == len(failures)
    assert afdb.sleeps == sorted(afdb.sleeps) and afdb.sleeps[0] >= 0.5  # grows


def test_5xx_gives_up_after_four_attempts(afdb: Afdb) -> None:
    afdb.server.fail_first[afdb.tar_path(PROTEOME_TAR)] = [503] * 10
    with pytest.raises(af.FetchError, match="after 4 attempts: HTTP 503"):
        af.fetch("P69905")
    assert len(tar_requests(afdb)) == 4
    assert len(afdb.sleeps) == 3


@pytest.mark.parametrize("status", [400, 403, 404, 429])
def test_4xx_is_never_retried(afdb: Afdb, status: int) -> None:
    afdb.server.fail_first[afdb.tar_path(PROTEOME_TAR)] = [status]
    with pytest.raises(af.FetchError, match=f"HTTP {status}") as caught:
        af.fetch("P69905")
    assert caught.value.status == status
    assert len(tar_requests(afdb)) == 1
    assert afdb.sleeps == []


def test_a_connection_cut_mid_body_is_retried(afdb: Afdb) -> None:
    afdb.server.truncate.add(afdb.tar_path(PROTEOME_TAR))
    assert af.fetch("P69905") == afdb.text("AF-P69905-F1")
    assert len(tar_requests(afdb)) == 2


def test_a_timeout_is_retried(afdb: Afdb, monkeypatch: pytest.MonkeyPatch) -> None:
    real = urllib.request.urlopen
    calls: list[float] = []

    def flaky(request: urllib.request.Request, timeout: float = 0.0) -> object:
        calls.append(timeout)
        if len(calls) < 3:
            raise TimeoutError("timed out")
        return real(request, timeout=timeout)

    monkeypatch.setattr(urllib.request, "urlopen", flaky)
    assert af.fetch("P69905", timeout=7.5) == afdb.text("AF-P69905-F1")
    assert calls == [7.5, 7.5, 7.5]
    assert len(afdb.sleeps) == 2


def test_a_refused_connection_is_retried_then_reported(afdb: Afdb) -> None:
    dead = afdb.root / "dead.sqlite"
    afdb.write_index(dead)
    conn = sqlite3.connect(dead)
    conn.execute("UPDATE tars SET url = 'http://127.0.0.1:1/v6/' || name")
    conn.execute("UPDATE meta SET value = 'http://127.0.0.1:1/v6' WHERE key = 'base_url'")
    conn.commit()
    conn.close()
    afdb.use(dead)
    with pytest.raises(af.FetchError, match="after 4 attempts"):
        af.fetch("P69905")
    assert len(afdb.sleeps) == 3


# -- fetch_many ----------------------------------------------------------------------------------------


def test_fetch_many_all_good(afdb: Afdb) -> None:
    ids = ["P69905", "AF-Q00001-F2", "p12345-2", "P11111", "P69905"]
    result = af.fetch_many(ids)
    assert list(result) == ["P69905", "AF-Q00001-F2", "p12345-2", "P11111"]  # as given, in order, deduplicated
    assert result["P69905"] == afdb.text("AF-P69905-F1")
    assert result["AF-Q00001-F2"] == afdb.text("AF-Q00001-F2")
    assert result["p12345-2"] == afdb.text("AF-P12345-2-F1")
    assert result["P11111"] == afdb.text("AF-P11111-F1", tar=SWISSPROT_TAR)
    assert len(afdb.server.requests) == 4
    assert af.fetch_many([]) == {}
    assert af.fetch_many("P69905", fmt="pdb") == {"P69905": afdb.text("AF-P69905-F1", "pdb")}


def test_fetch_many_reports_failures_without_losing_successes(afdb: Afdb) -> None:
    ids = ["P69905", "A0A999ZZZ9", "P12345", "Q00001", "not an id", "P11111"]
    with pytest.raises(af.FetchManyError) as caught:
        af.fetch_many(ids, workers=3)
    error = caught.value
    assert error.results == {
        "P69905": afdb.text("AF-P69905-F1"),
        "P12345": afdb.text("AF-P12345-F1"),
        "P11111": afdb.text("AF-P11111-F1", tar=SWISSPROT_TAR),
    }
    assert set(error.errors) == {"A0A999ZZZ9", "Q00001", "not an id"}
    assert isinstance(error.errors["A0A999ZZZ9"], af.NotInBulkSetError)
    assert isinstance(error.errors["Q00001"], af.AmbiguousEntryError)
    assert isinstance(error.errors["not an id"], af.InvalidIdError)
    assert "3 of 6 ids failed" in str(error) and "A0A999ZZZ9" in str(error)
    assert isinstance(error, af.AlphaFoldError)


def test_fetch_many_skip_returns_only_the_successes(afdb: Afdb) -> None:
    afdb.server.etags[afdb.tar_path(SWISSPROT_TAR)] = '"changed"'
    result = af.fetch_many(["P69905", "A0A999ZZZ9", "P11111"], on_error="skip")
    assert result == {"P69905": afdb.text("AF-P69905-F1")}
    with pytest.raises(af.FetchManyError) as caught:
        af.fetch_many(["P69905", "P11111"])
    assert isinstance(caught.value.errors["P11111"], af.StaleIndexError)
    assert list(caught.value.results) == ["P69905"]


def test_fetch_many_argument_checks(afdb: Afdb) -> None:
    with pytest.raises(ValueError, match="on_error"):
        af.fetch_many(["P69905"], on_error="ignore")
    with pytest.raises(ValueError, match="'cif' or 'pdb'"):
        af.fetch_many(["P69905"], fmt="xyz")
    assert afdb.server.requests == []


def test_fetch_many_runs_in_parallel_up_to_workers(afdb: Afdb) -> None:
    afdb.server.delay = 0.2
    ids = ["P69905", "P12345", "P12345-2", "AF-Q00001-F1", "AF-Q00001-F2", "AF-Q00001-F10", "A0A000ISO1-2", "A0A000ISO1-3"]
    started = time.monotonic()
    assert len(af.fetch_many(ids)) == 8
    parallel = time.monotonic() - started
    assert parallel < 1.0, f"8 requests of 0.2 s each took {parallel:.2f} s with 8 workers"
    started = time.monotonic()
    assert len(af.fetch_many(ids[:4], workers=1)) == 4
    assert time.monotonic() - started >= 0.8


# -- download --------------------------------------------------------------------------------------------


def test_download_writes_the_decompressed_file(afdb: Afdb, tmp_path: Path) -> None:
    out = tmp_path / "structures" / "nested"
    path = af.download("p69905", out)
    assert path == out / "AF-P69905-F1-model_v6.cif"
    assert path.read_text() == afdb.text("AF-P69905-F1")
    assert af.download("P69905", out, fmt="pdb") == out / "AF-P69905-F1-model_v6.pdb"
    assert af.download("P12345-2", out).name == "AF-P12345-2-F1-model_v6.cif"
    assert af.download("AF-Q00001-F10", out).name == "AF-Q00001-F10-model_v6.cif"
    assert sorted(p.name for p in out.iterdir()) == [
        "AF-P12345-2-F1-model_v6.cif",
        "AF-P69905-F1-model_v6.cif",
        "AF-P69905-F1-model_v6.pdb",
        "AF-Q00001-F10-model_v6.cif",
    ]
    if sys.platform != "win32":
        assert oct(path.stat().st_mode & 0o777) == "0o644"


def test_download_compressed(afdb: Afdb, tmp_path: Path) -> None:
    path = af.download("P69905", tmp_path, decompress=False)
    assert path == tmp_path / "AF-P69905-F1-model_v6.cif.gz"
    assert path.read_bytes() == af.fetch_bytes("P69905")
    assert gzip.decompress(path.read_bytes()).decode() == afdb.text("AF-P69905-F1")


def test_download_default_directory_is_the_current_one(afdb: Afdb, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    assert af.download("P69905") == Path("AF-P69905-F1-model_v6.cif")
    assert (tmp_path / "AF-P69905-F1-model_v6.cif").is_file()


def test_download_is_atomic(afdb: Afdb, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    target = tmp_path / "AF-P69905-F1-model_v6.cif"
    target.write_text("the previous complete file")
    seen: list[tuple[str, str, str]] = []
    real_replace = os.replace

    def watching_replace(src: str, dst: str) -> None:
        # At the moment of the rename: the temp file is complete, the target is still the old one.
        seen.append((Path(src).read_text(), Path(dst).read_text(), Path(src).parent.name))
        real_replace(src, dst)

    monkeypatch.setattr(os, "replace", watching_replace)
    af.download("P69905", tmp_path)
    assert seen == [(afdb.text("AF-P69905-F1"), "the previous complete file", tmp_path.name)]
    assert target.read_text() == afdb.text("AF-P69905-F1")
    assert [p.name for p in tmp_path.iterdir() if p.name.startswith(".")] == []


def test_a_failed_download_leaves_the_old_file_and_no_temp(afdb: Afdb, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    out = tmp_path / "out"
    out.mkdir()
    target = out / "AF-P69905-F1-model_v6.cif"
    target.write_text("the previous complete file")

    def failing_replace(src: str, dst: str) -> None:
        raise OSError("disk went away")

    with monkeypatch.context() as patched:
        patched.setattr(os, "replace", failing_replace)
        with pytest.raises(OSError, match="disk went away"):
            af.download("P69905", out)
    assert [p.name for p in out.iterdir()] == [target.name]
    assert target.read_text() == "the previous complete file"


def test_a_failed_fetch_writes_nothing(afdb: Afdb, tmp_path: Path) -> None:
    out = tmp_path / "out"
    with pytest.raises(af.NotInBulkSetError):
        af.download("A0A999ZZZ9", out)
    afdb.server.fail_first[afdb.tar_path(PROTEOME_TAR)] = [403]
    with pytest.raises(af.FetchError):
        af.download("P69905", out)
    assert not out.exists() or list(out.iterdir()) == []


# -- live: the real EBI server -------------------------------------------------------------------------------


@pytest.mark.live
def test_live_one_range_request_against_ebi(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Index the first member of a small real tar from its first 512 bytes, then fetch it."""
    name = "UP000000429_85962_HELPY_v6.tar"
    url = f"https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/{name}"
    request = urllib.request.Request(url, headers={"Range": "bytes=0-511", "User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=60) as response:
        assert response.status == 206
        header = response.read()
        total = int(response.headers["Content-Range"].rsplit("/", 1)[1])
        etag, modified = response.headers["ETag"], response.headers["Last-Modified"]
    info = tarfile.TarInfo.frombuf(header, "utf-8", "surrogateescape")
    entry_id, _, rest = info.name.partition("-model_v6.")
    fmt = rest.split(".")[0]
    assert entry_id.startswith("AF-") and fmt in ("cif", "pdb") and rest.endswith(".gz")

    path = tmp_path / "index-v6.sqlite"
    conn = sqlite3.connect(path)
    from conftest import SCHEMA

    conn.executescript(SCHEMA)
    conn.execute(
        "INSERT INTO meta VALUES ('schema_version', '1'), ('afdb_version', '6'), (?, ?)",
        ("base_url", url.rsplit("/", 1)[0]),
    )
    conn.execute(
        "INSERT INTO tars VALUES (1, ?, ?, ?, ?, ?, 'proteome', 'UP000000429', 85962, NULL, NULL, 1)",
        (name, url, total, modified, etag),
    )
    conn.execute("INSERT INTO members VALUES (?, ?, 1, 512, ?)", (entry_id, fmt, info.size))
    conn.commit()
    conn.close()
    monkeypatch.setenv(_index.ENV_INDEX, str(path))
    monkeypatch.setenv(_index.ENV_CACHE, str(tmp_path / "cache"))
    _index.reset()
    try:
        text = af.fetch(entry_id, fmt=fmt)
    finally:
        _index.reset()
    if fmt == "cif":
        assert text.startswith(f"data_{entry_id}")
    else:
        assert "\nATOM " in text
    assert _fetch.DEFAULT_WORKERS == 4


# -- review fixes ----------------------------------------------------------------------------------------------


def test_file_name_id_reads_its_own_format(afdb: Afdb, tmp_path: Path) -> None:
    assert af.fetch("AF-P69905-F1-model_v6.pdb.gz") == afdb.text("AF-P69905-F1", "pdb")
    assert af.fetch("AF-P69905-F1-model_v6.pdb", fmt="pdb") == afdb.text("AF-P69905-F1", "pdb")
    assert gzip.decompress(af.fetch_bytes("af-p69905-f1-model_v6.pdb.gz")).decode() == afdb.text("AF-P69905-F1", "pdb")
    assert af.download("AF-P69905-F1-model_v6.pdb.gz", tmp_path).name == "AF-P69905-F1-model_v6.pdb"
    assert af.fetch_many(["AF-P69905-F1-model_v6.pdb.gz", "P12345"]) == {
        "AF-P69905-F1-model_v6.pdb.gz": afdb.text("AF-P69905-F1", "pdb"),
        "P12345": afdb.text("AF-P12345-F1"),
    }


def test_file_name_id_conflicts_make_no_request(afdb: Afdb) -> None:
    with pytest.raises(af.InvalidIdError, match="is a pdb file name but fmt='cif'"):
        af.fetch("AF-P69905-F1-model_v6.pdb.gz", fmt="cif")
    with pytest.raises(af.InvalidIdError, match="names AlphaFold DB v4, but the index describes v6"):
        af.fetch("AF-P69905-F1-model_v4.cif")
    with pytest.raises(af.FetchManyError) as caught:
        af.fetch_many(["AF-P69905-F1-model_v6.pdb.gz", "P12345"], fmt="cif")
    assert isinstance(caught.value.errors["AF-P69905-F1-model_v6.pdb.gz"], af.InvalidIdError)
    assert list(caught.value.results) == ["P12345"]
    assert len(afdb.server.requests) == 1


def test_download_name_cannot_be_steered_by_the_index(afdb: Afdb, tmp_path: Path) -> None:
    out = tmp_path / "a" / "b" / "c" / "out"
    afdb.use(afdb.edit_index("traversal.sqlite", "UPDATE meta SET value = '6/../../../escaped' WHERE key = 'afdb_version'"))
    with pytest.raises(af.IndexSchemaError, match="afdb_version"):
        af.download("P69905", out)
    assert not out.exists()
    assert [p for p in tmp_path.rglob("*escaped*")] == []
    assert afdb.server.requests == []


def test_hostile_tar_url_is_refused_for_every_reader(afdb: Afdb, tmp_path: Path) -> None:
    afdb.server.files["/internal/admin"] = afdb.tars[PROTEOME_TAR]
    afdb.use(afdb.edit_index("hostile.sqlite", f"UPDATE tars SET url = '{afdb.server.base}/internal/admin'"))
    for call in (
        lambda: af.fetch("P69905"),
        lambda: af.fetch_bytes("P69905"),
        lambda: af.download("P69905", tmp_path / "out"),
    ):
        with pytest.raises(af.IndexSchemaError, match="outside its own base_url"):
            call()
    with pytest.raises(af.FetchManyError) as caught:
        af.fetch_many(["P69905"])
    assert isinstance(caught.value.errors["P69905"], af.IndexSchemaError)
    assert afdb.server.requests == []
    assert not (tmp_path / "out").exists()


# --- the gzip header names the file it holds ------------------------------------------


def _gz_named(name: str, payload: bytes = b"data_x\n") -> bytes:
    import io

    buf = io.BytesIO()
    with gzip.GzipFile(filename=name, mode="wb", fileobj=buf, mtime=0) as fh:
        fh.write(payload)
    return buf.getvalue()


def _member(entry_id: str = "AF-P69905-F1", fmt: str = "cif") -> af.Member:
    return af.Member(entry_id=entry_id, fmt=fmt, tar="t.tar", url="https://x/t.tar", offset=512, length=1)


def test_stored_name_is_read_from_the_gzip_header() -> None:
    assert _fetch._stored_name(_gz_named("AF-P69905-F1-model_v6.cif")) == "AF-P69905-F1-model_v6.cif"
    assert _fetch._stored_name(gzip.compress(b"x", mtime=0)) is None
    assert _fetch._stored_name(b"") is None
    assert _fetch._stored_name(b"\x1f\x8b\x08\x08" + b"\x00" * 6 + b"no-terminator") is None


def test_matching_stored_name_is_accepted() -> None:
    data = _gz_named("AF-P69905-F1-model_v6.cif", b"data_AF-P69905-F1\n")
    assert _fetch._gunzip(_member(), data) == b"data_AF-P69905-F1\n"


@pytest.mark.parametrize(
    "stored",
    ["AF-P68871-F1-model_v6.cif", "AF-P69905-F1-model_v6.pdb", "AF-P69905-F10-model_v6.cif", "AF-P69905-2-F1-model_v6.cif"],
)
def test_another_members_bytes_raise_stale_instead_of_returning_it(stored: str) -> None:
    with pytest.raises(af.StaleIndexError, match="refresh_index"):
        _fetch._gunzip(_member(), _gz_named(stored))
