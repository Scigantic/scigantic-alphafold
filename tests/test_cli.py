from __future__ import annotations

import os
import subprocess
import sys
import types
from pathlib import Path
from typing import Any

import pytest

import scigantic_alphafold as af
from scigantic_alphafold import _index, cli

from conftest import PROTEOME_TAR, SWISSPROT_TAR, Afdb


def run(capsys: pytest.CaptureFixture[str], *argv: str) -> tuple[int, str, str]:
    code = cli.main(list(argv))
    captured = capsys.readouterr()
    return code, captured.out, captured.err


def one_line_error(err: str) -> str:
    lines = err.splitlines()
    assert len(lines) == 1, err
    assert lines[0].startswith("scigantic-alphafold: error: ")
    assert "Traceback" not in err
    return lines[0]


# -- get ----------------------------------------------------------------------------------------------


def test_get_writes_files_and_prints_their_paths(afdb: Afdb, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    out = tmp_path / "out"
    code, stdout, stderr = run(capsys, "get", "P69905", "p12345-2", "-o", str(out))
    assert (code, stderr) == (0, "")
    assert stdout.splitlines() == [str(out / "AF-P69905-F1-model_v6.cif"), str(out / "AF-P12345-2-F1-model_v6.cif")]
    assert (out / "AF-P69905-F1-model_v6.cif").read_text() == afdb.text("AF-P69905-F1")


def test_get_pdb_into_the_current_directory(afdb: Afdb, tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    code, stdout, _ = run(capsys, "get", "P69905", "-f", "pdb")
    assert code == 0
    assert stdout.strip() == "AF-P69905-F1-model_v6.pdb"
    assert (tmp_path / "AF-P69905-F1-model_v6.pdb").read_text() == afdb.text("AF-P69905-F1", "pdb")


def test_get_stdout_writes_the_decompressed_text_and_no_file(afdb: Afdb, tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    code, stdout, stderr = run(capsys, "get", "P69905", "--stdout")
    assert (code, stderr) == (0, "")
    assert stdout == afdb.text("AF-P69905-F1")
    assert list(tmp_path.glob("AF-*")) == []
    code, stdout, _ = run(capsys, "get", "P69905", "AF-Q00001-F2", "--stdout", "--format", "cif")
    assert stdout == afdb.text("AF-P69905-F1") + afdb.text("AF-Q00001-F2")


def test_get_missing_id_exits_1_with_one_line(afdb: Afdb, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    code, stdout, stderr = run(capsys, "get", "A0A999ZZZ9", "-o", str(tmp_path / "out"))
    assert code == 1 and stdout == ""
    line = one_line_error(stderr)
    assert "A0A999ZZZ9" in line and "not in any AlphaFold DB v6 bulk tar" in line


def test_get_keeps_going_after_a_failure(afdb: Afdb, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    out = tmp_path / "out"
    code, stdout, stderr = run(capsys, "get", "A0A999ZZZ9", "P69905", "Q00001", "-o", str(out))
    assert code == 1
    assert stdout.splitlines() == [str(out / "AF-P69905-F1-model_v6.cif")]
    lines = stderr.splitlines()
    assert len(lines) == 2 and all(line.startswith("scigantic-alphafold: error: ") for line in lines)
    assert "A0A999ZZZ9" in lines[0] and "AF-Q00001-F2" in lines[1]


@pytest.mark.parametrize(
    ("break_it", "expect"),
    [("etag", "refresh"), ("403", "HTTP 403"), ("200", "HTTP 200 instead of 206")],
)
def test_get_server_problems_exit_1(afdb: Afdb, tmp_path: Path, capsys: pytest.CaptureFixture[str], break_it: str, expect: str) -> None:
    path = afdb.tar_path(PROTEOME_TAR)
    if break_it == "etag":
        afdb.server.etags[path] = '"changed"'
    elif break_it == "403":
        afdb.server.fail_first[path] = [403]
    else:
        afdb.server.ignore_range.add(path)
    code, stdout, stderr = run(capsys, "get", "P69905", "-o", str(tmp_path))
    assert code == 1 and stdout == ""
    assert expect in one_line_error(stderr)


def test_get_rejects_an_unknown_format(afdb: Afdb, capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as caught:
        cli.main(["get", "P69905", "-f", "mmtf"])
    assert caught.value.code == 2
    assert "invalid choice" in capsys.readouterr().err


# -- locate ---------------------------------------------------------------------------------------------


def test_locate_prints_every_copy(afdb: Afdb, capsys: pytest.CaptureFixture[str]) -> None:
    code, stdout, stderr = run(capsys, "locate", "P69905")
    assert (code, stderr) == (0, "")
    rows = [line.split("\t") for line in stdout.splitlines()]
    assert rows[0] == ["entry_id", "fmt", "tar", "offset", "length", "url"]
    assert [(r[0], r[1], r[2]) for r in rows[1:]] == [
        ("AF-P69905-F1", "cif", PROTEOME_TAR),
        ("AF-P69905-F1", "cif", SWISSPROT_TAR),
        ("AF-P69905-F1", "pdb", PROTEOME_TAR),
    ]
    offset, length = afdb.where[PROTEOME_TAR][("AF-P69905-F1", "cif")]
    assert rows[1][3:] == [str(offset), str(length), afdb.tar_url(PROTEOME_TAR)]
    assert afdb.server.requests == []  # locate is index only


def test_locate_lists_fragments_and_filters_by_format(afdb: Afdb, capsys: pytest.CaptureFixture[str]) -> None:
    code, stdout, _ = run(capsys, "locate", "Q00001")
    assert code == 0
    assert [line.split("\t")[0] for line in stdout.splitlines()[1:]] == ["AF-Q00001-F1", "AF-Q00001-F2", "AF-Q00001-F10"]
    code, stdout, _ = run(capsys, "locate", "P69905", "-f", "pdb")
    assert [line.split("\t")[1] for line in stdout.splitlines()[1:]] == ["pdb"]


def test_locate_missing_exits_1(afdb: Afdb, capsys: pytest.CaptureFixture[str]) -> None:
    code, stdout, stderr = run(capsys, "locate", "A0A999ZZZ9")
    assert code == 1 and stdout == ""
    assert "not in any AlphaFold DB v6 bulk tar" in one_line_error(stderr)
    code, stdout, stderr = run(capsys, "locate", "P11111", "-f", "pdb")
    assert code == 1
    assert "only as cif, not as pdb" in one_line_error(stderr)
    code, _, stderr = run(capsys, "locate", "not an id")
    assert code == 1
    assert "is not a UniProt accession" in one_line_error(stderr)


# -- info and refresh --------------------------------------------------------------------------------------


def test_info(afdb: Afdb, capsys: pytest.CaptureFixture[str]) -> None:
    code, stdout, stderr = run(capsys, "info")
    assert (code, stderr) == (0, "")
    shown = dict(line.split(": ", 1) for line in stdout.splitlines())
    assert shown["path"] == str(afdb.index_path)
    assert shown["schema_version"] == "1"
    assert shown["afdb_version"] == "6"
    assert shown["tars"] == "2"
    assert shown["members"] == "11"


def test_info_with_no_published_index_exits_1(afdb: Afdb, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(_index, "DEFAULT_INDEX_URL", f"{afdb.server.base}/releases/index-v6.sqlite.gz")
    afdb.use(None)
    code, stdout, stderr = run(capsys, "info")
    assert code == 1 and stdout == ""
    line = one_line_error(stderr)
    assert "may not be published yet or may be mid-update" in line and "SCIGANTIC_ALPHAFOLD_INDEX" in line


def test_info_with_an_unknown_schema_exits_1(afdb: Afdb, capsys: pytest.CaptureFixture[str]) -> None:
    afdb.use(afdb.write_index(afdb.root / "future.sqlite", schema_version="9"))
    code, _, stderr = run(capsys, "info")
    assert code == 1
    assert "schema_version 9" in one_line_error(stderr)


def test_refresh(afdb: Afdb, capsys: pytest.CaptureFixture[str]) -> None:
    path = "/index/index-v6.sqlite.gz"
    afdb.use(afdb.server.base + path)
    code, stdout, stderr = run(capsys, "refresh")
    assert (code, stderr) == (0, "")
    cached = Path(stdout.strip())
    assert cached.parent == afdb.cache_dir and cached.is_file()
    code, stdout, _ = run(capsys, "refresh")
    assert code == 0 and Path(stdout.strip()) == cached
    assert afdb.server.count(path) == 2
    afdb.server.fail_first[path] = [404]
    code, _, stderr = run(capsys, "refresh")
    assert code == 1
    assert "HTTP 404" in one_line_error(stderr)
    assert cached.is_file()  # a failed refresh keeps the copy that worked


# -- build-part and merge delegate to _build -----------------------------------------------------------------


@pytest.fixture
def fake_build(monkeypatch: pytest.MonkeyPatch) -> list[tuple[str, tuple[Any, ...], dict[str, Any]]]:
    calls: list[tuple[str, tuple[Any, ...], dict[str, Any]]] = []
    module = types.ModuleType("scigantic_alphafold._build")

    def build_part(*args: Any, **kwargs: Any) -> Path:
        calls.append(("build_part", args, kwargs))
        if args[0] == "missing_v6.tar":
            raise OSError("HTTP Error 404: Not Found")
        return Path(args[1]) / f"{args[0]}.tsv.gz"

    def merge_parts(*args: Any, **kwargs: Any) -> Path:
        calls.append(("merge_parts", args, kwargs))
        return Path(args[1])

    def _filtered_metadata(only: list[str]) -> list[dict[str, Any]]:
        return [{"archive_name": name} for name in only]

    module._filtered_metadata = _filtered_metadata  # type: ignore[attr-defined]

    module.build_part = build_part  # type: ignore[attr-defined]
    module.merge_parts = merge_parts  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "scigantic_alphafold._build", module)
    return calls


def test_build_part_delegates(fake_build: list[Any], capsys: pytest.CaptureFixture[str]) -> None:
    code, stdout, stderr = run(capsys, "build-part", "UP000000429_85962_HELPY_v6.tar", "-o", "parts")
    assert (code, stderr) == (0, "")
    assert stdout.strip() == str(Path("parts") / "UP000000429_85962_HELPY_v6.tar.tsv.gz")
    assert fake_build == [("build_part", ("UP000000429_85962_HELPY_v6.tar", "parts"), {})]
    fake_build.clear()
    code, _, _ = run(capsys, "build-part", "x_v6.tar", "--base-url", "http://127.0.0.1:9/v6/")
    assert code == 0
    assert fake_build == [("build_part", ("x_v6.tar", "."), {"base_url": "http://127.0.0.1:9/v6/"})]


def test_build_part_failure_exits_1(fake_build: list[Any], capsys: pytest.CaptureFixture[str]) -> None:
    code, stdout, stderr = run(capsys, "build-part", "missing_v6.tar")
    assert code == 1 and stdout == ""
    assert "404" in one_line_error(stderr)


def test_merge_delegates(fake_build: list[Any], capsys: pytest.CaptureFixture[str]) -> None:
    code, stdout, stderr = run(capsys, "merge", "parts", "-o", "index-v6.sqlite")
    assert (code, stderr) == (0, "")
    assert stdout.strip() == "index-v6.sqlite"
    assert fake_build == [("merge_parts", ("parts", "index-v6.sqlite"), {})]


def test_merge_requires_an_output(fake_build: list[Any], capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as caught:
        cli.main(["merge", "parts"])
    assert caught.value.code == 2
    assert fake_build == []


# -- the program itself ----------------------------------------------------------------------------------------


def test_no_command_and_version(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as caught:
        cli.main([])
    assert caught.value.code == 2
    capsys.readouterr()
    with pytest.raises(SystemExit) as caught:
        cli.main(["--version"])
    assert caught.value.code == 0
    assert capsys.readouterr().out.strip() == f"scigantic-alphafold {af.__version__}"


def test_help_lists_every_subcommand(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit):
        cli.main(["--help"])
    text = capsys.readouterr().out
    for command in ("get", "locate", "info", "refresh", "build-part", "merge"):
        assert command in text


def test_module_entry_point_end_to_end(afdb: Afdb, tmp_path: Path) -> None:
    """The real process: exit codes and streams as a shell sees them."""
    env = {**os.environ, _index.ENV_INDEX: str(afdb.index_path), _index.ENV_CACHE: str(afdb.cache_dir)}
    base = [sys.executable, "-m", "scigantic_alphafold.cli"]
    good = subprocess.run([*base, "get", "P69905", "--stdout"], capture_output=True, text=True, env=env)
    assert (good.returncode, good.stderr) == (0, "")
    assert good.stdout == afdb.text("AF-P69905-F1")
    bad = subprocess.run([*base, "get", "A0A999ZZZ9", "-o", str(tmp_path)], capture_output=True, text=True, env=env)
    assert bad.returncode == 1 and bad.stdout == ""
    assert bad.stderr.startswith("scigantic-alphafold: error: A0A999ZZZ9") and bad.stderr.count("\n") == 1


# -- review fixes ----------------------------------------------------------------------------------------------

REPO = Path(__file__).resolve().parents[1]


def test_get_honors_a_file_name_id(afdb: Afdb, tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    code, stdout, _ = run(capsys, "get", "AF-P69905-F1-model_v6.pdb.gz", "-o", str(tmp_path))
    assert code == 0
    assert stdout.strip() == str(tmp_path / "AF-P69905-F1-model_v6.pdb")
    code, stdout, stderr = run(capsys, "get", "AF-P69905-F1-model_v6.pdb.gz", "-f", "cif", "-o", str(tmp_path))
    assert code == 1 and stdout == ""
    assert "is a pdb file name but fmt='cif'" in one_line_error(stderr)
    code, _, stderr = run(capsys, "get", "AF-P69905-F1-model_v4.cif", "--stdout")
    assert code == 1
    assert "the index describes v6" in one_line_error(stderr)
    code, stdout, _ = run(capsys, "locate", "AF-P69905-F1-model_v6.pdb.gz")
    assert code == 0
    assert [line.split("\t")[1] for line in stdout.splitlines()[1:]] == ["pdb"]


def test_relative_index_path_from_the_shell(afdb: Afdb) -> None:
    env = {**os.environ, _index.ENV_INDEX: "index-v6.sqlite", _index.ENV_CACHE: os.path.join("rel", "cache")}
    done = subprocess.run(
        [sys.executable, "-m", "scigantic_alphafold.cli", "get", "P69905", "--stdout"],
        capture_output=True, text=True, env=env, cwd=afdb.root,
    )  # fmt: skip
    assert (done.returncode, done.stderr) == (0, "")
    assert done.stdout == afdb.text("AF-P69905-F1")


def test_merge_only_delegates_a_filtered_tar_list(fake_build: list[Any], capsys: pytest.CaptureFixture[str]) -> None:
    code, stdout, stderr = run(capsys, "merge", "parts", "-o", "index-v6.sqlite", "--only", "A_v6.tar", "B_v6.tar")
    assert (code, stderr) == (0, "")
    assert stdout.strip() == "index-v6.sqlite"
    assert fake_build == [
        ("merge_parts", ("parts", "index-v6.sqlite"), {"metadata": [{"archive_name": "A_v6.tar"}, {"archive_name": "B_v6.tar"}]})
    ]


def test_the_partial_build_advice_works_as_written(afdb: Afdb, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch) -> None:
    """build-part one tar, merge --only that tar, point the env var at the result (a relative path), fetch."""
    from conftest import build_tar, structure_text
    from scigantic_alphafold import _build

    name = "UP000000002_2_MINI_v6.tar"
    members = [("AF-M00001-F1", "cif"), ("AF-M00001-F1", "pdb"), ("AF-M00002-F1", "cif"), ("AF-M00002-F1", "pdb")]
    afdb.server.files[f"/v6/{name}"], _ = build_tar(members, name)
    listing = [
        {"archive_name": name, "num_predicted_structures": 2, "size_bytes": len(afdb.server.files[f"/v6/{name}"]),
         "type": "proteome", "species": "Minimus", "common_name": "mini", "reference_proteome": "UP000000002", "label": "mini"},
        {"archive_name": "UP000000003_3_NEVERBUILT_v6.tar", "num_predicted_structures": 9, "size_bytes": 1,
         "type": "proteome", "species": "Other", "common_name": "other", "reference_proteome": "UP000000003", "label": "other"},
    ]  # fmt: skip
    monkeypatch.setattr(_build, "fetch_download_metadata", lambda **kwargs: listing)
    monkeypatch.chdir(afdb.root)

    code, _, stderr = run(capsys, "build-part", name, "-o", "parts", "--base-url", f"{afdb.server.base}/v6")
    assert code == 0, stderr
    # Without --only, every listed tar needs a part file: this is what the old advice ran into.
    code, _, stderr = run(capsys, "merge", "parts", "-o", "full.sqlite")
    assert code == 1 and "scigantic-alphafold: error:" in stderr
    assert not (afdb.root / "full.sqlite").exists()
    code, stdout, stderr = run(capsys, "merge", "parts", "-o", "index-v6-mini.sqlite", "--only", name)
    assert code == 0, stderr

    afdb.use("index-v6-mini.sqlite")
    assert [t.name for t in af.tars()] == [name]
    assert af.fetch("M00002", fmt="pdb") == structure_text("AF-M00002-F1", "pdb", name)
    assert af.fetch("m00001") == structure_text("AF-M00001-F1", "cif", name)


def test_version_has_one_source() -> None:
    from importlib.metadata import version

    assert version("scigantic-alphafold") == af.__version__
    pyproject = REPO / "pyproject.toml"
    if not pyproject.is_file():
        pytest.skip("not running from a source tree")
    text = pyproject.read_text(encoding="utf-8")
    project = text.split("[project]", 1)[1].split("\n[", 1)[0]
    assert 'dynamic = ["version"]' in project
    assert "\nversion =" not in project
    assert 'version = { attr = "scigantic_alphafold._version.__version__" }' in text
    assert 'license = "MIT-0"' in project and "License ::" not in text


def test_sdist_manifest_ships_the_test_fixtures() -> None:
    manifest = REPO / "MANIFEST.in"
    if not (REPO / "pyproject.toml").is_file():
        pytest.skip("not running from a source tree")
    lines = manifest.read_text(encoding="utf-8").splitlines()
    assert "include tests/conftest.py" in lines and "prune build" in lines


def publish_workflow() -> str:
    path = REPO / ".github" / "workflows" / "publish.yml"
    if not path.is_file():
        pytest.skip("the workflows are not part of an sdist")
    return path.read_text(encoding="utf-8")


def test_publish_workflow_is_locked_down() -> None:
    text = publish_workflow()
    assert "\npermissions:\n  contents: read\n" in text
    assert "workflow_dispatch" not in text.split("\njobs:", 1)[0].split("\non:", 1)[1]
    assert "tags: ['v*']" in text
    # The tag check comes before anything is built or uploaded.
    assert text.index("The tag must be v<version>") < text.index("python -m build") < text.index("twine upload")
    assert "needs: [test, typecheck]" in text


@pytest.mark.skipif(sys.platform == "win32", reason="runs the workflow's bash step")
@pytest.mark.parametrize(
    ("ref_type", "ref_name", "ok"),
    [("tag", f"v{af.__version__}", True), ("tag", "v999.0.0", False), ("tag", af.__version__, False), ("branch", f"v{af.__version__}", False), ("branch", "main", False)],
)
def test_publish_tag_check_step(ref_type: str, ref_name: str, ok: bool) -> None:
    import textwrap

    text = publish_workflow()
    step = text.split("- name: The tag must be v<version>", 1)[1].split("run: |\n", 1)[1].split("      - name:", 1)[0]
    script = textwrap.dedent(step)
    env = {**os.environ, "GITHUB_REF_TYPE": ref_type, "GITHUB_REF_NAME": ref_name}
    done = subprocess.run(["bash", "-euo", "pipefail", "-c", script], capture_output=True, text=True, env=env, cwd=REPO)
    assert (done.returncode == 0) == ok, done.stdout + done.stderr
    if not ok:
        assert "nothing was published" in done.stdout


def test_ci_has_a_windows_leg() -> None:
    path = REPO / ".github" / "workflows" / "ci.yml"
    if not path.is_file():
        pytest.skip("the workflows are not part of an sdist")
    text = path.read_text(encoding="utf-8")
    assert "os: windows-latest" in text and "runs-on: ${{ matrix.os }}" in text
