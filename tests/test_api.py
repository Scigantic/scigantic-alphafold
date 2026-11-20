"""Tests for the two EBI passthroughs in ``_api.py``.

The offline tests talk to an ``http.server`` thread on localhost. Tests marked
``live`` talk to EBI and run only when ``SCIGANTIC_ALPHAFOLD_LIVE=1``.
"""

from __future__ import annotations

import json
import os
import threading
import time
from collections.abc import Callable, Iterator
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any

import pytest

from scigantic_alphafold import _api
from scigantic_alphafold._version import USER_AGENT

Reply = tuple[int, bytes, float]  # status, body, seconds to wait before answering

_skip_unless_live = pytest.mark.skipif(
    os.environ.get("SCIGANTIC_ALPHAFOLD_LIVE") != "1",
    reason="set SCIGANTIC_ALPHAFOLD_LIVE=1 to run tests against EBI",
)


def live(test: Callable[..., None]) -> Callable[..., None]:
    marked: Callable[..., None] = pytest.mark.live(_skip_unless_live(test))
    return marked


HBA = {
    "entryId": "AF-P69905-F1",
    "uniprotAccession": "P69905",
    "gene": "HBA1",
    "taxId": 9606,
    "globalMetricValue": 98.06,
    "latestVersion": 6,
}

DOWNLOAD_METADATA = [
    {
        "archive_name": "UP000000625_83333_ECOLI_v6.tar",
        "species": "Escherichia coli",
        "common_name": "E. coli",
        "latin_common_name": True,
        "reference_proteome": "UP000000625",
        "num_predicted_structures": 4370,
        "size_bytes": 477871104,
        "type": "proteome",
    },
    {
        # A row without reference_proteome: the name has to supply it.
        "archive_name": "UP000000429_85962_HELPY_v6.tar",
        "num_predicted_structures": 1540,
        "size_bytes": 173476864,
        "type": "global_health",
    },
    {
        "archive_name": "UP000094526_86049_9EURO1_v6.tar",
        "reference_proteome": "UP000094526",
        "num_predicted_structures": 11170,
        "size_bytes": 1812907520,
        "type": "global_health",
    },
    {
        "archive_name": "swissprot_cif_v6.tar",
        "label": "Swiss-Prot (CIF files)",
        "num_predicted_structures": 550122,
        "size_bytes": 40062921728,
        "type": "swissprot",
    },
]


class FakeEbi:
    """A local server that replays scripted replies and records what it was asked."""

    def __init__(self) -> None:
        self.replies: dict[str, list[Reply]] = {}
        self.requests: list[tuple[str, dict[str, str]]] = []
        fake = self

        class Handler(BaseHTTPRequestHandler):
            def do_GET(self) -> None:
                fake.requests.append((self.path, {k.lower(): v for k, v in self.headers.items()}))
                queue = fake.replies.get(self.path, [])
                status, body, delay = (
                    (queue.pop(0) if len(queue) > 1 else queue[0])
                    if queue
                    else (404, b'{"detail":"Not Found"}', 0.0)
                )
                if delay:
                    time.sleep(delay)
                try:
                    self.send_response(status)
                    self.send_header("Content-Type", "application/json")
                    self.send_header("Content-Length", str(len(body)))
                    self.end_headers()
                    self.wfile.write(body)
                except OSError:
                    pass  # the client gave up waiting, which the timeout test expects

            def log_message(self, format: str, *args: Any) -> None:
                pass

        self.server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self.server.daemon_threads = True
        self.base = f"http://127.0.0.1:{self.server.server_address[1]}"
        self.thread = threading.Thread(
            target=self.server.serve_forever, kwargs={"poll_interval": 0.01}, daemon=True
        )
        self.thread.start()

    def script(self, path: str, *replies: tuple[int, Any] | tuple[int, Any, float]) -> None:
        """Queue replies for ``path``. The last one repeats once the queue runs down."""
        queue: list[Reply] = []
        for reply in replies:
            payload = reply[1]
            body = payload if isinstance(payload, bytes) else json.dumps(payload).encode()
            queue.append((reply[0], body, reply[2] if len(reply) == 3 else 0.0))
        self.replies[path] = queue

    def paths(self) -> list[str]:
        return [path for path, _ in self.requests]

    def close(self) -> None:
        self.server.shutdown()
        self.server.server_close()


@pytest.fixture()
def ebi(monkeypatch: pytest.MonkeyPatch) -> Iterator[FakeEbi]:
    fake = FakeEbi()
    monkeypatch.setattr(_api, "API_BASE_URL", fake.base + "/api")
    monkeypatch.setattr(_api, "FTP_BASE_URL", fake.base + "/pub/databases/alphafold")
    monkeypatch.setattr(_api, "_BACKOFF_SECONDS", 0.0)
    yield fake
    fake.close()


# --- metadata() ---------------------------------------------------------------


def test_metadata_success(ebi: FakeEbi) -> None:
    ebi.script("/api/prediction/P69905", (200, [HBA]))
    assert _api.metadata("P69905") == [HBA]
    assert ebi.paths() == ["/api/prediction/P69905"]


def test_metadata_sends_the_honest_user_agent(ebi: FakeEbi) -> None:
    ebi.script("/api/prediction/P69905", (200, [HBA]))
    _api.metadata("P69905")
    headers = ebi.requests[0][1]
    assert headers["user-agent"] == USER_AGENT
    assert USER_AGENT.startswith("scigantic-alphafold/")
    assert "Mozilla" not in USER_AGENT


def test_metadata_accepts_a_model_id_and_strips_whitespace(ebi: FakeEbi) -> None:
    ebi.script("/api/prediction/AF-P69905-F1", (200, [HBA]))
    assert _api.metadata(" AF-P69905-F1\n") == [HBA]
    assert ebi.paths() == ["/api/prediction/AF-P69905-F1"]


def test_metadata_missing_accession_is_an_empty_list(ebi: FakeEbi) -> None:
    # What EBI sends for a well-formed accession with no prediction: 404 and {}.
    ebi.script("/api/prediction/Z9Z9Z9", (404, {}))
    assert _api.metadata("Z9Z9Z9") == []
    assert len(ebi.requests) == 1


def test_metadata_route_not_found_is_an_error_not_an_empty_list(ebi: FakeEbi) -> None:
    # What EBI sends for a path it does not route. If the API ever moves, this
    # must not read as "no prediction".
    ebi.script("/api/prediction/P69905", (404, {"detail": "Not Found"}))
    with pytest.raises(_api.ApiError) as caught:
        _api.metadata("P69905")
    assert caught.value.status == 404
    assert "Not Found" in str(caught.value)
    assert len(ebi.requests) == 1


def test_metadata_400_is_not_retried_and_carries_the_server_message(ebi: FakeEbi) -> None:
    message = "Invalid identifier format. Please use a UniProt accession or a supported AlphaFold DB ID."
    ebi.script("/api/prediction/A0A000000000", (400, {"error": message}))
    with pytest.raises(_api.ApiError) as caught:
        _api.metadata("A0A000000000")
    assert caught.value.status == 400
    assert message in str(caught.value)
    assert caught.value.url.endswith("/api/prediction/A0A000000000")
    assert len(ebi.requests) == 1


def test_metadata_429_is_not_retried(ebi: FakeEbi) -> None:
    ebi.script("/api/prediction/P69905", (429, b"slow down"))
    with pytest.raises(_api.ApiError) as caught:
        _api.metadata("P69905")
    assert caught.value.status == 429
    assert len(ebi.requests) == 1


def test_metadata_5xx_then_success(ebi: FakeEbi) -> None:
    ebi.script("/api/prediction/P69905", (503, b"busy"), (502, b"busy"), (200, [HBA]))
    assert _api.metadata("P69905") == [HBA]
    assert len(ebi.requests) == 3


def test_metadata_backs_off_between_attempts(ebi: FakeEbi, monkeypatch: pytest.MonkeyPatch) -> None:
    waits: list[float] = []
    monkeypatch.setattr(_api, "_BACKOFF_SECONDS", 1.0)
    monkeypatch.setattr(time, "sleep", waits.append)
    ebi.script("/api/prediction/P69905", (500, b""), (500, b""), (200, [HBA]))
    assert _api.metadata("P69905") == [HBA]
    assert waits == [1.0, 2.0]


def test_metadata_gives_up_after_three_5xx(ebi: FakeEbi) -> None:
    ebi.script("/api/prediction/P69905", (500, b"down"))
    with pytest.raises(_api.ApiError) as caught:
        _api.metadata("P69905")
    assert caught.value.status == 500
    assert "3 attempts" in str(caught.value)
    assert len(ebi.requests) == 3


def test_metadata_timeout_is_retried_then_raised(ebi: FakeEbi) -> None:
    ebi.script("/api/prediction/P69905", (200, [HBA], 1.0))
    started = time.monotonic()
    with pytest.raises(_api.ApiError) as caught:
        _api.metadata("P69905", timeout=0.1)
    assert caught.value.status is None
    assert "no response within 0.1 s" in str(caught.value)
    assert len(ebi.requests) == 3
    assert time.monotonic() - started < 1.0


def test_metadata_timeout_then_success(ebi: FakeEbi) -> None:
    ebi.script("/api/prediction/P69905", (200, [HBA], 1.0), (200, [HBA]))
    assert _api.metadata("P69905", timeout=0.2) == [HBA]
    assert len(ebi.requests) == 2


def test_metadata_connection_refused_is_not_retried(monkeypatch: pytest.MonkeyPatch) -> None:
    fake = FakeEbi()
    fake.close()  # the port is now closed
    waits: list[float] = []
    monkeypatch.setattr(_api, "API_BASE_URL", fake.base + "/api")
    monkeypatch.setattr(time, "sleep", waits.append)
    # Windows retries the SYN to a closed loopback port for about 2 s before it reports the
    # refusal, so a 2 s timeout races it and the call is retried as a timeout. Leave room.
    with pytest.raises(_api.ApiError) as caught:
        _api.metadata("P69905", timeout=30.0)
    assert caught.value.status is None
    assert waits == []


@pytest.mark.parametrize(
    "qualifier",
    [
        "P69905/../../files",
        "../openapi.json",
        "P69905?include_complexes=true",
        "P69905#x",
        "P69905%2F..",
        "P69905 Q9Y6K9",
        "P69905\\x",
        "P69905;x",
        "@evil.example/P69905",
        "..",
        ".",
        "-P69905",
        "",
        "   ",
        "P" * 65,
        "P69905é",
    ],
)
def test_metadata_rejects_path_characters_without_a_request(ebi: FakeEbi, qualifier: str) -> None:
    with pytest.raises(ValueError):
        _api.metadata(qualifier)
    assert ebi.requests == []


def test_metadata_rejects_non_strings(ebi: FakeEbi) -> None:
    with pytest.raises(TypeError):
        _api.metadata(69905)  # type: ignore[arg-type]
    assert ebi.requests == []


@pytest.mark.parametrize("body", [b"<html>maintenance</html>", b'{"a": 1}', b"[1, 2]", b'"x"'])
def test_metadata_unexpected_200_body_is_an_error(ebi: FakeEbi, body: bytes) -> None:
    ebi.script("/api/prediction/P69905", (200, body))
    with pytest.raises(_api.ApiError):
        _api.metadata("P69905")
    assert len(ebi.requests) == 1


# --- proteomes() --------------------------------------------------------------

DM_PATH = "/pub/databases/alphafold/download_metadata.json"


def test_proteomes_adds_url_tax_id_and_reference_proteome(ebi: FakeEbi) -> None:
    ebi.script(DM_PATH, (200, DOWNLOAD_METADATA))
    rows = _api.proteomes()
    assert [r["archive_name"] for r in rows] == [r["archive_name"] for r in DOWNLOAD_METADATA]
    ecoli, helpy, euro, swissprot = rows

    assert ecoli["url"] == ebi.base + "/pub/databases/alphafold/v6/UP000000625_83333_ECOLI_v6.tar"
    assert ecoli["tax_id"] == 83333
    assert ecoli["reference_proteome"] == "UP000000625"
    assert ecoli["species"] == "Escherichia coli"
    assert ecoli["size_bytes"] == 477871104

    assert helpy["tax_id"] == 85962
    assert helpy["reference_proteome"] == "UP000000429"  # parsed from the name

    assert euro["tax_id"] == 86049  # the mnemonic 9EURO1 starts with a digit
    assert euro["reference_proteome"] == "UP000094526"

    assert swissprot["url"].endswith("/v6/swissprot_cif_v6.tar")
    assert swissprot["tax_id"] is None
    assert swissprot["reference_proteome"] is None
    assert swissprot["label"] == "Swiss-Prot (CIF files)"

    assert ebi.requests[0][1]["user-agent"] == USER_AGENT


def test_proteomes_keeps_values_the_server_already_has(ebi: FakeEbi) -> None:
    row = {"archive_name": "UP000000625_83333_ECOLI_v6.tar", "tax_id": 562, "reference_proteome": "UPX"}
    ebi.script(DM_PATH, (200, [row]))
    (out,) = _api.proteomes()
    assert out["tax_id"] == 562
    assert out["reference_proteome"] == "UPX"


def test_proteomes_url_follows_the_version_in_the_name(ebi: FakeEbi) -> None:
    ebi.script(DM_PATH, (200, [{"archive_name": "UP000000625_83333_ECOLI_v7.tar"}]))
    (out,) = _api.proteomes()
    assert out["url"].endswith("/v7/UP000000625_83333_ECOLI_v7.tar")


def test_proteomes_does_not_mutate_between_calls(ebi: FakeEbi) -> None:
    ebi.script(DM_PATH, (200, DOWNLOAD_METADATA))
    first = _api.proteomes()
    first[0]["url"] = "changed"
    assert _api.proteomes()[0]["url"].endswith("_ECOLI_v6.tar")


@pytest.mark.parametrize(
    "name", ["../../etc/passwd_v6.tar", "a/b_v6.tar", "ECOLI.tar", "ECOLI_v6.tar.gz", "", None, 7]
)
def test_proteomes_rejects_archive_names_that_are_not_tar_file_names(ebi: FakeEbi, name: Any) -> None:
    ebi.script(DM_PATH, (200, [{"archive_name": name}]))
    with pytest.raises(_api.ApiError):
        _api.proteomes()


def test_proteomes_5xx_then_success(ebi: FakeEbi) -> None:
    ebi.script(DM_PATH, (503, b""), (200, DOWNLOAD_METADATA))
    assert len(_api.proteomes()) == 4
    assert len(ebi.requests) == 2


def test_proteomes_404_is_an_error_and_not_retried(ebi: FakeEbi) -> None:
    with pytest.raises(_api.ApiError) as caught:
        _api.proteomes()
    assert caught.value.status == 404
    assert len(ebi.requests) == 1


def test_proteomes_unexpected_body_is_an_error(ebi: FakeEbi) -> None:
    ebi.script(DM_PATH, (200, {"archives": []}))
    with pytest.raises(_api.ApiError):
        _api.proteomes()


# --- live ---------------------------------------------------------------------


@live
def test_live_metadata_p69905() -> None:
    rows = _api.metadata("P69905")
    assert len(rows) == 1
    row = rows[0]
    assert row["entryId"] == "AF-P69905-F1"
    assert row["uniprotAccession"] == "P69905"
    assert row["gene"] == "HBA1"
    assert row["taxId"] == 9606
    assert row["latestVersion"] == 6
    assert len(row["sequence"]) == 142
    assert _api.metadata("AF-P69905-F1") == rows


@live
def test_live_metadata_missing_accession() -> None:
    assert _api.metadata("Z9Z9Z9") == []


@live
def test_live_metadata_malformed_identifier() -> None:
    with pytest.raises(_api.ApiError) as caught:
        _api.metadata("A0A000000000")
    assert caught.value.status == 400
    assert "Invalid identifier format" in str(caught.value)


@live
def test_live_proteomes_totals() -> None:
    rows = _api.proteomes()
    assert len(rows) == 48
    assert all(r["archive_name"].endswith("_v6.tar") for r in rows)
    assert all(r["url"] == f"{_api.FTP_BASE_URL}/v6/{r['archive_name']}" for r in rows)
    by_type: dict[str, int] = {}
    for r in rows:
        by_type[r["type"]] = by_type.get(r["type"], 0) + 1
    assert by_type == {"proteome": 16, "global_health": 30, "swissprot": 2}
    assert 153.0e9 < sum(r["size_bytes"] for r in rows) < 153.5e9
    for r in rows:
        if r["type"] == "swissprot":
            assert r["tax_id"] is None and r["reference_proteome"] is None
            assert r["num_predicted_structures"] == 550122
        else:
            assert isinstance(r["tax_id"], int)
            assert r["archive_name"].startswith(f"{r['reference_proteome']}_{r['tax_id']}_")
    ecoli = next(r for r in rows if r["archive_name"] == "UP000000625_83333_ECOLI_v6.tar")
    assert ecoli["tax_id"] == 83333
    assert ecoli["num_predicted_structures"] == 4370


@live
def test_live_proteome_url_exists() -> None:
    import urllib.request

    helpy = next(r for r in _api.proteomes() if "_HELPY_" in r["archive_name"])
    request = urllib.request.Request(helpy["url"], method="HEAD", headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=30) as response:
        assert response.status == 200
        assert int(response.headers["Content-Length"]) == helpy["size_bytes"]
