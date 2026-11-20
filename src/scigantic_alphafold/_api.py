"""Thin clients for two public EBI resources that need no index.

``metadata()`` passes through the AlphaFold DB prediction API and ``proteomes()``
passes through the FTP area's ``download_metadata.json``. Neither function reads
structure files; that is the job of the index and the range reader.

Behavior of the prediction API as observed on 2026-11-20:

* a known accession or model id returns 200 and a JSON list of predictions,
* a well-formed accession with no prediction returns 404 and the body ``{}``,
* a malformed identifier returns 400 and ``{"error": "Invalid identifier format. ..."}``,
* a path the API does not route returns 404 and ``{"detail": "Not Found"}``.
"""

from __future__ import annotations

import http.client
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any

from ._errors import AlphaFoldError, InvalidIdError
from ._version import USER_AGENT

__all__ = ["ApiError", "metadata", "proteomes"]

# Read at call time, so tests can point them at a local server.
API_BASE_URL = "https://alphafold.ebi.ac.uk/api"
FTP_BASE_URL = "https://ftp.ebi.ac.uk/pub/databases/alphafold"

_ATTEMPTS = 3  # one request plus two retries, on timeouts and 5xx only
_BACKOFF_SECONDS = 1.0  # doubled after every failed attempt

# UniProt accessions (P69905, A0A385XJ53) and AlphaFold DB model ids
# (AF-P69905-F1, AF-F4HVG8-2-F1, AF-0000000078834380) use only these characters.
_QUALIFIER_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]{0,63}")
_ARCHIVE_NAME_RE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\.tar")
_ARCHIVE_VERSION_RE = re.compile(r"_v(\d+)\.tar")
_PROTEOME_NAME_RE = re.compile(r"(UP\d+)_(\d+)_[A-Za-z0-9]+_v\d+\.tar")


class ApiError(AlphaFoldError, RuntimeError):
    """An EBI endpoint answered with an error, an unexpected body, or not at all.

    ``status`` is the HTTP status code, or ``None`` when no response arrived
    (a timeout or a connection failure). ``url`` is the URL that was requested.
    """

    def __init__(self, message: str, *, url: str, status: int | None = None) -> None:
        super().__init__(message)
        self.url = url
        self.status = status


class _HttpStatusError(Exception):
    def __init__(self, status: int, body: bytes) -> None:
        super().__init__(status)
        self.status = status
        self.body = body


def _is_timeout(exc: BaseException) -> bool:
    if isinstance(exc, TimeoutError):
        return True
    return isinstance(exc, urllib.error.URLError) and isinstance(exc.reason, TimeoutError)


def _get(url: str, timeout: float) -> bytes:
    """GET ``url`` and return the body of a 2xx response.

    Timeouts and 5xx responses are retried with exponential backoff. Any other
    HTTP status raises ``_HttpStatusError`` at once, and any other failure
    raises ``ApiError`` at once.
    """
    request = urllib.request.Request(
        url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"}
    )
    failure = ""
    status: int | None = None
    for attempt in range(_ATTEMPTS):
        if attempt:
            time.sleep(_BACKOFF_SECONDS * 2 ** (attempt - 1))
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                body: bytes = response.read()
                return body
        except urllib.error.HTTPError as exc:
            error_body = exc.read()
            if exc.code < 500:
                raise _HttpStatusError(exc.code, error_body) from None
            failure, status = f"HTTP {exc.code}", exc.code
        except (urllib.error.URLError, OSError, http.client.HTTPException) as exc:
            if not _is_timeout(exc):
                raise ApiError(f"GET {url} failed: {exc}", url=url) from exc
            failure, status = f"no response within {timeout:g} s", None
    raise ApiError(f"GET {url} failed after {_ATTEMPTS} attempts: {failure}", url=url, status=status)


def _parse_json(body: bytes, url: str) -> Any:
    try:
        return json.loads(body)
    except ValueError as exc:
        raise ApiError(f"GET {url} did not return JSON", url=url, status=200) from exc


def _rows(payload: Any, url: str) -> list[dict[str, Any]]:
    if not isinstance(payload, list) or not all(isinstance(row, dict) for row in payload):
        raise ApiError(f"GET {url} returned JSON that is not a list of objects", url=url, status=200)
    return payload


def _error_text(body: bytes) -> str:
    """The server's own explanation, if the error body carries one."""
    try:
        payload = json.loads(body)
    except ValueError:
        return ""
    if isinstance(payload, dict):
        for key in ("error", "detail"):
            if isinstance(payload.get(key), str):
                return str(payload[key])
    return ""


def metadata(qualifier: str, *, timeout: float = 30.0) -> list[dict[str, Any]]:
    """Return EBI's prediction records for a UniProt accession or a model id.

    This is a passthrough of ``GET https://alphafold.ebi.ac.uk/api/prediction/<qualifier>``.
    Each record is a dict as the API sends it (``entryId``, ``uniprotAccession``,
    ``gene``, ``taxId``, ``globalMetricValue``, ``sequence``, ``latestVersion``, ...).
    The API covers all of AlphaFold DB, which is far more than the bulk set this
    library can read, so a record here does not mean ``fetch()`` will find the entry.

    Returns ``[]`` when the identifier is well formed but AlphaFold DB has no
    prediction for it. The API answers that case with 404 and an empty JSON
    object, and an empty list is the same statement in this function's return type.

    Raises ``InvalidIdError`` (a ``ValueError``), without sending a request, when ``qualifier`` contains
    anything other than letters, digits, ``-`` and ``_`` or is longer than 64
    characters, so the value cannot change the URL path. Raises ``ApiError`` when
    the server rejects the identifier (400), answers 404 with anything other than
    an empty body, keeps failing with 5xx or timing out after three attempts, or
    returns something that is not a JSON list of objects.
    """
    if not isinstance(qualifier, str):
        raise TypeError(f"qualifier must be a str, not {type(qualifier).__name__}")
    cleaned = qualifier.strip()
    if not _QUALIFIER_RE.fullmatch(cleaned):
        raise InvalidIdError(
            f"not a UniProt accession or AlphaFold DB model id: {qualifier!r} "
            "(expected something like 'P69905' or 'AF-P69905-F1')"
        )
    url = f"{API_BASE_URL}/prediction/{urllib.parse.quote(cleaned, safe='')}"
    try:
        body = _get(url, timeout)
    except _HttpStatusError as exc:
        if exc.status == 404:
            try:
                payload = json.loads(exc.body) if exc.body.strip() else {}
            except ValueError:
                payload = None
            if payload in ({}, []):
                return []
        detail = _error_text(exc.body)
        raise ApiError(
            f"GET {url} returned HTTP {exc.status}" + (f": {detail}" if detail else ""),
            url=url,
            status=exc.status,
        ) from None
    return _rows(_parse_json(body, url), url)


def proteomes(*, timeout: float = 30.0) -> list[dict[str, Any]]:
    """Return one row per bulk tar that EBI publishes, from ``download_metadata.json``.

    Every row keeps EBI's fields (``archive_name``, ``num_predicted_structures``,
    ``size_bytes``, ``type``, and for proteome tars ``species``, ``common_name``,
    ``latin_common_name`` and ``reference_proteome``; the two Swiss-Prot rows carry
    ``label`` instead) and gains three:

    * ``url``: where the tar is, ``<FTP>/v<N>/<archive_name>`` with ``N`` taken
      from the ``_v<N>.tar`` suffix of the name (``v6`` for every row today),
    * ``tax_id``: the NCBI taxonomy id parsed from names such as
      ``UP000000625_83333_ECOLI_v6.tar``, or ``None`` for the Swiss-Prot rows,
    * ``reference_proteome``: parsed from the same name when EBI's row lacks it,
      or ``None`` for the Swiss-Prot rows.

    Raises ``ApiError`` when the file cannot be fetched (4xx at once, 5xx and
    timeouts after three attempts) or does not have the expected shape.
    """
    url = f"{FTP_BASE_URL}/download_metadata.json"
    try:
        body = _get(url, timeout)
    except _HttpStatusError as exc:
        raise ApiError(f"GET {url} returned HTTP {exc.status}", url=url, status=exc.status) from None
    out: list[dict[str, Any]] = []
    for source in _rows(_parse_json(body, url), url):
        name = source.get("archive_name")
        version = _ARCHIVE_VERSION_RE.search(name) if isinstance(name, str) else None
        if not isinstance(name, str) or not _ARCHIVE_NAME_RE.fullmatch(name) or version is None:
            raise ApiError(
                f"GET {url} has a row whose archive_name is not a versioned tar name: {name!r}",
                url=url,
                status=200,
            )
        row = dict(source)
        row["url"] = f"{FTP_BASE_URL}/v{version.group(1)}/{name}"
        parsed = _PROTEOME_NAME_RE.fullmatch(name)
        if row.get("tax_id") is None:
            row["tax_id"] = int(parsed.group(2)) if parsed else None
        if row.get("reference_proteome") is None:
            row["reference_proteome"] = parsed.group(1) if parsed else None
        out.append(row)
    return out
