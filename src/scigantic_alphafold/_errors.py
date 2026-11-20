"""Error types raised by scigantic-alphafold.

Every error the library raises on purpose derives from :class:`AlphaFoldError`, so
``except AlphaFoldError`` catches all of them. Each one also derives from the closest
built-in type (``LookupError``, ``RuntimeError``, ``ValueError``) so that generic
handlers keep working.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence


class AlphaFoldError(Exception):
    """Base class of every error this library raises on purpose."""


class InvalidIdError(AlphaFoldError, ValueError):
    """The text is not a UniProt accession or an AlphaFold DB entry id."""


class NotInBulkSetError(AlphaFoldError, LookupError):
    """The accession or entry is not in any AlphaFold DB v6 bulk tar."""


class AmbiguousEntryError(AlphaFoldError, LookupError):
    """A bare accession has several entries and none is the single canonical one.

    ``candidates`` lists the entry ids that exist; ask for one of them by its full id.
    """

    def __init__(self, message: str, candidates: Sequence[str] = ()) -> None:
        super().__init__(message)
        self.candidates: tuple[str, ...] = tuple(candidates)


class StaleIndexError(AlphaFoldError, RuntimeError):
    """The tar on the server no longer matches the index (or the bytes do not gunzip)."""


class IndexUnavailableError(AlphaFoldError, RuntimeError):
    """The index could not be found, downloaded or opened."""


class IndexSchemaError(AlphaFoldError, RuntimeError):
    """The index has a schema version this release does not understand."""


class FetchError(AlphaFoldError, RuntimeError):
    """The range request failed: an HTTP error, a server that ignored Range, or retries ran out.

    ``status`` is the HTTP status when there was one.
    """

    def __init__(self, message: str, status: int | None = None) -> None:
        super().__init__(message)
        self.status = status


class FetchManyError(AlphaFoldError):
    """Some ids of a ``fetch_many`` call failed. Nothing that succeeded is lost.

    ``results`` maps each id that succeeded to its text, exactly as the return value
    would have. ``errors`` maps each id that failed to the exception it raised.
    """

    def __init__(self, results: Mapping[str, str], errors: Mapping[str, Exception]) -> None:
        self.results: dict[str, str] = dict(results)
        self.errors: dict[str, Exception] = dict(errors)
        shown = "; ".join(f"{key}: {exc}" for key, exc in list(self.errors.items())[:3])
        more = len(self.errors) - 3
        if more > 0:
            shown += f"; and {more} more"
        super().__init__(
            f"{len(self.errors)} of {len(self.errors) + len(self.results)} ids failed "
            f"({len(self.results)} succeeded and are in .results; failures are in .errors): {shown}"
        )
