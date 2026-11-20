"""Read single AlphaFold DB v6 structures with one HTTP range request.

    import scigantic_alphafold as af

    text = af.fetch("P69905")                  # mmCIF text
    path = af.download("P69905", "structures", fmt="pdb")

The bytes come straight from EBI's public bulk tars; a prebuilt offset index says where
each structure sits. Importing this package touches neither the network nor the disk:
the index is located and opened on the first call that needs it.
"""

from __future__ import annotations

from ._api import ApiError, metadata, proteomes
from ._errors import (
    AlphaFoldError,
    AmbiguousEntryError,
    FetchError,
    FetchManyError,
    IndexSchemaError,
    IndexUnavailableError,
    InvalidIdError,
    NotInBulkSetError,
    StaleIndexError,
)
from ._fetch import download, fetch, fetch_bytes, fetch_many
from ._index import Member, Tar, entries, has, index_info, locate, refresh_index, tars
from ._version import __version__

__all__ = [
    "AlphaFoldError",
    "AmbiguousEntryError",
    "ApiError",
    "FetchError",
    "FetchManyError",
    "IndexSchemaError",
    "IndexUnavailableError",
    "InvalidIdError",
    "Member",
    "NotInBulkSetError",
    "StaleIndexError",
    "Tar",
    "__version__",
    "download",
    "entries",
    "fetch",
    "fetch_bytes",
    "fetch_many",
    "has",
    "index_info",
    "locate",
    "metadata",
    "proteomes",
    "refresh_index",
    "tars",
]
