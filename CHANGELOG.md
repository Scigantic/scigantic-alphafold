# Changelog

All notable changes to this project are recorded here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 0.1.0 - 2026-11-20

First release.

### Added

- `fetch`, `fetch_bytes`, `download` and `fetch_many`: read single AlphaFold DB v6 structures (mmCIF or PDB) with one HTTP range request each against EBI's public bulk tars at `https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/`.
- `locate`, `has`, `entries`, `tars` and `index_info`: answer from the offset index without a request to EBI.
- The offset index `index-v6.sqlite` (schema version 1), downloaded on first use from this repository's `index-v6` release and cached in `~/.cache/scigantic-alphafold`. `SCIGANTIC_ALPHAFOLD_INDEX` and `SCIGANTIC_ALPHAFOLD_CACHE` override the source and the cache directory. `refresh_index()` downloads it again.
- `NotInBulkSetError` for accessions outside the bulk set (Swiss-Prot plus 46 proteome tars) and `StaleIndexError` when a tar on the server no longer matches the index.
- `metadata()`: passthrough of EBI's prediction API, `https://alphafold.ebi.ac.uk/api/prediction/<qualifier>`. Returns `[]` for an accession with no prediction, raises `InvalidIdError` for input that is not an accession or model id, and raises `ApiError` for HTTP errors.
- `proteomes()`: passthrough of EBI's `download_metadata.json`, with `url`, `tax_id` and `reference_proteome` added to each row.
- The `scigantic-alphafold` command: `get`, `locate`, `info`, `refresh`, and the index builders `build-part` and `merge`.
- `examples/quickstart.ipynb`.

### Notes

- Coverage is the v6 bulk set only, coordinates only. There are no PAE files and no MSAs.
- No runtime dependencies. Python 3.10 or newer.
- Every request sends the User-Agent `scigantic-alphafold/<version> (+https://github.com/Scigantic/scigantic-alphafold)`.
