<h1 align="center">scigantic-alphafold</h1>

<p align="center">
    <a href="https://github.com/Scigantic/scigantic-alphafold/actions/workflows/ci.yml">
        <img alt="CI" src="https://github.com/Scigantic/scigantic-alphafold/actions/workflows/ci.yml/badge.svg" /></a>
    <a href="https://pypi.org/project/scigantic-alphafold/">
        <img alt="PyPI" src="https://img.shields.io/pypi/v/scigantic-alphafold" /></a>
    <a href="https://pypi.org/project/scigantic-alphafold/">
        <img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/scigantic-alphafold" /></a>
    <a href="https://github.com/Scigantic/scigantic-alphafold/blob/main/LICENSE">
        <img alt="License" src="https://img.shields.io/github/license/Scigantic/scigantic-alphafold" /></a>
</p>

Read a single [AlphaFold DB](https://alphafold.ebi.ac.uk/) v6 structure with one HTTP range request against EBI's public bulk tars, using a prebuilt byte-offset index. No mirror, no account, no key, and no multi-gigabyte tar download to get one protein.

```python
import scigantic_alphafold as af

cif = af.fetch("P69905")                                 # hemoglobin alpha, mmCIF text
pdb = af.fetch("P69905", fmt="pdb")                      # the same model as PDB text
path = af.download("P69905", "structures")               # writes the .cif, returns a pathlib.Path
many = af.fetch_many(["P0A7G6", "P0A6F5", "P00350"])     # dict of id to text, 4 connections by default

af.has("P69905")                                         # True: the accession is in the bulk set
af.locate("P69905")                                      # [Member(entry_id, fmt, tar, url, offset, length), ...]
ids = af.entries(proteome="UP000000625")                 # iterator over the entry ids of the E. coli tar

af.metadata("P69905")                                    # EBI prediction API records, no index needed
af.proteomes()                                           # EBI's list of bulk tars, with url and tax_id
```

A longer walkthrough, including per-residue pLDDT from the B-factor column, is in [examples/quickstart.ipynb](https://github.com/Scigantic/scigantic-alphafold/blob/main/examples/quickstart.ipynb).

**Coverage is limited to the v6 bulk set: Swiss-Prot plus 46 proteomes, coordinates only.** Most of AlphaFold DB is outside that set and cannot be read with this library. See [What is covered](#what-is-covered) before you depend on it.

## Installation

```console
$ pip install scigantic-alphafold
```

Python 3.10 or newer. There are no runtime dependencies; the package uses only the standard library.

## Why this exists

EBI publishes AlphaFold DB v6 in bulk as 48 tar files at `https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/`, 153.2 GB in total. The smallest is 173 MB (*Helicobacter pylori*), the human proteome is 5.2 GB and the Swiss-Prot mmCIF tar is 40.1 GB. A tar has no table of contents, so finding one protein in it normally means downloading the whole file.

Inside each tar every structure is its own small gzip file, and EBI's server honors HTTP `Range` requests. If you know the byte offset and length of a member, you can read that member alone. This package ships that knowledge as an index and does the range read.

On 2026-11-20 we observed that EBI's per-file URLs (`https://alphafold.ebi.ac.uk/files/AF-<accession>-F1-model_v6.cif`) returned HTTP 403 to scripted requests, while the prediction API and the FTP area answered normally. We know of no announcement about this and do not know whether it is intended or permanent. This library does not use the per-file endpoint. It reads only the bulk files that EBI publishes for download.

## What is covered

The index covers the 48 tars listed in EBI's [`download_metadata.json`](https://ftp.ebi.ac.uk/pub/databases/alphafold/download_metadata.json) (counts and sizes checked on 2026-11-20):

| set | tars | structures | size |
|---|---|---|---|
| Swiss-Prot, mmCIF | 1 | 550,122 | 40.1 GB |
| Swiss-Prot, PDB | 1 | 550,122 | 28.6 GB |
| Model organism proteomes | 16 | 328,602 | 49.9 GB |
| Global health proteomes | 30 | 227,244 | 34.7 GB |

The proteome tars hold both formats for every entry. The 16 model organisms include human (UP000005640), mouse, rat, zebrafish, fruit fly, *C. elegans*, yeast, *E. coli* and *Arabidopsis*. `af.proteomes()` returns the full list with species, taxonomy id, reference proteome, structure count, size and URL. A protein can be in both Swiss-Prot and a proteome tar, so the counts above overlap.

What is not covered:

- **Any accession outside those tars.** AlphaFold DB holds predictions for far more UniProt entries than the bulk set. For an accession that is in no v6 bulk tar, `fetch`, `fetch_bytes` and `download` raise `NotInBulkSetError` (a `LookupError`). `af.locate(accession)` returns an empty list and `af.has(accession)` returns `False`, both without raising.
- **PAE.** The bulk tars hold coordinates only. There are no predicted aligned error files.
- **MSAs.** There are no alignment files.
- **Other files EBI links from an entry page**, such as the per-residue confidence JSON, binary CIF and AlphaMissense annotations. Per-residue pLDDT is still available, because it is stored in the B-factor column of the coordinates.
- **Older versions.** The index describes v6 only.

`af.metadata()` is the exception to the coverage limit. It calls EBI's prediction API, which answers for all of AlphaFold DB, so a record from `metadata()` does not mean the structure is in the bulk set.

According to EBI's [download page](https://alphafold.ebi.ac.uk/download), proteins longer than 2,700 residues are predicted as overlapping 1,400-residue fragments (`F1`, `F2`, ...), and those fragments are in the human proteome tar. An entry id such as `AF-Q8WZ42-F2` names one fragment.

## Public API

| function | returns |
|---|---|
| `fetch(id, fmt=None)` | the decompressed text of one structure (`str`) |
| `fetch_bytes(id, fmt=None)` | the gzip bytes exactly as stored in the tar |
| `download(id, out_dir, fmt=None, decompress=True)` | `pathlib.Path` of the file written |
| `fetch_many(ids, fmt=None, workers=4, on_error="raise")` | `dict[str, str]` keyed by the ids as given, fetched on a thread pool |
| `locate(id, fmt=None)` | `list[Member]`; a `Member` has `entry_id`, `fmt`, `tar`, `url`, `offset`, `length` |
| `has(id)` | `bool`, answered from the index without a request to EBI |
| `entries(proteome=..., fmt="cif")`, `entries(tar=..., fmt="cif")` | iterator over the entry ids of one tar |
| `tars()` | `list[Tar]`, the tars the index covers |
| `index_info()` | `dict` describing the index in use |
| `refresh_index()` | `pathlib.Path` of the freshly downloaded index |
| `metadata(qualifier, timeout=30.0)` | `list[dict]` from EBI's prediction API |
| `proteomes(timeout=30.0)` | `list[dict]` from EBI's `download_metadata.json` |

`id` is a UniProt accession (`P69905`) or an AlphaFold DB entry id (`AF-P69905-F1`). `fmt` is `"cif"` or `"pdb"`; left as `None` it means mmCIF. A member file name such as `AF-P69905-F1-model_v6.pdb.gz` is also accepted as an id and selects its own format; a file name for another version, or one whose extension conflicts with an explicit `fmt`, raises `InvalidIdError`.

How an id resolves:

- A bare accession resolves to `AF-<accession>-F1` when that is its only fragment.
- An accession that is split into several fragments (`F1`, `F2`, ...) raises `AmbiguousEntryError` listing the entry ids, because `F1` alone is only the first part of the protein. Pass the entry id of the fragment you want. `locate()` returns all of them.
- An isoform is its own accession (`P12345-2`).
- When the same entry is in two tars (a Swiss-Prot protein that is also in a proteome tar), the tar whose name sorts first is read, which in practice is the proteome tar. `locate()` lists that copy first.

`fetch_many` attempts every id. If any fail it raises `FetchManyError` after the rest have finished, with the successes in `.results` and the failures in `.errors`. Pass `on_error="skip"` to get only the successes back.

Errors:

- `NotInBulkSetError` (a `LookupError`): the accession is in no v6 bulk tar.
- `AmbiguousEntryError`: a bare accession has several fragments, or only isoform entries.
- `StaleIndexError` (a `RuntimeError`): EBI replaced a tar after the index was built. Call `refresh_index()`.
- `FetchError`: EBI answered a range request with an HTTP error. It carries `status`.
- `FetchManyError`: one or more ids in `fetch_many` failed. It carries `results` and `errors`.
- `IndexUnavailableError`, `IndexSchemaError`: the index could not be downloaded, or was built for a newer version of this library.
- `InvalidIdError` (a `ValueError`): the id is not an accession or entry id.
- `ApiError` (a `RuntimeError`): `metadata()` or `proteomes()` got an HTTP error, a timeout or an unexpected body. It carries `status` and `url`.

### `metadata()` and `proteomes()`

Both are thin passthroughs. They need no index, and they return what EBI sends.

`metadata("P69905")` calls `https://alphafold.ebi.ac.uk/api/prediction/P69905` and returns a list of records with keys such as `entryId`, `uniprotAccession`, `gene`, `taxId`, `organismScientificName`, `sequence`, `globalMetricValue` (mean pLDDT) and `latestVersion`. The API answers a well-formed accession that has no prediction with HTTP 404 and an empty JSON object; `metadata()` returns `[]` for that case. A value that is not made of letters, digits, `-` and `_` raises `InvalidIdError` (a `ValueError`) before any request is sent. An identifier the server rejects (HTTP 400) raises `ApiError` with the server's message.

`proteomes()` returns one dict per bulk tar with EBI's fields (`archive_name`, `type`, `num_predicted_structures`, `size_bytes`, and for proteome tars `species`, `common_name` and `reference_proteome`) plus `url`, the location of the tar, and `tax_id`, parsed from names such as `UP000000625_83333_ECOLI_v6.tar`. `tax_id` and `reference_proteome` are `None` for the two Swiss-Prot rows.

Both functions retry timeouts and 5xx responses twice with backoff, and never retry a 4xx.

## How it works

1. **The index.** One SQLite file, `index-v6.sqlite`, has a row per structure file: entry id, format, tar, byte offset and length of the gzip member. It also records each tar's size, `Last-Modified` and `ETag` at build time.
2. **Where the index is hosted.** It is a release asset of this repository: `https://github.com/Scigantic/scigantic-alphafold/releases/download/index-v6/index-v6.sqlite.gz`. EBI does not host or maintain it.
3. **First use.** The first call that needs the index downloads it to the cache directory, `~/.cache/scigantic-alphafold` by default. The download goes to a temporary name and is renamed into place, so an interrupted download does not leave a partial index. The download is checked against the SHA-256 published next to it in `index-v6.json`; a mismatch is refused, which can happen for a moment while the index is being republished. Later calls open the cached file. `import scigantic_alphafold` alone does not download anything.
4. **One range read.** `fetch` looks up the entry, sends one `GET` with a `Range` header to the tar on `ftp.ebi.ac.uk`, and gunzips the bytes that come back. A member is small: the mmCIF of a 217-residue *H. pylori* protein (`AF-G5EB01-F1`) is 47,367 bytes gzipped, read out of a 173 MB tar.
5. **Staleness check.** If the `Last-Modified` or `ETag` of the response differs from what the index recorded, the bytes do not gunzip, or the file name recorded in the gzip header is not the entry that was asked for, the offsets can no longer be trusted and the call raises `StaleIndexError`. `af.refresh_index()` downloads the current index again.

These environment variables change where the index comes from and what it may point at:

| variable | meaning |
|---|---|
| `SCIGANTIC_ALPHAFOLD_INDEX` | a local path or a URL of an index to use instead of the default, as `.sqlite` or `.sqlite.gz` |
| `SCIGANTIC_ALPHAFOLD_CACHE` | the cache directory, instead of `~/.cache/scigantic-alphafold` |
| `SCIGANTIC_ALPHAFOLD_ALLOW_ANY_HOST` | set to `1` to accept an index whose tars are not on `https://ftp.ebi.ac.uk`, for a private mirror or for tests |

By default the library only sends range requests to `https://ftp.ebi.ac.uk`, whatever the index says, so a wrong or tampered index cannot point it at another host.

`af.index_info()` reports which index is in use. The index is built by reading every tar once and recording where each member starts; [builders/README.md](https://github.com/Scigantic/scigantic-alphafold/blob/main/builders/README.md) describes the build and how to reproduce it.

## Being a polite client of EBI

The FTP area is a shared public service that EMBL-EBI runs at its own cost. This library is meant for reading some structures, not for copying the database.

- Every request carries an honest `User-Agent`: `scigantic-alphafold/<version> (+https://github.com/Scigantic/scigantic-alphafold)`. Do not change it to look like a browser.
- `fetch_many` uses 4 connections by default. Keep `workers` near that, and do not run several processes in parallel against the FTP.
- Requests are retried with backoff only after timeouts and 5xx responses.
- If you need a whole proteome or all of Swiss-Prot, do not loop over `fetch`. Download EBI's tar once: `af.proteomes()` lists the URLs. For large transfers EBI recommends [Globus](https://app.globus.org/file-manager?origin_id=47772002-3e5b-4fd3-b97c-18cee38d6df2&origin_path=%2Fpub%2Fdatabases%2Falphafold%2F) (endpoint "EMBL-EBI Public Data", path `/pub/databases/alphafold/`), as described in the [FTP README](https://ftp.ebi.ac.uk/pub/databases/alphafold/README.txt).
- If EBI limits or closes range access to the tars, this library will stop working, and it will not try to get around that.

## Command line

`scigantic-alphafold` wraps the same functions.

```console
$ scigantic-alphafold get P69905                      # writes the mmCIF file to the current directory
$ scigantic-alphafold get P69905 P68871 -f pdb -o structures
$ scigantic-alphafold get P69905 --stdout | head
$ scigantic-alphafold locate P69905                   # tar, offset and length
$ scigantic-alphafold info                            # the index in use
$ scigantic-alphafold refresh                         # download the index again
```

Two more subcommands, `build-part <tar>` and `merge <parts dir> -o index-v6.sqlite`, build the index; `merge --only <tar>...` builds a partial index from the tars you streamed. They are documented in [builders/README.md](https://github.com/Scigantic/scigantic-alphafold/blob/main/builders/README.md).

## Data license and citation

This project is not affiliated with, endorsed by or supported by EMBL-EBI or Google DeepMind. It redistributes no AlphaFold DB data: the structures are read from EBI's servers at the time of the call, and the index holds only file names, byte offsets and lengths.

The [AlphaFold DB site](https://alphafold.ebi.ac.uk/) states: "Data is available for academic and commercial use, under a [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/) licence." It also states: "AlphaFold Data Copyright (2022) DeepMind Technologies Limited." (both read on 2026-11-20). CC BY 4.0 requires attribution. The site's disclaimer says the data "is for theoretical modelling only" and has "not been validated for, and are not approved for, any clinical use"; read it before relying on a structure. Use of AlphaFold DB is also subject to the [EMBL-EBI Terms of Use](https://www.ebi.ac.uk/about/terms-of-use).

If you use structures obtained through this library, cite AlphaFold DB and AlphaFold, not this package. The papers the site asked users to cite on 2026-11-20:

- Bertoni D. et al. *AlphaFold Protein Structure Database 2025: a redesigned interface and updated structural coverage*. Nucleic Acids Research (2025). [academic.oup.com/nar/article/54/D1/D358/8340156](https://academic.oup.com/nar/article/54/D1/D358/8340156). Listed on the home page for use of the resource.
- Fleming J. et al. *AlphaFold Protein Structure Database and 3D-Beacons: New Data and Capabilities*. Journal of Molecular Biology (2025). [sciencedirect.com/science/article/pii/S0022283625000336](https://www.sciencedirect.com/science/article/pii/S0022283625000336). Listed on the download page for use of a prediction.
- Jumper J. et al. *Highly accurate protein structure prediction with AlphaFold*. Nature (2021). [doi:10.1038/s41586-021-03819-2](https://doi.org/10.1038/s41586-021-03819-2). Listed on the download page for predictions from the UniProt dataset.

The site's guidance can change. Check [alphafold.ebi.ac.uk](https://alphafold.ebi.ac.uk/) for the current text before you publish.

This package's code is MIT-0. That does not extend to the data.

## Testing

```console
$ pip install -e ".[dev]"
$ pytest
$ mypy --strict src/scigantic_alphafold
```

The default run is offline: tests use local files and an HTTP server on localhost. Tests marked `live` send a few small requests to EBI and are skipped unless `SCIGANTIC_ALPHAFOLD_LIVE=1` is set.

## License

MIT-0 for the code. See [Data license and citation](#data-license-and-citation) for the data.
