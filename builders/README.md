# Building the index

`scigantic-alphafold` reads one structure with one HTTP range request. That
only works because it knows, for every model in every AlphaFold DB v6 bulk
tar at EBI, where the model's bytes start and how long they are. This page
describes how that index is produced and kept current. None of it is needed
to use the library.

The code is `src/scigantic_alphafold/_build.py` (standard library only) and
the workflow is `.github/workflows/build-index.yml`.

## What gets built

EBI publishes 48 v6 tars (153 GB) under
`https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/` and lists them, with the
number of predicted structures in each, in
`https://ftp.ebi.ac.uk/pub/databases/alphafold/download_metadata.json`.

Every member is an individually gzipped model named
`AF-<accession>-F<n>-model_v6.cif.gz` or `...pdb.gz`. A tar stores each member
as a 512-byte header followed by the member's bytes, so the gzip stream of one
model is a plain byte range of the tar. The index records that range.

The build has three steps.

1. **One part file per tar.** `build_part` streams the tar over HTTP, walks it
   with `tarfile` in stream mode and writes `<tar name>.tsv.gz`. The tar is
   never written to disk.
2. **Merge.** `merge_parts` checks the parts against `download_metadata.json`
   and writes `index-v6.sqlite` and `index-v6.sqlite.gz`.
3. **Spot check.** `verify_index` reads random members back from EBI with the
   offsets in the new index and checks that each one is the right model.

## Part file format

`<tar name>.tsv.gz` is gzipped text:

```
#{"etag": "\"a570c00-63d1a3a8c9640\"", "last_modified": "Sun, 24 Aug 2025 10:52:01 GMT", "name": "UP000000429_85962_HELPY_v6.tar", "size_bytes": 173476864, "url": "https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/UP000000429_85962_HELPY_v6.tar"}
AF-G5EB01-F1	cif	512	47367
AF-G5EB01-F1	pdb	48640	33194
...
#end	{"bytes_read": 173476864, "end_of_members": 173475840, "n_cif": 1540, "n_pdb": 1540, "n_skipped": 0, "reconnects": 0}
```

* Line 1 is `#` plus a JSON object taken from the HTTP response: the tar's
  name, URL, Content-Length, Last-Modified and ETag.
* Each row is `entry_id`, `fmt`, `offset`, `length`, tab separated, in tar
  order. `offset` is the position of the member's data (not its header) and
  `length` is the size of the `.gz`. Requesting
  `Range: bytes=<offset>-<offset+length-1>` returns exactly that file.
* A member whose name does not match the model naming is written as
  `#skipped<TAB><member name>` and is not indexed. Nothing is dropped
  silently.
* The last line is `#end` plus the counts. A part without it was cut short
  and the merge rejects it. The part is also written under a temporary name
  and renamed when complete, so a killed run leaves no file that looks
  finished.

## Dropped connections

The stream that `tarfile` reads from reconnects by itself. When the connection
drops, times out or the server answers 5xx, it sends
`Range: bytes=<position>-` (with `If-Range`) and continues from the exact byte
it had reached, so nothing is downloaded twice and the parser is not
restarted. Each reconnect checks that the ETag, Last-Modified and total size
are unchanged; if the tar was replaced mid-build the offsets already recorded
are worthless and the build stops. It gives up after 8 consecutive failed
attempts (backoff 2, 4, 8 ... 60 seconds) or 500 reconnects in total. 4xx
responses are never retried.

`tarfile` stops quietly at the first block it cannot parse, so after the last
member the builder reads to the end of the tar and fails if anything other
than zero padding follows.

## Checks in the merge

`merge_parts` prints a table with one line per tar and writes nothing unless
all of these hold:

* every tar in `download_metadata.json` has a part, and every part is listed;
* every part is readable, ends with `#end`, and its rows agree with the counts
  on that line;
* no part has zero rows or the same entry and format twice;
* the member count per format equals the tar's `num_predicted_structures`.

The count rule was worked out from real tars on 2026-11-20:

| tar | type | num_predicted_structures | cif members | pdb members | skipped |
| --- | --- | --- | --- | --- | --- |
| UP000000429_85962_HELPY_v6.tar | global_health | 1540 | 1540 | 1540 | 0 |
| UP000000799_192222_CAMJE_v6.tar | global_health | 1620 | 1620 | 1620 | 0 |
| UP000000579_71421_HAEIN_v6.tar | global_health | 1660 | 1660 | 1660 | 0 |
| UP000000806_272631_MYCLE_v6.tar | global_health | 1602 | 1602 | 1602 | 0 |
| UP000000805_243232_METJA_v6.tar | proteome | 1773 | 1773 | 1773 | 0 |

So a proteome or global health tar must hold exactly
`num_predicted_structures` cif members and the same number of pdb members.
The two Swiss-Prot tars hold one format each (`swissprot_cif_v6.tar` only cif,
`swissprot_pdb_v6.tar` only pdb), and the expected count is
`num_predicted_structures` (550122) for that format and zero for the other.
Their member naming was checked by reading windows at the start, at 25, 50 and
75 percent and at the end of each tar, and it is the same as in the proteome
tars. Their counts have not been measured, because that takes a full pass over
40 GB and 29 GB; the first full build measures them. The 43 other tars are
also unmeasured until then.

If a tar turns out to differ, the merge fails and reports the exact numbers,
for example `12 pdb members, expected 13 (num_predicted_structures=13,
difference -1)`. The way to accept a real difference is an explicit entry in
`KNOWN_COUNT_EXCEPTIONS` in `_build.py`, with a comment saying what was
measured, not a looser check.

A difference between the Content-Length of a tar and `size_bytes` in
`download_metadata.json` is shown in the table but does not fail the merge.
The index stores what the server sent.

In the index, `tars.n_members` is the number of indexed rows for that tar
(cif plus pdb), so it is twice `num_predicted_structures` for a proteome tar.

## Refresh

The workflow runs on demand and every Monday at 04:17 UTC.

1. `plan` reads `download_metadata.json` and the `index-v6.json` currently
   attached to the `index-v6` release, then sends one HEAD request per tar.
   If every tar's size, Last-Modified and ETag match and the set of tars is
   the same, the scheduled run stops there. A manual run rebuilds unless
   `force` is unticked.
2. `build` is a matrix with one job per tar, eight at a time, each uploading
   its part file as an artifact. Jobs that fail can be re-run on their own.
3. `publish` downloads the parts, merges, spot checks 20 random members plus 2
   from every tar, writes `index-v6.json` (build time, per-tar name, size,
   Last-Modified, ETag and member count, total rows, SHA-256 of the `.gz`) and
   uploads both files to the `index-v6` release. It does not run unless every
   part was built.

Things to know:

* A hosted runner job is limited to 6 hours and the build jobs are set to 350
  minutes. `swissprot_cif_v6.tar` is 40,062,921,728 bytes, so it needs a
  sustained **1.9 MB/s**; `swissprot_pdb_v6.tar` (28,580,050,944 bytes) needs
  1.4 MB/s. Parsing is not the limit: the stream walk runs at about 340 MB/s
  against a local server. A home connection measured 1.1 to 2.0 MB/s from EBI,
  which would not be enough. The rate from a GitHub runner has not been
  measured yet. A job that hits the time limit starts that tar again from the
  beginning on re-run: resume works within a run, not across runs.
* Runners have about 14 GB of free disk. The build writes only part files
  (about 35 KB for a 180 MB tar, an estimated 6 MB for a Swiss-Prot tar).
* GitHub disables scheduled workflows in a repository with no activity for 60
  days. Re-enable it from the Actions tab if that happens.
* A reader holding an older index notices a replaced tar by itself, because
  the ETag of the range response no longer matches, and asks for a refresh.

## Size

The local index of the five tars above has 16,390 rows in 557,056 bytes
(181,035 gzipped), 34 bytes per row. All 48 tars come to 2,211,936 rows, which
projects to about 77 MB, and roughly 25 to 31 MB gzipped. The gzipped figure
is an estimate from a synthetic index of that size; the real one is known
after the first full build.

## Run one tar locally

From the repository root, with nothing installed:

```
export PYTHONPATH=src
T=UP000000429_85962_HELPY_v6.tar

python -m scigantic_alphafold._build build-part $T -o build/local/parts
python -m scigantic_alphafold._build merge build/local/parts -o build/local/index-v6.sqlite --only $T
python -m scigantic_alphafold._build verify build/local/index-v6.sqlite -n 20
python -m scigantic_alphafold._build summary build/local/index-v6.sqlite -o build/local/index-v6.json
```

`--only` restricts the tar list so the merge does not fail on the 47 tars you
did not build; list every tar you have parts for. The same functions are
available from Python:

```python
from scigantic_alphafold._build import build_part, merge_parts, verify_index

build_part("UP000000429_85962_HELPY_v6.tar", "build/local/parts")
```

Then point the library at the result:

```
export SCIGANTIC_ALPHAFOLD_INDEX=$PWD/build/local/index-v6.sqlite
```

Use the small tars for this (HELPY, CAMJE, HAEIN, MYCLE and METJA are about
180 MB each) and keep to a few connections at a time.

## Tests

`pytest tests/test_build.py` is offline: it serves a synthetic tar from a
local HTTP server that supports Range and can drop connections, answer 503 or
change its ETag. Three tests marked `live` talk to EBI and only run with
`SCIGANTIC_ALPHAFOLD_LIVE=1`; one of them streams the 173 MB HELPY tar.
