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

1. **Part files.** `build_part` streams a tar over HTTP, walks it with
   `tarfile` in stream mode and writes a part file. The tar is never written
   to disk. A small tar is one part; a large tar is built as several
   byte-range segments, one part each, by separate jobs.
2. **Merge.** `merge_parts` checks the parts against `download_metadata.json`,
   proves that the segments of a split tar cover it exactly, and writes
   `index-v6.sqlite` and `index-v6.sqlite.gz`.
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
  name, URL, total size, Last-Modified and ETag.
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

## Segments

A hosted runner gets about 1.2 MB/s per stream from EBI (ECOLI, 456 MB, took
6 min 44 s on 2026-11-20). At that rate `swissprot_cif_v6.tar` (40 GB) needs
more than 9 hours, and a job may run for 6. So a tar larger than
`SEGMENT_THRESHOLD_BYTES` (6 GiB) is split into equal byte ranges of at most
`SEGMENT_BYTES` (4 GiB, about an hour each), and each range is indexed by its
own job:

```
python -m scigantic_alphafold._build build-part swissprot_cif_v6.tar --segment 3/10
```

writes `swissprot_cif_v6.tar.seg03of10.tsv.gz`. Its first line has one more
key and its last line three more:

```
#{..., "segment": {"index": 3, "count": 10, "start": 8012584448, "end": 12018876928}}
...
#end	{..., "first_header": <offset>, "stop": <offset>, "reached_end": false}
```

How a segment is read:

* It starts at a 512-aligned offset that is almost always in the middle of a
  member. Tar headers sit on 512-byte boundaries, so it tries every block
  from there as a header and takes the first one that has the ustar magic and
  a valid header checksum, is a regular file whose name matches the model
  naming, whose data starts with the gzip magic, and which is followed, right
  after the member, by another valid header or the zero block that ends the
  archive. That offset is `first_header`.
* From there it parses member by member exactly as a whole-tar build does,
  with the same reconnect logic, and stops at the first model member header
  at or past the end of its range. That offset is `stop`; it is the next
  segment's first member. A member belongs to the segment in which its header
  starts.
* The last segment runs into the end-of-archive marker instead. It then
  checks that only zero padding follows, sets `stop` to where the marker
  starts and `reached_end` to true.
* A segment reads its own range plus at most the one member that straddles
  its end and a 64 KB read buffer. If no header turns up within 64 MB past
  its end it fails.
* A segment that lies entirely inside one member has no rows and
  `first_header == stop`. A segment that starts inside the last member or in
  the padding finds no header before the end of the file and records
  `first_header` and `stop` as null.

The header test is strict, but the merge does not rely on it. See the next
section.

Segments need members with plain headers. A tar that puts a pax extended
header in front of every member would not resynchronize onto the same offsets
and the merge would reject it; the EBI tars do not do that.

## Checks in the merge

`merge_parts` prints a table with one line per tar and writes nothing unless
all of these hold:

* every tar in `download_metadata.json` has a part (or a full set of
  segments), and every part belongs to a listed tar;
* every part is readable, ends with `#end`, and its rows agree with the counts
  on that line;
* no tar has zero rows or the same entry and format twice, and its rows are in
  increasing tar order without overlap;
* the member count per format equals the tar's `num_predicted_structures`;
* for a tar built in segments, the segments tile it exactly.

**Tiling.** All segments of a tar must be there, from the same layout
(`...of10`) and the same version of the file (same size, ETag and
Last-Modified). Segment 1 must have started at header offset 0. Every
segment's `stop` must equal the next segment's `first_header`. The last
segment must have reached the end of the archive. A whole-tar part and
segment parts for the same tar may not be mixed.

This is what makes a wrong resynchronization impossible to publish. A
segment's `stop` was reached by walking the tar member by member from its
first header, so if the first header is real, so is the stop. Segment 1
starts at offset 0, which is real. Requiring each next segment to have
started at exactly the previous stop therefore proves, one segment after the
other, that every segment walked the one true chain of members. A segment
that had locked onto a forged header inside some member's data would have
started before the previous segment's stop, and the merge reports it with
both offsets, in this form:

```
segment 1/2 stopped at the member header at byte <A> but segment 2/2
starts at the header at byte <B>: an overlap of <A - B> bytes
```

The tests build exactly that forgery (a block inside a member with valid
magic, checksum and name, gzip magic after it and a size that leads to the
real next header) and check that the merge refuses it, along with gaps,
overlaps, missing segments and a last segment that stops early.

**Counts.** The count rule was worked out from real tars on 2026-11-20:

| tar | type | num_predicted_structures | cif members | pdb members | skipped |
| --- | --- | --- | --- | --- | --- |
| UP000000429_85962_HELPY_v6.tar | global_health | 1540 | 1540 | 1540 | 0 |
| UP000000799_192222_CAMJE_v6.tar | global_health | 1620 | 1620 | 1620 | 0 |
| UP000000579_71421_HAEIN_v6.tar | global_health | 1660 | 1660 | 1660 | 0 |
| UP000000806_272631_MYCLE_v6.tar | global_health | 1602 | 1602 | 1602 | 0 |
| UP000000805_243232_METJA_v6.tar | proteome | 1773 | 1773 | 1773 | 0 |
| UP000000625_83333_ECOLI_v6.tar | proteome | 4370 | 4370 | 4370 | 0 |
| UP000000559_237561_CANAL_v6.tar | proteome | 5973 | 5973 | 5973 | 0 |

The first five were built on a laptop, the last two by the hosted runner. So
a proteome or global health tar must hold exactly `num_predicted_structures`
cif members and the same number of pdb members. The two Swiss-Prot tars hold
one format each (`swissprot_cif_v6.tar` only cif, `swissprot_pdb_v6.tar` only
pdb), and the expected count is `num_predicted_structures` (550122) for that
format and zero for the other. Their member naming was checked by reading
windows at the start, at 25, 50 and 75 percent and at the end of each tar,
and it is the same as in the proteome tars. Their counts are not measured
until the first full build.

If a tar turns out to differ, the merge fails and reports the exact numbers,
for example `12 pdb members, expected 13 (num_predicted_structures=13,
difference -1)`. The way to accept a real difference is an explicit entry in
`KNOWN_COUNT_EXCEPTIONS` in `_build.py`, with a comment saying what was
measured, not a looser check.

A difference between the size of a tar on the server and `size_bytes` in
`download_metadata.json` is shown in the table but does not fail the merge.
The index stores what the server sent.

In the index, `tars.n_members` is the number of indexed rows for that tar
(cif plus pdb), so it is twice `num_predicted_structures` for a proteome tar.

## Dropped connections

The stream that `tarfile` reads from reconnects by itself. When the connection
drops, times out or the server answers 5xx, it sends
`Range: bytes=<position>-` (with `If-Range`) and continues from the exact byte
it had reached, so nothing is downloaded twice and the parser is not
restarted. Each reconnect checks that the ETag, Last-Modified and total size
are unchanged; if the tar was replaced mid-build the offsets already recorded
are worthless and the build stops. It gives up after 8 consecutive failed
attempts (backoff 2, 4, 8 ... 60 seconds) or 500 reconnects in total. 4xx
responses are never retried. This works the same inside a segment.

`tarfile` stops quietly at the first block it cannot parse, so when parsing
ends the builder reads to the end of the tar and fails if anything other than
zero padding follows.

## Refresh

The workflow runs on demand and every Monday at 04:17 UTC.

1. `plan` reads `download_metadata.json` and the `index-v6.json` currently
   attached to the `index-v6` release, then sends one HEAD request per tar.
   If every tar's size, Last-Modified and ETag match and the set of tars is
   the same, the scheduled run stops there. Otherwise it emits the list of
   build jobs: one per tar, or one per segment for the tars above 6 GiB. For
   the 48 tars today that is 64 jobs: 45 whole tars, 2 segments for SOYBN
   (7.6 GB), 10 for `swissprot_cif_v6.tar` and 7 for `swissprot_pdb_v6.tar`.
   The largest job is the human proteome, 5.2 GB.
2. `build` runs those jobs, eight at a time, each uploading its part file as
   an artifact. Jobs that fail can be re-run on their own.
3. `publish` downloads the parts, merges, spot checks 20 random members plus 2
   from every tar, writes `index-v6.json` (build time, per-tar name, size,
   Last-Modified, ETag and member count, total rows, SHA-256 of the `.gz`) and
   replaces both files on the `index-v6` release. It does not run unless every
   part was built.

**The first build is manual.** With no `index-v6.json` on the release, the
scheduled run does nothing. Otherwise a build that keeps failing would stream
all 153 GB from EBI again every week. A manual run rebuilds unless `force` is
unticked.

**Reusing an earlier run.** A manual run takes an optional `reuse_run`, the id
of an earlier run of this workflow. `plan` downloads that run's part
artifacts (`gh run download <id> --pattern 'part-*'`, which needs the
`actions: read` permission) and leaves out every job whose part is already
there, complete, and built from the file the server has now (same size, ETag
and Last-Modified). A complete whole-tar part stands in for all segments of
that tar. The reused parts are uploaded again as the `part-reused` artifact of
the new run, so `publish` finds them next to the new ones and a later run can
reuse this one in turn. Artifacts are kept for 7 days. The artifacts of a
run's finished jobs can be downloaded while the run is still going, and they
stay when the run fails or is cancelled.

**Replacing the release assets.** `gh release upload --clobber` deletes the
existing asset before it uploads the new one, so a failed upload would leave
the index URL returning 404. The publish step instead uploads to
`<name>.new`, downloads that back and compares it byte for byte with the file
it built, renames the current asset to `<name>.old`, renames `<name>.new` to
`<name>` and deletes `<name>.old`. The live name is absent only between two
rename calls, never during an upload, and a failed upload or a mismatch stops
the job with the previous asset untouched. If a run dies between the two
renames, the next run puts `<name>.old` back first. The index is replaced
before the summary, so the summary never describes an index that is not
there yet.

Things to know:

* A hosted runner job is limited to 6 hours and the build jobs are set to 350
  minutes. At 1.2 MB/s a 4 GiB segment takes about an hour and the human
  proteome about 72 minutes. Parsing is not the limit: the stream walk runs
  at about 340 MB/s against a local server.
* All 64 jobs together are 153 GB. If eight parallel streams each keep
  1.2 MB/s, that is about four and a half hours. Whether the rate holds with
  eight streams has not been measured.
* A job that hits the time limit starts its tar or segment again from the
  beginning on re-run: resume works within a job, not across jobs.
* Runners have about 14 GB of free disk. The build writes only part files
  (about 35 KB for a 180 MB tar, an estimated 6 MB for a whole Swiss-Prot
  tar).
* GitHub disables scheduled workflows in a repository with no activity for 60
  days. Re-enable it from the Actions tab if that happens.
* A reader holding an older index notices a replaced tar by itself, because
  the ETag of the range response no longer matches, and asks for a refresh.

## Size

The local index of five small tars has 16,390 rows in 557,056 bytes (181,035
gzipped), 34 bytes per row. All 48 tars come to 2,211,936 rows, which projects
to about 77 MB, and roughly 25 to 31 MB gzipped. The gzipped figure is an
estimate from a synthetic index of that size; the real one is known after the
first full build.

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
did not build; list every tar you have parts for. To try the segmented path,
build the same tar in pieces into another folder and merge that:

```
for k in 1 2 3 4; do
  python -m scigantic_alphafold._build build-part $T -o build/local/seg --segment $k/4
done
python -m scigantic_alphafold._build merge build/local/seg -o build/local/seg-index/index-v6.sqlite --only $T
```

The rows come out identical to the whole-tar build. `plan --force` prints the
jobs the workflow would run without starting anything.

The same functions are available from Python:

```python
from scigantic_alphafold._build import build_part, merge_parts, verify_index

build_part("UP000000429_85962_HELPY_v6.tar", "build/local/parts")
build_part("UP000000429_85962_HELPY_v6.tar", "build/local/seg", segment=(2, 4))
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
change its ETag. The segment tests index one tar whole and in many layouts
(boundaries inside member data, exactly on a header, inside a header block,
in the end padding, segments with no header at all) and require the merged
rows to be identical. Tests marked `live` talk to EBI and only run with
`SCIGANTIC_ALPHAFOLD_LIVE=1`; one of them streams the 173 MB HELPY tar and
four resynchronize at offsets deep inside the two Swiss-Prot tars.
