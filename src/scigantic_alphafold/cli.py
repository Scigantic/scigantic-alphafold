"""Command line: ``scigantic-alphafold get|locate|info|refresh|build-part|merge``."""

from __future__ import annotations

import argparse
import os
import sys
from collections.abc import Sequence

from . import _fetch, _index
from ._errors import AlphaFoldError
from ._version import __version__

PROG = "scigantic-alphafold"


def _fail(message: str) -> None:
    print(f"{PROG}: error: {' '.join(str(message).split())}", file=sys.stderr)


def _get(args: argparse.Namespace) -> int:
    failed = 0
    for id in args.ids:
        try:
            if args.stdout:
                sys.stdout.write(_fetch.fetch(id, args.format))
            else:
                print(_fetch.download(id, args.out, args.format))
        except BrokenPipeError:
            raise
        except (AlphaFoldError, OSError) as exc:
            failed += 1
            _fail(f"{id}: {exc}")
    return 1 if failed else 0


def _locate(args: argparse.Namespace) -> int:
    members = _index.locate(args.id, args.format)
    if not members:
        # Same wording as a failed fetch.
        _index.get_index().resolve(args.id, args.format or "cif")
    print("entry_id\tfmt\ttar\toffset\tlength\turl")
    for m in members:
        print(f"{m.entry_id}\t{m.fmt}\t{m.tar}\t{m.offset}\t{m.length}\t{m.url}")
    return 0


def _info(args: argparse.Namespace) -> int:
    for key, value in _index.index_info().items():
        print(f"{key}: {value}")
    return 0


def _refresh(args: argparse.Namespace) -> int:
    print(_index.refresh_index())
    return 0


def _build_part(args: argparse.Namespace) -> int:
    from ._build import build_part

    try:
        if args.base_url:
            print(build_part(args.tar, args.out, base_url=args.base_url))
        else:
            print(build_part(args.tar, args.out))
    except (RuntimeError, ValueError) as exc:
        # The builder's own errors (a tar that changed mid-read, a name that is not a v6 tar).
        _fail(str(exc))
        return 1
    return 0


def _merge(args: argparse.Namespace) -> int:
    from ._build import merge_parts

    try:
        print(merge_parts(args.parts_dir, args.out))
    except (RuntimeError, ValueError) as exc:
        _fail(str(exc))
        return 1
    return 0


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=PROG,
        description="Read single AlphaFold DB v6 structures with one HTTP range request "
        "against EBI's public bulk tars.",
    )
    parser.add_argument("--version", action="version", version=f"{PROG} {__version__}")
    sub = parser.add_subparsers(dest="command", required=True, metavar="command")

    get = sub.add_parser("get", help="download structures, or print them with --stdout")
    get.add_argument("ids", nargs="+", metavar="id", help="UniProt accession or entry id (P69905, AF-P69905-F1)")
    get.add_argument("-f", "--format", choices=_index.FORMATS, default="cif", help="file format (default: cif)")
    get.add_argument("-o", "--out", default=".", metavar="DIR", help="output directory (default: .)")
    get.add_argument("--stdout", action="store_true", help="write the decompressed text to standard output")
    get.set_defaults(run=_get)

    locate = sub.add_parser("locate", help="show where an accession is stored: tar, offset, length")
    locate.add_argument("id", help="UniProt accession or entry id")
    locate.add_argument("-f", "--format", choices=_index.FORMATS, default=None, help="only this format")
    locate.set_defaults(run=_locate)

    info = sub.add_parser("info", help="describe the index in use")
    info.set_defaults(run=_info)

    refresh = sub.add_parser("refresh", help="download the index again")
    refresh.set_defaults(run=_refresh)

    part = sub.add_parser("build-part", help="index one bulk tar into a <tar>.tsv.gz part file")
    part.add_argument("tar", help="tar name, for example UP000000429_85962_HELPY_v6.tar")
    part.add_argument("-o", "--out", default=".", metavar="DIR", help="directory for the part file (default: .)")
    part.add_argument("--base-url", default=None, help="where the tars are served (default: EBI's v6 directory)")
    part.set_defaults(run=_build_part)

    merge = sub.add_parser("merge", help="merge part files into one SQLite index")
    merge.add_argument("parts_dir", metavar="parts-dir", help="directory holding the <tar>.tsv.gz part files")
    merge.add_argument("-o", "--out", required=True, metavar="FILE", help="index to write, for example index-v6.sqlite")
    merge.set_defaults(run=_merge)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        code: int = args.run(args)
    except BrokenPipeError:
        # The reader of our output went away (`... get --stdout | head`). Point stdout at
        # the null device so the interpreter's exit-time flush does not complain.
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        return 1
    except (AlphaFoldError, OSError) as exc:
        _fail(str(exc))
        return 1
    except KeyboardInterrupt:
        return 130
    return code


if __name__ == "__main__":
    sys.exit(main())
