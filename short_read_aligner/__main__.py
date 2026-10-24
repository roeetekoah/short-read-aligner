"""Command-line entry point; write output only after all records succeed."""

import argparse
import csv
from pathlib import Path
import tempfile

from . import __version__
from .core import AlignContext, align_read, build_index, build_minimizer_index
from .io import read_fasta, read_fastq

HEADER = ["read_id", "reference", "start", "end", "strand", "mismatches", "n_best", "tie_scope"]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="Ungapped Hamming seed-and-extend mapping. See docs/contract.md.")
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--reads", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path, help="New TSV path; existing files are refused")
    parser.add_argument("--k", type=int, default=12)
    parser.add_argument("--max-mismatches", type=int, default=2)
    parser.add_argument("--max-hits-per-seed", type=int, default=None, help="Skip seeds exceeding this count (default: no cap)")
    parser.add_argument("--min-seed-votes", type=int, default=1)
    parser.add_argument("--minimizer-window", type=int, default=None)
    parser.add_argument("--early-exit", action="store_true", help="Stop at first exact hit; exact-hit ties are not counted")
    args = parser.parse_args(argv)
    temporary = None
    try:
        if args.output.resolve() in (args.reference.resolve(), args.reads.resolve()):
            raise ValueError("output must differ from input paths")
        if args.output.exists():
            raise ValueError("output already exists; choose a new path")
        name, reference = read_fasta(args.reference)
        ctx = AlignContext({}, reference, args.k, args.max_mismatches, args.max_hits_per_seed,
                           args.early_exit, args.minimizer_window, args.min_seed_votes)
        ctx.index = build_index(reference, ctx.k) if ctx.minimizer_w is None else build_minimizer_index(reference, ctx.k, ctx.minimizer_w)
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="", suffix=".tsv.tmp",
                                         dir=args.output.parent, delete=False) as stream:
            temporary = Path(stream.name)
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(HEADER)
            for rid, seq, _quality in read_fastq(args.reads):
                pos, strand, nm, n_best = align_read(seq, ctx)
                scope = "exact-hit-truncated" if args.early_exit and pos is not None and nm == 0 else "retrieved-candidates"
                writer.writerow([rid, name, "." if pos is None else pos,
                                 "." if pos is None else pos + len(seq), strand or ".", nm, n_best, scope])
        # Exclusive creation also protects against a destination appearing during the run.
        with args.output.open("xb") as destination, temporary.open("rb") as source:
            import shutil
            shutil.copyfileobj(source, destination)
    except (ValueError, OSError, UnicodeError) as error:
        parser.exit(2, f"error: {error}\n")
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
