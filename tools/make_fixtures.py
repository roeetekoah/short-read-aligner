"""Deterministically generate synthetic inputs and independent golden fixture.

Run --check to compare existing files without rewriting them.
No coursework data, genome downloads or aligner imports.
"""

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import random
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from validation.oracle import best_hits


def reverse(seq):
    return "".join({"A": "T", "C": "G", "G": "C", "T": "A"}[c] for c in reversed(seq))


def mutate(seq, positions):
    result = list(seq)
    for i in positions:
        result[i] = "ACGT"[("ACGT".index(result[i]) + 1) % 4]
    return "".join(result)


def fasta(name, reference):
    return f">{name}\n" + "\n".join(reference[i:i + 70] for i in range(0, len(reference), 70)) + "\n"


def fastq(reads):
    return "".join(f"@{rid}\n{seq}\n+\n{'I' * len(seq)}\n" for rid, seq in reads)


def generated_files():
    rng = random.Random(20261011)
    dna = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
    repeat = dna(40)
    reference = dna(180) + "N" * 12 + repeat + "N" * 12 + repeat
    reads = [("forward", reference[20:60]), ("reverse", reverse(reference[90:126])),
             ("one_substitution", mutate(reference[30:70], [7])),
             ("two_substitutions", mutate(reference[110:146], [4, 25])),
             ("over_budget", mutate(reference[50:90], [0, 15, 30])),
             ("repeat_tie", repeat), ("unknown_bases", "N" * 36)]
    files = {"examples/reference.fa": fasta("synthetic", reference), "examples/reads.fastq": fastq(reads)}
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(["read_id", "reference", "start", "end", "strand", "mismatches", "n_best", "tie_scope"])
    for rid, seq in reads:
        hits = best_hits(reference, seq, 2)
        if hits:
            pos, strand, nm = sorted(hits)[0]
            writer.writerow([rid, "synthetic", pos, pos + len(seq), strand, nm, len(hits), "retrieved-candidates"])
        else:
            writer.writerow([rid, "synthetic", ".", ".", ".", 3, 0, "retrieved-candidates"])
    files["examples/expected.tsv"] = stream.getvalue()
    workloads = []
    for name, ref in (("random", dna(2048)), ("repeat-rich", dna(1024) + "A" * 512 + "ACGT" * 128)):
        bench_reads = []
        for i in range(160):
            length = (36, 40)[i % 2]
            start = rng.randrange(len(ref) - length + 1)
            seq = ref[start:start + length]
            seq = mutate(seq, rng.sample(range(length), i % 4))
            if i % 3 == 0:
                seq = reverse(seq)
            bench_reads.append((f"{name}_{i:03d}", seq))
        files[f"benchmarks/data/{name}.fa"] = fasta(name, ref)
        files[f"benchmarks/data/{name}.fastq"] = fastq(bench_reads)
        workloads.append({"name": name, "reference_bases": len(ref), "reads": len(bench_reads), "lengths": [36, 40], "mutation_counts": [0, 1, 2, 3]})
    stress = [("exact36", "A" * 36), ("exact40", "A" * 40),
              ("overlap36", mutate("A" * 36, [12, 23])), ("overlap40", mutate("A" * 40, [13, 25])),
              ("three_errors", mutate("A" * 40, [0, 15, 30]))]
    files["benchmarks/data/seed-disruption.fa"] = fasta("seed-disruption", "A" * 40)
    files["benchmarks/data/seed-disruption.fastq"] = fastq(stress)
    workloads.append({"name": "seed-disruption", "reference_bases": 40, "reads": len(stress), "lengths": [36, 40], "mutation_counts": [0, 2, 3]})
    manifest = {"generator_seed": 20261011, "origin": "new synthetic portfolio data; no course reads or genome data", "workloads": workloads,
                "sha256": {path: hashlib.sha256(text.encode("utf-8")).hexdigest() for path, text in sorted(files.items())}}
    files["benchmarks/data/manifest.json"] = json.dumps(manifest, indent=2) + "\n"
    return files


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    for relative, content in generated_files().items():
        path = ROOT / relative
        if args.check:
            if not path.exists() or path.read_bytes() != content.encode("utf-8"):
                raise SystemExit(f"fixture differs: {relative}")
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content.encode("utf-8"))
    print("Synthetic fixtures verified" if args.check else "Synthetic fixtures generated")
