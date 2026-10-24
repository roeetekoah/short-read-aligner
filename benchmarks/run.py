"""Fixed synthetic Hamming benchmarks with one exact-exit ablation.

Times are single-threaded in-process compute with preloaded inputs. Traced
memory is measured separately so instrumentation does not distort timing.
"""

import argparse
from collections import Counter
from datetime import datetime, timezone
import gc
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import statistics
import sys
import time
import tracemalloc

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from short_read_aligner.core import AlignContext, align_read, build_index, build_minimizer_index
from short_read_aligner.io import read_fasta, read_fastq
from validation.oracle import best_hits

CONFIGS = {
    "full-k12": dict(k=12),
    "full-k15": dict(k=15),
    "minimizer-k12-w4": dict(k=12, minimizer_w=4),
    "full-k12-early-exit": dict(k=12, early_exit=True),
}


def execute(reference, reads, config):
    start = time.perf_counter()
    ctx = AlignContext({}, reference, **config)
    ctx.index = build_index(reference, ctx.k) if ctx.minimizer_w is None else build_minimizer_index(reference, ctx.k, ctx.minimizer_w)
    built = time.perf_counter()
    results = [align_read(read, ctx) for read in reads]
    done = time.perf_counter()
    return results, ctx, built - start, done - built


def assess(reference, reads, results, truths):
    counts = Counter(oracle_mapped=sum(bool(t) for t in truths), total=len(reads))
    for read, (pos, strand, nm, n_best), truth in zip(reads, results, truths):
        if pos is None:
            counts["false_unmapped"] += bool(truth)
            continue
        counts["mapped"] += 1
        query = read if strand == "+" else "".join({"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}[c] for c in reversed(read))
        sound = 0 <= pos <= len(reference) - len(read) and nm == sum(a != b for a, b in zip(query, reference[pos:pos + len(read)])) and nm <= 2
        counts["invalid_mapping"] += not sound
        if truth and (pos, strand, nm) in truth:
            counts["optimal_mapping"] += 1
            counts["tie_count_disagreements"] += n_best != len(truth)
            counts["oracle_ambiguous_mapped"] += len(truth) > 1
        else:
            counts["suboptimal_mapping"] += 1
    return {key: counts[key] for key in ("total", "oracle_mapped", "mapped", "false_unmapped", "invalid_mapping", "optimal_mapping", "suboptimal_mapping", "tie_count_disagreements", "oracle_ambiguous_mapped")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("benchmark-output"))
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if args.repeats < 3:
        parser.error("at least three timing repeats are required")
    manifest = json.loads((ROOT / "benchmarks/data/manifest.json").read_text())
    for rel, expected in manifest["sha256"].items():
        if hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() != expected:
            raise SystemExit(f"input fingerprint mismatch: {rel}")
    record = {
        "prepared_date": "2026-10-11", "generated_at_utc": datetime.now(timezone.utc).isoformat(), "python": platform.python_version(), "implementation": platform.python_implementation(),
        "platform": platform.platform(), "processor": platform.processor(), "logical_cpus": os.cpu_count(),
        "single_thread": True, "repeats": args.repeats, "timing": "preloaded inputs; build and alignment separate; no tracing or file I/O; one warm-up per variant",
        "memory": "separate fresh build+alignment; peak additional traced Python allocations; excludes preloaded inputs/interpreter; NOT peak RSS",
        "data_manifest_sha256": hashlib.sha256((ROOT / "benchmarks/data/manifest.json").read_bytes()).hexdigest(),
        "code_sha256": {rel: hashlib.sha256((ROOT / rel).read_bytes()).hexdigest() for rel in ("short_read_aligner/core.py", "validation/oracle.py", "benchmarks/run.py")},
        "workloads": [],
    }
    args.output.mkdir(parents=True, exist_ok=True)
    lines = ["# Fixed synthetic benchmark", "", "All rows use the same preloaded reads, both strands, literal Hamming budget 2, no hit cap and one required vote. Wall times are medians of three or more untraced runs after warm-up. Memory is a separate peak of additional traced Python allocations, not RSS. Oracle output uses all best hits; aligner output uses one representative plus a count. Timings are implementation measurements, not asymptotic or production claims.", "",
             "| Workload | Method | Build ms | Align ms | Cold compute ms | Python peak KiB | Mapped / oracle | False unmapped | Wrong score / invalid | Tie count disagreements |", "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for workload in manifest["workloads"]:
        name = workload["name"]
        _, reference = read_fasta(ROOT / f"benchmarks/data/{name}.fa")
        reads = [seq for _, seq, _ in read_fastq(ROOT / f"benchmarks/data/{name}.fastq")]
        truths = [best_hits(reference, read, 2) for read in reads]  # oracle warm-up and truth
        oracle_times = []
        for _ in range(args.repeats):
            start = time.perf_counter()
            repeated = [best_hits(reference, read, 2) for read in reads]
            oracle_times.append(time.perf_counter() - start)
            assert repeated == truths
        del repeated
        gc.collect()
        tracemalloc.start()
        memory_truth = [best_hits(reference, read, 2) for read in reads]
        _, oracle_peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        del memory_truth
        defaults = dict(max_mismatches=2, max_hits_per_seed=None, min_seed_votes=1, early_exit=False, minimizer_w=None)
        variants = {method: {"config": defaults | config, "build_seconds": [], "alignment_seconds": [], "cold_compute_seconds": []} for method, config in CONFIGS.items()}
        for config in CONFIGS.values():
            warm = execute(reference, reads, config)
            del warm
        for repeat in range(args.repeats):
            methods = list(CONFIGS)
            random.Random(20261011 + repeat).shuffle(methods)
            for method in methods:
                gc.collect()
                results, ctx, build, align = execute(reference, reads, CONFIGS[method])
                variants[method]["build_seconds"].append(build)
                variants[method]["alignment_seconds"].append(align)
                variants[method]["cold_compute_seconds"].append(build + align)
                counts = assess(reference, reads, results, truths)
                if "correctness" in variants[method]:
                    assert variants[method]["correctness"] == counts
                variants[method]["correctness"] = counts
                del ctx, results
        for method, config in CONFIGS.items():
            gc.collect()
            tracemalloc.start()
            results, ctx, _, _ = execute(reference, reads, config)
            _, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            item = variants[method]
            item["peak_additional_python_bytes"] = peak
            item["index_keys"] = len(ctx.index)
            item["index_positions"] = sum(map(len, ctx.index.values()))
            for field in ("build_seconds", "alignment_seconds", "cold_compute_seconds"):
                item["median_" + field] = statistics.median(item[field])
            # Save outputs for inspection and later comparison; no hidden truth tie selection.
            (args.output / f"{name}-{method}-alignments.json").write_text(json.dumps(results, indent=2) + "\n")
            c = item["correctness"]
            lines.append(f"| {name} | {method} | {1000*item['median_build_seconds']:.3f} | {1000*item['median_alignment_seconds']:.3f} | {1000*item['median_cold_compute_seconds']:.3f} | {peak/1024:.1f} | {c['mapped']}/{c['oracle_mapped']} | {c['false_unmapped']} | {c['suboptimal_mapping']}/{c['invalid_mapping']} | {c['tie_count_disagreements']} |")
            del results, ctx
        oracle = {"alignment_seconds": oracle_times, "median_alignment_seconds": statistics.median(oracle_times), "peak_additional_python_bytes": oracle_peak}
        lines.append(f"| {name} | exhaustive oracle | 0 | {1000*oracle['median_alignment_seconds']:.3f} | {1000*oracle['median_alignment_seconds']:.3f} | {oracle_peak/1024:.1f} | {sum(bool(t) for t in truths)}/{sum(bool(t) for t in truths)} | 0 | 0/0 | 0 |")
        record["workloads"].append({"name": name, "reference_bases": len(reference), "reads": len(reads), "oracle": oracle, "variants": variants})
        # Default completeness must hold on these ACGT reads of length >=36.
        baseline = variants["full-k12"]["correctness"]
        assert not any(baseline[key] for key in ("false_unmapped", "invalid_mapping", "suboptimal_mapping", "tie_count_disagreements"))
        for item in variants.values():
            assert item["correctness"]["invalid_mapping"] == 0
        print(f"Verified {name}: {len(reads)} reads", flush=True)
    (args.output / "results.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    (args.output / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Results: {args.output}")


if __name__ == "__main__":
    main()
