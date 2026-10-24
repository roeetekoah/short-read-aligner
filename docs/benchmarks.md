# Fixed benchmarks and one bounded ablation

The benchmark asks how retrieval strategy changes correctness, index storage and compute cost on small controlled inputs. It does not compare biological mapping quality or production aligners. The reference oracle is deliberately simple exhaustive Python, intended for correctness; its elapsed time is context, not a competitive implementation baseline.

## Matched contract and fixed inputs

All methods receive the identical preloaded reference/read strings and use both strands, full-length literal Hamming distance, budget 2, no hit cap and one required vote. No paired-end model, quality weighting, clipping, gaps, unmapped exclusion or multimapper exclusion is applied. The oracle returns all minimum-score `(start,strand)` hits; the aligner returns one representative and a count. Correctness checks compare membership in the oracle optimum and tie counts, so a different equally optimal representative is not an error. Oracle-mappable counts, false-unmapped reads, suboptimal scores, invalid returned mappings and tie-count disagreements are separate metrics.

Inputs are committed synthetic strings from seed 20261011, with an exact-byte [manifest](../benchmarks/data/manifest.json). The generator uses a private random generator and independent oracle, with no course data or genome snippets. Regeneration is checked before release. Every benchmark run verifies the input fingerprints before measurement.

| Workload | Reference | Reads | Purpose |
|---|---:|---:|---|
| random | 2,048 ACGT bases | 160, alternating 36/40 bases | Planted 0/1/2/3 substitutions with deterministic strand flips |
| repeat-rich | 1,024 random + 512 A + 512 ACGT-repeat bases | 160, alternating 36/40 bases | Identical sampling/mutation scheme in a repeat-rich reference |
| seed-disruption | 40 A bases | 5 | Exact repeats, two overlapping-seed failures and one three-error read |

Planted origin is not treated as truth: each read is re-evaluated against **all** reference windows on both strands. For example, a three-substitution planted read can still be mappable elsewhere. This is why the repeat-rich workload has 121 oracle-mappable reads rather than the planted scheme's nominal 120.

## Results and interpretation

The [complete measured table](../evidence/benchmark/summary.md) includes three-repeat timing medians, cold compute, memory and counts. [Raw evidence](../evidence/benchmark/results.json) includes every timing repetition, machine/interpreter details, exact configuration, source fingerprints and index counts. Per-read outputs for every variant/workload are saved beside it.

| Method | Random recovered / oracle | Repeat-rich recovered / oracle | Stress recovered / oracle | Random stored positions | Repeat-rich stored positions |
|---|---:|---:|---:|---:|---:|
| Full k=12 | 120/120 | 121/121 | 4/4 | 2,037 | 2,037 |
| Full k=15 | 109/120 | 115/121 | 2/4 | 2,034 | 2,034 |
| Minimizer k=12,w=4 | 120/120 | 121/121 | 4/4 | 868 | 1,083 |
| Full k=12 + exact exit | 120/120 | 121/121 | 4/4 | 2,037 | 2,037 |

Every returned mapping was in bounds with its reported Hamming score, and all mapped results attained an oracle-optimal score in these workloads. The full-k12 baseline had zero false-unmapped, suboptimal, invalid or miscounted-tie results, as expected under its completeness assumptions. Full-k15 loses 11/120, 6/121 and 2/4 mappable reads despite no cap; the experiment isolates seed length and resulting geometry. The small alignment-time saving on some fixtures comes with observed sensitivity loss. It does not establish a generally favorable tradeoff.

Minimizers store 57.4% fewer positions on the random fixture and 46.8% fewer on repeat-rich input. Separate traced-allocation peaks are also lower on these two fixtures, but building and aligning are slower because the simple implementation scans every minimizer window in Python. Minimizer success here is empirical: k=12,w=4,e=2 has sufficient length 45, which these 36/40-base reads do not meet. The test suite includes a different minimizer failure example. Smaller indexes do not establish universal sensitivity or speed.

The single bounded ablation toggles perfect-match exit on the full-k12 baseline, changing no retrieval parameter. It preserves the observed best score and mapping count. On repeat-rich input it reduces alignment time but miscounts ties for 20 of the 49 mapped reads with multiple oracle-optimal alignments; it also miscounts one stress read. The other ambiguous reads do not exit at a perfect hit and retain their counts. The random fixture has no best-hit ambiguity and little work to save. This is an ambiguity/runtime tradeoff, not evidence that count=1 means uniqueness.

## Measurement limits

- One warm-up per variant, then three untraced single-threaded compute repeats in deterministically shuffled variant order. The saved run was repeated without the test suite running concurrently. Loaded input parsing and disk output are excluded for every method. Build includes context validation/index creation; alignment includes read validation, candidate generation, verification and materialized results. Cold compute is build+alignment measured in each repeat, not a mixture of cached and uncached timings. It excludes process startup and file I/O.
- Python memory is measured separately with `tracemalloc` during a fresh index build and alignment, excluding preloaded input objects and interpreter overhead. It is **peak additional traced Python allocation**, not peak RSS or whole-machine memory. The oracle materializes acceptable-hit lists and best-hit sets; repeat multiplicity increases that storage. Different internal/output representations limit memory comparisons.
- The small inputs make sub-millisecond differences noisy, particularly the five-read stress case. These are correctness/stress fixtures, not representative genome-scale workloads. Larger k has almost the same position count here and does not reproduce the historical chromosome-wide dictionary-key growth. There is no hardware scaling, real-data sensitivity, Bowtie2 speed ranking, confidence interval or sequencing-error-model validation.
- Historical paired-end Bowtie2 results and course-run timing claims are documented only as provenance. Their scoring, pairing, output, source edition and postprocessing do not satisfy this matched contract.
