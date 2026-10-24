# Provenance and changes

This portfolio release derives from Roee Tekoah's computational genomics coursework aligner. It was prepared separately on 11 October 2026. Original coursework sources, data, reports and analysis outputs were preserved locally. This release owns read mapping; downstream coverage, peak finding, ATAC interpretation and the historical Bowtie2 analysis remain outside its boundary.

## Contribution boundaries

| Material | Origin and release treatment |
|---|---|
| Core spaced seeding, k-mer index, both-strand candidates, votes, Hamming verification, minimizers | Derived from the implemented coursework `src/aligner.py`; retained algorithmic structure, with later cleanup and validation |
| Reverse-complement helper and instructional wiring | Original source explicitly marks trivial helpers as “provided”; preserve this contribution distinction |
| Instructional reference implementation | Original `aligner_reference.py` explicitly says it was written by Claude; not included, imported or used as an oracle |
| Exact line authorship of scaffolding, optimization and report prose | Surviving comments/recovery notes show AI-assisted instruction and development; per-line ownership cannot be fully reconstructed. Do not describe every line as solely authored by Roee |
| Package, strict loaders/CLI, oracle, exhaustive evidence, synthetic fixtures, new benchmarks, theorem write-up and release docs | Later AI-assisted portfolio additions made during this preparation, separate from the original course submission |
| Original hg19/ATAC data, assignment/presentation, report PDFs, genome-browser images, caches and BAMs | Supplied/external/course evidence, excluded from the distributable release |

The independent oracle enumerates every full-length reference window and both strands without importing any aligner code. It shares the stated scoring contract, not the candidate-generation algorithm. The tests do not use the course test suite's skip-on-`NotImplementedError` hook.

## Source editions reconciled

The implemented `aligner.py` is the selected algorithm source. `old_align_read.py` is unfinished and excluded. The instructional reference has position-only candidates and premature exact-hit exit and cannot replace the working edition. The long root Part A report is byte-identical to the archived `report/old/part_a_report.tex`; other full/trimmed editions differ. The combined `report/ex1_report.tex` contains downstream material and repeats a broader guarantee than its seed geometry supports. The supplement references a missing report filename, and the original README/benchmark route no longer matches the relocated directory structure. This release's docs and commands replace those routes; historical reports are not silently repaired or copied into the package.

Original selected source SHA-256: `fb25759b42cc0ef4e4c9ba822765683ad9b1e19926370c36c6d74f42496e4e6a`.

## Deliberate behavior changes from the working source

- Defaults: k=12 instead of 15; no hit cap instead of 200; one vote and perfect-match exit off. The latter matches the core default but differs from the old benchmark's fast configuration. The default is complete for ACGT reads/aligned regions of length ≥36 with budget 2.
- Standard-library infinity replaces NumPy's only core use. There are no runtime dependencies.
- Added alphabet/parameter validation and errors for reads shorter than the seed span. The old path could silently fail or concatenate FASTA records; the CLI now validates one record, complete FASTQ structure and matching qualities.
- Preserved candidate vote ordering, first-encountered score ties, literal N scoring, minimizer order and leftmost minimizer ties. Clarified that the selection key is base-ordered bytes rather than randomized hashing.
- Omitted cascading/dedup wrappers from the supported standalone API. `align_all` has no progress printing, while the CLI streams outputs. Invalid input leaves no output; existing destinations are refused.

## Historical performance claims

The original fast-run JSON reports 665,651/677,996 mapped R01 reads and 658,115/677,996 R02 reads, with 35.44+34.36 seconds for alignment. Those are historical mapping counts, not verified sensitivity or new release results. The `final_results_best` JSON is identical to the early-exit run metadata and still names that directory. There is no source-version fingerprint tying every run to a specific historical implementation.

The surviving source always vote-orders candidates; the benchmark's voting toggle controls threshold/metadata, not sorting. Therefore historical voting/no-voting speedups cannot be presented as a reproducible controlled ablation of the current code. Memory labels also conflated end-of-phase RSS with peak RSS; the new evidence measures traced Python allocations separately and names that metric explicitly.

The historical Bowtie2 comparison used paired reads, gapped/scoring behavior and different reporting/postprocessing; its mapping-rate gap cannot isolate the cause of missed mappings or establish biological equivalence. It is not used as the release correctness oracle or performance baseline. The new fixed benchmark compares the same preloaded inputs under the same Hamming/strand/multimapping contract, with stated output and implementation differences.

## Licensing status

No upstream LICENSE or explicit redistribution grant was found for the exercise code, scaffolding, supplied data or reports. No permissive license has been invented during preparation. Original datasets and instructor material are excluded; new synthetic fixtures are labeled as later additions. Publication does not add a permissive license or change those contribution boundaries. See [LICENSE.md](../LICENSE.md).
