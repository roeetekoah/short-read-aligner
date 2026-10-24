# Fixed synthetic benchmark

All rows use the same preloaded reads, both strands, literal Hamming budget 2, no hit cap and one required vote. Wall times are medians of three or more untraced runs after warm-up. Memory is a separate peak of additional traced Python allocations, not RSS. Oracle output uses all best hits; aligner output uses one representative plus a count. Timings are implementation measurements, not asymptotic or production claims.

| Workload | Method | Build ms | Align ms | Cold compute ms | Python peak KiB | Mapped / oracle | False unmapped | Wrong score / invalid | Tie count disagreements |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| random | full-k12 | 0.353 | 0.705 | 1.058 | 442.2 | 120/120 | 0 | 0/0 | 0 |
| random | full-k15 | 0.362 | 0.691 | 1.054 | 447.7 | 109/120 | 11 | 0/0 | 0 |
| random | minimizer-k12-w4 | 3.240 | 12.227 | 15.493 | 382.1 | 120/120 | 0 | 0/0 | 0 |
| random | full-k12-early-exit | 0.350 | 0.711 | 1.059 | 442.3 | 120/120 | 0 | 0/0 | 0 |
| random | exhaustive oracle | 0 | 936.643 | 936.643 | 49.1 | 120/120 | 0 | 0/0 | 0 |
| repeat-rich | full-k12 | 0.311 | 32.934 | 33.252 | 405.0 | 121/121 | 0 | 0/0 | 0 |
| repeat-rich | full-k15 | 0.311 | 30.364 | 30.683 | 397.6 | 115/121 | 6 | 0/0 | 0 |
| repeat-rich | minimizer-k12-w4 | 2.974 | 69.171 | 72.160 | 290.9 | 121/121 | 0 | 0/0 | 0 |
| repeat-rich | full-k12-early-exit | 0.311 | 22.752 | 23.073 | 404.7 | 121/121 | 0 | 0/0 | 20 |
| repeat-rich | exhaustive oracle | 0 | 923.720 | 923.720 | 2797.4 | 121/121 | 0 | 0/0 | 0 |
| seed-disruption | full-k12 | 0.012 | 0.074 | 0.085 | 7.1 | 4/4 | 0 | 0/0 | 0 |
| seed-disruption | full-k15 | 0.011 | 0.050 | 0.061 | 6.6 | 2/4 | 2 | 0/0 | 0 |
| seed-disruption | minimizer-k12-w4 | 0.043 | 0.487 | 0.530 | 12.8 | 4/4 | 0 | 0/0 | 0 |
| seed-disruption | full-k12-early-exit | 0.010 | 0.063 | 0.073 | 7.1 | 4/4 | 0 | 0/0 | 1 |
| seed-disruption | exhaustive oracle | 0 | 0.048 | 0.048 | 4.2 | 4/4 | 0 | 0/0 | 0 |
