# Seed-and-extend design and its actual guarantees

The index turns an exact substring query into a list of possible locations. A seed starting at offset o in an oriented read and found at h in the reference proposes the read start h−o. Several seeds may propose the same start, so candidates are deduplicated by `(start,strand)` and accumulate votes. Verification scans the full read, rejecting boundaries and over-budget Hamming distance. This trades exhaustive reference-window comparison for index storage and a smaller candidate set.

## Spaced-seed theorem

Assume integer k≥1, e≥0, L≥k; the index was built with the same normalized reference and k; read and aligned reference region contain only ACGT; no hit cap removes seeds; minimum votes is one. The implementation's e+1 length-k seed intervals give a universal guarantee of retrieving **every** full-length alignment with at most e substitutions exactly when:

**L ≥ (e+1)k.**

For e=0, offset 0 suffices. For e>0, adjacent offsets `floor(i*(L-k)/e)` differ by at least `floor((L-k)/e)`. All e+1 intervals are disjoint precisely when L−k≥ek. With disjoint intervals, e errors cannot hit all e+1 seeds. At least one seed is exact and indexed at its corresponding reference coordinate; subtracting its read offset generates the true start and strand. This reasoning applies to every acceptable alignment independently, including reverse-complement and repeat hits.

Conversely, if two seed intervals overlap, one substitution in their intersection destroys both, and one in each remaining interval destroys all seeds using at most e substitutions. Set the reference to A repeated L times and substitute C at those positions. Every forward seed now contains C; every reverse-complement seed contains T/G and cannot match the all-A index. The read is unmapped despite a valid Hamming alignment. Thus overlapping e+1 equal-length seeds do **not** merely weaken a universal guarantee by a small amount: they invalidate it. Particular reads may still map correctly.

| Read length L | k | e | Actual offsets | Substitution positions in A^L | Result |
|---:|---:|---:|---|---|---|
| 36 | 15 | 2 | 0,10,21 | 12,23 | Misses the valid two-mismatch alignment |
| 40 | 15 | 2 | 0,12,25 | 13,25 | Misses the valid two-mismatch alignment |
| 36 or 40 | 12 | 2 | 0,12,24 or 0,14,28 | Same positions | Recovers the alignment without suppression |

Positions are zero-based. These counterexamples are executed in the test suite and represented in the benchmark stress fixture. The theorem was separately reviewed as the requested bounded Xhigh question and checked against an exhaustive oracle; it is a code-specific derivation, not a borrowed biological claim.

For a particular alignment the useful condition is more local: enough selected seeds must be exact, N-free, and retained by the hit cap to reach the minimum vote threshold. Disjoint seeds and m≤e substitutions leave at least e+1−m exact seeds before filtering, but across all m≤e only one is guaranteed. Equal N characters can invalidate retrieval even when literal Hamming distance is zero: N-only exact reads have no usable seeds.

## Minimizer guarantee

Let W=k+w−1, the base span of one window of w k-mers. An identical W-base stretch containing at least one N-free k-mer gives a corresponding minimizer occurrence on both sequences: each examines the same ordered list, chooses the same leftmost minimum, and retains that position after deduplication. This generates the coordinate-consistent start unless its hit list or vote is suppressed.

For ACGT alignments with at most e substitutions, **L≥(e+1)W** is sufficient to force an exact W-base run. With m≤e substitutions, m+1 exact runs contain L−m bases. If every run were at most W−1, then L≤(m+1)W−1. The length bound therefore forces a shared exact window. This bound is sharp for forcing such a run, but is not asserted necessary for actual minimizer retrieval.

At k=15,w=10,e=2, the sufficient length is 72; the original 36/40-base reads do not meet it. A more specific failure exists even when ordinary k-mer seeds are disjoint: reference A^45; read `A*14 + C + A*14 + C + A*15`. Full k-mer seeding retrieves `(0,+)` with distance 2. Minimizer seeding misses it. The read's only exact A^15 is at offset 30, while leftmost ties cause the reference minimizer index to store A^15 only at positions 0–21. Tests execute this example.

The actual selection key is base-ordered byte encoding. A-rich words can be favored. The common random-hash density approximation 2/(w+1) is not a guarantee for this implementation; benchmark evidence reports actual stored positions.

## What verification establishes

Every returned mapping is a full-length in-bounds alignment of the reported strand and Hamming score ≤e. Full verification finds the minimum score and exact tie count among retained candidates. Under complete retrieval, these are also the globally best within-budget score and tie count. Early Hamming rejection is neutral. Perfect-match exit preserves the minimum score among retrieved candidates, but loses exact-match tie counting. A hit cap, vote filter, N-skipping, or incomplete seed geometry can affect retrieval, minimum score and ambiguity. `n_best` is not MAPQ.

The standalone route omits historical cascading and duplicate caching. The former changes the seed layout between passes and can stop at a suboptimal retained alignment when filtering hides a better alignment until a later pass. The latter can cache computation but is not biological duplicate removal. Omitting these modes makes the supported contract easier to test; their course implementations remain intact locally.

## Cost and tradeoffs

For reference length G, full indexing creates G−k+1 potential words, with O(Gk) character work/storage in this Python string-slicing implementation, plus integer posting storage. For fixed k it is often described as linear in G. Larger k reduces chance seed collisions but may increase distinct dictionary keys, memory, and overlap-related misses. With H retained seed hits across both strands, candidate construction is O(H), sorting C distinct candidates is O(C log C), and full verification is O(CL), excluding seed creation. Uncapped repeats can make C large.

The simple minimizer implementation explicitly scans w k-mers per window, so it uses O(Gwk) work before deduplication; it has no rolling hash or monotonic queue. It can store fewer positions but still run slower on small inputs. The benchmark is designed to expose those tradeoffs without asserting that they generalize to large genomes or optimized implementations.
