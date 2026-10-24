# Behavior contract

Version 0.1.0 is an educational, single-reference, ungapped aligner. It seeks the minimum Hamming distance among retrieved, retained candidates, not an unconditional best biological alignment.

## Input and entry routes

`python -m short_read_aligner --reference reference.fa --reads reads.fastq --output new.tsv` runs directly from the repository. An optional editable installation adds the `short-read-aligner` command. There are no runtime dependencies outside Python's standard library.

The CLI accepts ASCII plain-text files, one nonempty FASTA record (sequence may wrap), and four-line FASTQ records. Blank lines in FASTA are ignored. FASTQ sequence/quality lines may not wrap. Gzip, multi-contig FASTA, SAM/BAM, paired mapping, and ambiguous IUPAC codes other than N are unsupported. Both sequences are uppercased. FASTA uses the first whitespace-delimited header token as the reference name. FASTQ retains the entire name after `@`; duplicate names and reads remain separate rows, in input order. Names containing tabs are rejected. A `+` separator may repeat the complete read name; if present it must match. Quality length must equal read length; quality characters must be printable ASCII 33–126. Qualities are validated then ignored. An empty FASTQ produces a header-only output.

The public Python API is `AlignContext`, `build_index`, `build_minimizer_index`, `align_read`, `align_all`. Build the index from the **same normalized reference and k**, and the same minimizer window when applicable. The context does not scan an externally supplied index for consistency; changing its parameters after construction is outside the supported contract. Use `None` for an uncapped seed hit list, not zero. Reference and reads must be nonempty ACGTN strings. Integer parameters reject bools, negative budgets, nonpositive k/window/cap/vote thresholds. Reads shorter than k are errors; minimizer reads shorter than k+w−1 are errors. References shorter than the seed span have empty indexes and produce unmapped results for valid reads.

## Defaults and retrieval

| Parameter | Release default | Effect |
|---|---:|---|
| k | 12 | Indexed seed length |
| mismatch budget e | 2 | Maximum accepted Hamming distance |
| hit cap | None | No seed suppression |
| minimum seed votes | 1 | Verify every retrieved candidate |
| minimizer window | None | Full k-mer index and spaced read seeds |
| perfect-match early exit | false | Finish verification and count retrieved best ties |

The full index stores every N-free k-mer at ascending zero-based reference positions. A read uses e+1 spaced seeds per strand; for e>0 the offsets are `floor(i*(L-k)/e)`, i=0,…,e; for e=0 the offset is 0. Duplicate offsets are retained and therefore can contribute repeated votes for short overlapping configurations. The code searches the uppercase read on `+` and its reverse complement on `-`. A seed hit at h with read offset o proposes start h−o; negative starts are discarded. Candidates extending past the reference are rejected during verification.

A cap skips the **whole list** when its length is greater than the cap; a list exactly at the cap survives. Voting counts selected seed occurrences pointing to a `(start, strand)` pair. A minimum above one suppresses candidates. Neither is a biological confidence estimate.

In minimizer mode, both reference and read inspect every window of w consecutive k-mers. They select the leftmost N-free k-mer of minimum base-ordered byte key, deduplicate selected positions, and store reference positions sorted. This is lexicographic ordering for equal-length ACGT k-mers, not a randomized hash. Reference and query must both use this strategy.

## Extension, scoring and ties

“Extend” means compare the full oriented read with the equally long reference window; it does not perform dynamic programming. Score is literal Hamming distance: equal characters cost 0, unequal characters cost 1. N equals N at verification and costs 1 against A/C/G/T, but N-containing seeds are never indexed. There is no gap, clipping, quality, pair or mapping-quality model. Hamming counting stops after exceeding e; such candidates are rejected, so that optimization is neutral.

Candidates are verified by descending votes, with a stable sort preserving discovery order for equal votes. Discovery visits `+` before `-`, seeds in extraction order, and each reference posting in ascending order. The representative is the **first candidate encountered at the minimum score**; ties do not replace it. It need not have the smallest coordinate. A palindrome can count the same coordinate twice, once per strand. `n_best` counts distinct retained `(start, strand)` pairs at the best score. It is globally complete only under a completeness condition in [algorithm.md](algorithm.md). With early exit enabled, a zero-mismatch result always has `n_best=1`, and additional exact hits remain uncounted; this does not establish uniqueness.

`align_read` returns `(start, strand, mismatches, n_best)`. Unmapped is `(None, None, e+1, 0)`; e+1 is a sentinel, not a measured minimum. `align_all` returns a list of `(read_id, start, strand, mismatches, n_best)` and ignores quality values. The CLI streams records and writes the richer TSV below.

## TSV output

| Column | Meaning |
|---|---|
| read_id | Complete FASTQ name after `@` |
| reference | First token of the FASTA header, including on unmapped rows |
| start, end | Zero-based half-open reference interval `[start,end)`, with end=start+read length |
| strand | `+` or `-`; start is always the left reference coordinate, including on `-` |
| mismatches | Accepted full-read Hamming distance, or unmapped sentinel e+1 |
| n_best | Retained best `(start,strand)` count; 0 unmapped |
| tie_scope | `retrieved-candidates`, or `exact-hit-truncated` for an early-exit exact result |

Unmapped start/end/strand are `.`. `tie_scope` describes counting, not genome-wide certainty. There are no MAPQ or CIGAR fields. [The fixture](../examples/expected.tsv) is a downstream interface example; this repository does not compute coverage or call peaks.

Successful CLI runs exit 0. Input/parameter/file errors exit 2 and print an error to stderr. The destination parent must exist; existing outputs and paths colliding with inputs are refused. Records are validated into a temporary sibling before the destination is created, so late malformed input leaves no destination file. The final copy is not a crash-atomic transaction; interrupted output or disk errors may leave a partial destination. No index cache is loaded, serialized, or trusted.

## Assumptions and limits

The release defaults have a full-retrieval guarantee for N-free reads/aligned regions of length **at least 36** at e=2. Shorter accepted reads can have overlapping seeds and missed hits. See [the exact theorem and counterexamples](algorithm.md). Large repetitive references can create many candidates because the default removes the historical cap; this is intentional for transparent correctness, not a production scalability claim. Python string/hash-table indexes can consume substantial memory. No sensitivity, throughput, biological specificity, sequencing-error distribution or clinical claim is established beyond the specified synthetic and historical evidence.
