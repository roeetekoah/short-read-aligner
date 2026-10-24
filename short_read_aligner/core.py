"""Ungapped seed-and-extend mapping, derived from the coursework aligner.

See docs/provenance.md for inherited scaffolding and later portfolio additions.
The index must be built from the same reference, k and minimizer window as ctx.
"""

from collections import defaultdict
from dataclasses import dataclass
from math import inf
from typing import Iterable

K = 12
MAX_MISMATCHES = 2
MAX_HITS_PER_SEED = None
MINIMIZER_W = 10
_COMPLEMENT_TABLE = str.maketrans("ACGTN", "TGCAN")
_BASE_TO_BITS = str.maketrans({"A": "\x00", "C": "\x01", "G": "\x02", "T": "\x03", "N": "\x00"})


def _integer(name: str, value: int, minimum: int) -> None:
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def normalize_dna(seq: str, name: str = "sequence") -> str:
    seq = seq.upper()
    if not seq or set(seq) - set("ACGTN"):
        raise ValueError(f"{name} must be nonempty DNA containing only A, C, G, T, N")
    return seq


@dataclass
class AlignContext:
    index: dict[str, list[int]]
    ref: str
    k: int = K
    max_mismatches: int = MAX_MISMATCHES
    max_hits_per_seed: int | None = MAX_HITS_PER_SEED
    early_exit: bool = False
    minimizer_w: int | None = None
    min_seed_votes: int = 1

    def __post_init__(self):
        self.ref = normalize_dna(self.ref, "reference")
        _integer("k", self.k, 1)
        _integer("max_mismatches", self.max_mismatches, 0)
        _integer("min_seed_votes", self.min_seed_votes, 1)
        if self.max_hits_per_seed is not None:
            _integer("max_hits_per_seed", self.max_hits_per_seed, 1)
        if self.minimizer_w is not None:
            _integer("minimizer_w", self.minimizer_w, 1)


def reverse_complement(seq: str) -> str:
    return seq.upper().translate(_COMPLEMENT_TABLE)[::-1]


def valid_kmer(kmer: str) -> bool:
    return "N" not in kmer


def build_index(reference: str, k: int = K) -> dict[str, list[int]]:
    reference = normalize_dna(reference, "reference")
    _integer("k", k, 1)
    index = defaultdict(list)
    for i in range(len(reference) - k + 1):
        kmer = reference[i:i + k]
        if valid_kmer(kmer):
            index[kmer].append(i)
    return dict(index)


def count_mismatches(read: str, ref: str, start: int, max_mm=None) -> float:
    if start < 0 or start + len(read) > len(ref):
        return inf
    nm = 0
    for a, b in zip(read, ref[start:start + len(read)]):
        if a != b:
            nm += 1
            if max_mm is not None and nm > max_mm:
                return nm
    return nm


def extract_seeds(read: str, ctx: AlignContext) -> list[tuple[int, str]]:
    n_seeds = ctx.max_mismatches + 1
    offsets = [0] if n_seeds == 1 else [
        i * (len(read) - ctx.k) // (n_seeds - 1) for i in range(n_seeds)
    ]
    return [(o, read[o:o + ctx.k]) for o in offsets]


def extract_hits(seed: str, ctx: AlignContext) -> list[int]:
    hits = ctx.index.get(seed, [])
    if ctx.max_hits_per_seed is not None and len(hits) > ctx.max_hits_per_seed:
        return []
    return hits


def _minimizer_key(kmer: str) -> int:
    # Base-ordered byte encoding, not a randomized hash or packed 2-bit hash.
    return int.from_bytes(kmer.translate(_BASE_TO_BITS).encode("latin-1"), "big")


def build_minimizer_index(reference: str, k: int = K, w: int = MINIMIZER_W) -> dict[str, list[int]]:
    reference = normalize_dna(reference, "reference")
    _integer("k", k, 1)
    _integer("w", w, 1)
    index = defaultdict(set)
    for w_start in range(len(reference) - k - w + 2):
        best_key, best_pos, best_kmer = None, -1, None
        for i in range(w_start, w_start + w):
            kmer = reference[i:i + k]
            if not valid_kmer(kmer):
                continue
            key = _minimizer_key(kmer)
            if best_key is None or key < best_key:
                best_key, best_pos, best_kmer = key, i, kmer
        if best_kmer is not None:
            index[best_kmer].add(best_pos)
    return {kmer: sorted(positions) for kmer, positions in index.items()}


def extract_minimizers(seq: str, ctx: AlignContext) -> list[tuple[int, str]]:
    seen_offsets = set()
    out = []
    for w_start in range(len(seq) - ctx.k - ctx.minimizer_w + 2):
        best_key, best_off, best_kmer = None, -1, None
        for i in range(w_start, w_start + ctx.minimizer_w):
            kmer = seq[i:i + ctx.k]
            if not valid_kmer(kmer):
                continue
            key = _minimizer_key(kmer)
            if best_key is None or key < best_key:
                best_key, best_off, best_kmer = key, i, kmer
        if best_kmer is None or best_off in seen_offsets:
            continue
        seen_offsets.add(best_off)
        out.append((best_off, best_kmer))
    return out


def _collect_candidates(read: str, rc: str, ctx: AlignContext) -> dict[tuple[int, str], int]:
    candidates = {}
    for strand, seq in (("+", read), ("-", rc)):
        seeds = extract_minimizers(seq, ctx) if ctx.minimizer_w is not None else extract_seeds(seq, ctx)
        for offset, seed in seeds:
            for hit in extract_hits(seed, ctx):
                start = hit - offset
                if start >= 0:
                    key = (start, strand)
                    candidates[key] = candidates.get(key, 0) + 1
    return candidates


def _best_alignment(read: str, rc: str, candidates: dict[tuple[int, str], int], ctx: AlignContext):
    best_pos, best_strand, best_nm, n_best = None, None, ctx.max_mismatches + 1, 0
    # Stable sort retains discovery order for equal vote counts.
    for (pos, strand), votes in sorted(candidates.items(), key=lambda item: -item[1]):
        if votes < ctx.min_seed_votes:
            continue
        nm = count_mismatches(read if strand == "+" else rc, ctx.ref, pos, ctx.max_mismatches)
        if nm > ctx.max_mismatches:
            continue
        if nm < best_nm:
            best_pos, best_strand, best_nm, n_best = pos, strand, int(nm), 1
            if ctx.early_exit and best_nm == 0:
                break
        elif nm == best_nm:
            n_best += 1
    return best_pos, best_strand, best_nm, n_best


def align_read(read: str, ctx: AlignContext) -> tuple[int | None, str | None, int, int]:
    """Return (0-based left start, strand, Hamming distance, retrieved-best count).

    Unmapped: (None, None, budget + 1, 0). Completeness is conditional;
    see docs/algorithm.md. N compares literally but cannot seed.
    """
    read = normalize_dna(read, "read")
    minimum = ctx.k + (ctx.minimizer_w - 1 if ctx.minimizer_w is not None else 0)
    if len(read) < minimum:
        raise ValueError(f"read length {len(read)} is below seed span {minimum}")
    rc = reverse_complement(read)
    return _best_alignment(read, rc, _collect_candidates(read, rc, ctx), ctx)


def align_all(reads: Iterable[tuple[str, str, str]], ctx: AlignContext):
    """Return one five-field row per input record; qualities are ignored."""
    return [(rid, *align_read(seq, ctx)) for rid, seq, _quality in reads]
