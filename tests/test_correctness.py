"""Independent exhaustive comparisons and executable scope counterexamples."""

import itertools
import unittest

from short_read_aligner.core import (
    AlignContext, align_all, align_read, build_index, build_minimizer_index,
    _collect_candidates, count_mismatches, reverse_complement,
)
from validation.oracle import all_hits, best_hits


def strings(alphabet, length):
    return ("".join(chars) for chars in itertools.product(alphabet, repeat=length))


def context(ref, k, budget, **kwargs):
    w = kwargs.get("minimizer_w")
    index = build_index(ref, k) if w is None else build_minimizer_index(ref, k, w)
    return AlignContext(index, ref, k=k, max_mismatches=budget, **kwargs)


class CorrectnessTests(unittest.TestCase):
    comparisons = 0
    exhaustive_comparisons = 0

    def compare_complete(self, ref, read, ctx, exhaustive=True):
        # Check all acceptable hits, not only a planted location or the best score.
        expected = set(all_hits(ref, read, ctx.max_mismatches))
        rc = reverse_complement(read)
        retrieved = set()
        for pos, strand in _collect_candidates(read, rc, ctx):
            nm = count_mismatches(read if strand == "+" else rc, ref, pos)
            if nm <= ctx.max_mismatches:
                retrieved.add((pos, strand, int(nm)))
        self.assertEqual(retrieved, expected, (ref, read, ctx))
        best = best_hits(ref, read, ctx.max_mismatches)
        result = align_read(read, ctx)
        if best:
            self.assertIn(result[:3], best, (ref, read, ctx))
            self.assertEqual(result[3], len(best), (ref, read, ctx))
        else:
            self.assertEqual(result, (None, None, ctx.max_mismatches + 1, 0))
        type(self).comparisons += 1
        if exhaustive:
            type(self).exhaustive_comparisons += 1

    def test_exhaustive_acgt_k1(self):
        for ref_len in range(1, 5):
            for ref in strings("ACGT", ref_len):
                for read_len in range(1, min(3, ref_len) + 1):
                    for budget in range(read_len):
                        ctx = context(ref, 1, budget)
                        for read in strings("ACGT", read_len):
                            self.compare_complete(ref, read, ctx)

    def test_exhaustive_binary_longer_seeds(self):
        for k, budget, read_len, ref_lengths in ((2, 1, 4, (4, 5, 6)), (2, 2, 6, (6, 7))):
            for ref_len in ref_lengths:
                for ref in strings("AC", ref_len):
                    ctx = context(ref, k, budget)
                    for read in strings("AC", read_len):
                        self.compare_complete(ref, read, ctx)

    def test_exhaustive_minimizer_sufficient_length(self):
        # e=1, k+w-1=3; length 7 guarantees an exact 3-base run.
        for ref in strings("AC", 8):
            ctx = context(ref, 2, 1, minimizer_w=2)
            for read in strings("AC", 7):
                self.compare_complete(ref, read, ctx)

    def test_overlap_counterexamples_at_historical_read_lengths(self):
        for length, positions in ((36, (12, 23)), (40, (13, 25))):
            ref = "A" * length
            read = list(ref)
            for i in positions:
                read[i] = "C"
            read = "".join(read)
            self.assertEqual(best_hits(ref, read, 2), {(0, "+", 2)})
            self.assertEqual(align_read(read, context(ref, 15, 2, max_hits_per_seed=200)), (None, None, 3, 0))
            self.compare_complete(ref, read, context(ref, 12, 2), exhaustive=False)

    def test_hit_cap_and_vote_filter_can_remove_true_hit(self):
        ref, read = "AAAAAA", "AAAA"
        self.assertTrue(best_hits(ref, read, 1))
        self.assertEqual(align_read(read, context(ref, 2, 1, max_hits_per_seed=1)), (None, None, 2, 0))
        ref, read = "AAAA", "CAAA"
        self.assertEqual(align_read(read, context(ref, 2, 1, min_seed_votes=2)), (None, None, 2, 0))

    def test_minimizer_specific_failure_with_disjoint_regular_seeds(self):
        ref = "A" * 45
        read = "A" * 14 + "C" + "A" * 14 + "C" + "A" * 15
        self.assertEqual(best_hits(ref, read, 2), {(0, "+", 2)})
        self.assertEqual(align_read(read, context(ref, 15, 2)), (0, "+", 2, 1))
        self.assertEqual(align_read(read, context(ref, 15, 2, minimizer_w=10)), (None, None, 3, 0))

    def test_n_is_literal_but_not_seedable(self):
        self.assertEqual(best_hits("NNNN", "NNNN", 0), {(0, "+", 0), (0, "-", 0)})
        self.assertEqual(align_read("NNNN", context("NNNN", 2, 0)), (None, None, 1, 0))
        self.assertEqual(align_read("NAAA", context("NAAA", 1, 1)), (0, "+", 0, 1))
        self.assertEqual(align_read("NAAA", context("CAAA", 1, 1)), (0, "+", 1, 1))

    def test_early_exit_loses_ties_but_preserves_exact_score(self):
        ref, read = "ATAT", "AT"
        full = align_read(read, context(ref, 1, 0))
        fast = align_read(read, context(ref, 1, 0, early_exit=True))
        self.assertEqual(full, (0, "+", 0, 4))
        self.assertEqual(fast, (0, "+", 0, 1))

    def test_traversal_tie_rule_is_not_leftmost(self):
        # Larger vote count can precede a lower coordinate at the same score.
        from short_read_aligner.core import _best_alignment
        ctx = context("AAAAAA", 1, 0)
        self.assertEqual(_best_alignment("AA", "TT", {(0, "+"): 1, (3, "+"): 2}, ctx), (3, "+", 0, 2))

    def test_case_strand_boundaries_and_order(self):
        ctx = context("AACCGT", 1, 0)
        self.assertEqual(align_read("acgg", ctx), (2, "-", 0, 1))
        self.assertEqual(align_read("aaccgt", ctx), (0, "+", 0, 1))
        self.assertEqual(align_read("AACCGTA", ctx), (None, None, 1, 0))
        rows = align_all([("a", "AACC", "IIII"), ("a", "ACGG", "!!!!")], ctx)
        self.assertEqual([row[0] for row in rows], ["a", "a"])
        self.assertEqual(align_all([], ctx), [])

    def test_invalid_parameters_and_short_reads(self):
        for kwargs in ({"k": 0}, {"k": True}, {"max_mismatches": -1}, {"max_hits_per_seed": 0}, {"minimizer_w": 0}, {"min_seed_votes": 0}):
            with self.assertRaises(ValueError):
                AlignContext({}, "ACGT", **kwargs)
        for seq in ("", "ACGU", "AC R"):
            with self.assertRaises(ValueError):
                align_read(seq, context("ACGT", 1, 0))
        with self.assertRaises(ValueError):
            align_read("A", context("ACGT", 2, 0))
        with self.assertRaises(ValueError):
            align_read("AC", context("ACGT", 2, 0, minimizer_w=2))


if __name__ == "__main__":
    unittest.main()
