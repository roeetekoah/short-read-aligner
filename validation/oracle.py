"""Enumerate every full-length window on both strands; no aligner imports.

N is compared literally, matching the stated Hamming model. This oracle is
deliberately simple and has no seeds, filters, candidate votes or early exit.
"""


def all_hits(reference: str, read: str, budget: int):
    reference, read = reference.upper(), read.upper()
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    reverse = "".join(complement[c] for c in reversed(read))
    hits = []
    for start in range(len(reference) - len(read) + 1):
        window = reference[start:start + len(read)]
        for strand, query in (("+", read), ("-", reverse)):
            distance = sum(query[i] != window[i] for i in range(len(read)))
            if distance <= budget:
                hits.append((start, strand, distance))
    return hits


def best_hits(reference: str, read: str, budget: int):
    hits = all_hits(reference, read, budget)
    minimum = min((hit[2] for hit in hits), default=budget + 1)
    return {hit for hit in hits if hit[2] == minimum}
