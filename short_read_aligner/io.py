"""Strict, plain-text single-record FASTA and four-line FASTQ input."""

from pathlib import Path
from .core import normalize_dna


def read_fasta(path: str | Path) -> tuple[str, str]:
    name = None
    parts = []
    with open(path, encoding="ascii") as stream:
        for number, line in enumerate(stream, 1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    raise ValueError("FASTA must contain exactly one record")
                name = line[1:].split()[0] if line[1:].split() else ""
                if not name:
                    raise ValueError("FASTA header must contain a record name")
            else:
                if name is None:
                    raise ValueError(f"FASTA line {number}: sequence before header")
                parts.append(line)
    if name is None:
        raise ValueError("FASTA must contain exactly one record")
    return name, normalize_dna("".join(parts), "reference")


def read_fastq(path: str | Path):
    with open(path, encoding="ascii") as stream:
        record = 0
        while True:
            name = stream.readline()
            if not name:
                return
            record += 1
            seq, plus, qual = (stream.readline() for _ in range(3))
            if not seq or not plus or not qual:
                raise ValueError(f"FASTQ record {record}: incomplete four-line record")
            name, seq, plus, qual = (s.rstrip("\r\n") for s in (name, seq, plus, qual))
            if not name.startswith("@") or not name[1:] or "\t" in name:
                raise ValueError(f"FASTQ record {record}: invalid read name")
            if not plus.startswith("+") or (plus[1:] and plus[1:] != name[1:]):
                raise ValueError(f"FASTQ record {record}: invalid separator or repeated name")
            seq = normalize_dna(seq, f"FASTQ record {record}")
            if len(seq) != len(qual) or any(not 33 <= ord(c) <= 126 for c in qual):
                raise ValueError(f"FASTQ record {record}: invalid quality length or characters")
            yield name[1:], seq, qual
