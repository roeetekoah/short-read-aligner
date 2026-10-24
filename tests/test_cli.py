import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from short_read_aligner.io import read_fasta, read_fastq

ROOT = Path(__file__).resolve().parents[1]
TEMP_ROOT = ROOT / "output" / "test-temp"
TEMP_ROOT.mkdir(parents=True, exist_ok=True)


class InputAndCliTests(unittest.TestCase):
    def run_cli(self, *args):
        return subprocess.run([sys.executable, "-B", "-m", "short_read_aligner", *map(str, args)], cwd=ROOT, capture_output=True, text=True)

    def test_fixture_matches_oracle_verified_golden_output(self):
        with tempfile.TemporaryDirectory(dir=TEMP_ROOT) as directory:
            output = Path(directory) / "result.tsv"
            result = self.run_cli("--reference", "examples/reference.fa", "--reads", "examples/reads.fastq", "--output", output)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(output.read_bytes(), (ROOT / "examples/expected.tsv").read_bytes())
            self.assertEqual(result.stdout, "")
            self.assertEqual(self.run_cli("--reference", "examples/reference.fa", "--reads", "examples/reads.fastq", "--output", output).returncode, 2)

    def test_no_output_on_late_malformed_record_or_input_collision(self):
        with tempfile.TemporaryDirectory(dir=TEMP_ROOT) as directory:
            directory = Path(directory)
            reads = directory / "reads.fastq"
            reads.write_text("@good\nACGT\n+\nIIII\n@bad\nACGT\n+\nII\n")
            output = directory / "out.tsv"
            result = self.run_cli("--reference", "examples/reference.fa", "--reads", reads, "--output", output, "--k", "1")
            self.assertEqual(result.returncode, 2)
            self.assertFalse(output.exists())
            self.assertEqual(list(directory.glob("*.tmp")), [])
            before = reads.read_bytes()
            result = self.run_cli("--reference", "examples/reference.fa", "--reads", reads, "--output", reads, "--k", "1")
            self.assertEqual(result.returncode, 2)
            self.assertEqual(reads.read_bytes(), before)

    def test_early_exit_labels_truncated_exact_ties(self):
        with tempfile.TemporaryDirectory(dir=TEMP_ROOT) as directory:
            output = Path(directory) / "fast.tsv"
            result = self.run_cli("--reference", "examples/reference.fa", "--reads", "examples/reads.fastq", "--output", output, "--early-exit")
            self.assertEqual(result.returncode, 0, result.stderr)
            with output.open() as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            self.assertTrue(any(row["tie_scope"] == "exact-hit-truncated" for row in rows))

    def test_fasta_validation(self):
        with tempfile.TemporaryDirectory(dir=TEMP_ROOT) as directory:
            path = Path(directory) / "ref.fa"
            for text in ("", "ACGT\n", ">a\nACGT\n>b\nAAAA\n", ">\nACGT\n", ">a\n", ">a\nACGU\n"):
                path.write_text(text)
                with self.assertRaises(ValueError):
                    read_fasta(path)
            path.write_bytes(b">ref description\r\nacg\r\ntn\r\n")
            self.assertEqual(read_fasta(path), ("ref", "ACGTN"))

    def test_fastq_validation(self):
        with tempfile.TemporaryDirectory(dir=TEMP_ROOT) as directory:
            path = Path(directory) / "reads.fastq"
            for text in ("@a\nAC\n+\n", "a\nAC\n+\nII\n", "@a\nAC\n-\nII\n", "@a\nAC\n+b\nII\n", "@a\nAC\n+\nI\n", "@a\nAU\n+\nII\n", "@a\nAC\n+\nI \n", "@a\tx\nAC\n+\nII\n"):
                path.write_text(text)
                with self.assertRaises(ValueError):
                    list(read_fastq(path))
            path.write_bytes(b"@a description\r\nac\r\n+a description\r\nII")
            self.assertEqual(list(read_fastq(path)), [("a description", "AC", "II")])
            path.write_text("")
            self.assertEqual(list(read_fastq(path)), [])
