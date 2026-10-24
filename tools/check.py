"""Run tests and record the actual exhaustive-comparison count."""

import json
from pathlib import Path
import platform
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

suite = unittest.defaultTestLoader.discover(str(ROOT / "tests"), top_level_dir=str(ROOT))
result = unittest.TextTestRunner(verbosity=2).run(suite)
from tests.test_correctness import CorrectnessTests

record = {
    "python": platform.python_version(),
    "tests": result.testsRun,
    "failures": len(result.failures),
    "errors": len(result.errors),
    "skipped": len(result.skipped),
    "oracle_comparisons": CorrectnessTests.comparisons,
    "exhaustive_comparisons": CorrectnessTests.exhaustive_comparisons,
    "success": result.wasSuccessful(),
}
if len(sys.argv) == 2:
    path = Path(sys.argv[1])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
print(json.dumps(record))
raise SystemExit(0 if result.wasSuccessful() else 1)
