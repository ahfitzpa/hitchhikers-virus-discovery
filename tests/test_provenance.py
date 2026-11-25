"""Unit tests for the run provenance manifest writer (write_run_manifest.py)."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

P2_SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "phase2" / "scripts"
if str(P2_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(P2_SCRIPTS))

import write_run_manifest  # noqa: E402


class TestRunProvenanceManifest(unittest.TestCase):
    def test_write_manifest_records_sha256_and_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            manifest_csv = base / "manifest.csv"
            manifest_csv.write_text("sample_id,family\nSRR1,Iflaviridae\n", encoding="utf-8")
            out_json = base / "results" / "run_manifest.json"

            rc = write_run_manifest.main([
                "--workflow", "phase2",
                "--input", str(manifest_csv),
                "--output-dir", str(base / "results"),
                "--out", str(out_json),
                "--command", "snakemake -s Snakefile_phase2",
                "--log-level", "ERROR",
            ])
            self.assertEqual(rc, 0)
            data = json.loads(out_json.read_text(encoding="utf-8"))
            self.assertEqual(data["workflow"], "phase2")
            self.assertEqual(data["workflow_version"], write_run_manifest.__version__)
            self.assertEqual(data["command"], "snakemake -s Snakefile_phase2")
            self.assertEqual(
                data["input_sha256"]["manifest.csv"],
                write_run_manifest.sha256_file(manifest_csv),
            )

    def test_missing_input_file_raises(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaises(FileNotFoundError):
                write_run_manifest.build_manifest(
                    workflow="phase2",
                    input_files=[Path(td) / "missing.csv"],
                    output_dir=td,
                )


if __name__ == "__main__":
    unittest.main()
