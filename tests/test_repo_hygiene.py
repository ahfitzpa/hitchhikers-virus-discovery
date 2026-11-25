"""Fast repository hygiene and portability tests for hitchhikers-virus-discovery.

Run with:
    python3 -m unittest discover -s tests -v
"""

from __future__ import annotations

import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


class TestPublicWorkflowHygiene(unittest.TestCase):
    def test_required_top_level_files_present(self) -> None:
        for rel in (
            "README.md",
            "LICENSE",
            "CITATION.cff",
            ".gitignore",
            "workflow/phase1/Snakefile_phase1",
            "workflow/phase2/Snakefile_phase2",
            "workflow/phase2/setup_databases.sh",
            "workflow/phase2/SETUP.md",
        ):
            self.assertTrue((REPO_ROOT / rel).is_file(), f"Missing {rel}")

    def test_no_private_or_personal_paths_in_text_files(self) -> None:
        skip_suffixes = {".png", ".jpg", ".gif", ".svg", ".pyc"}
        self_path = Path(__file__).resolve()
        for path in REPO_ROOT.rglob("*"):
            if (
                not path.is_file()
                or ".git" in path.parts
                or "__pycache__" in path.parts
                or path.suffix in skip_suffixes
                or path.resolve() == self_path
            ):
                continue
            text = path.read_text(encoding="utf-8", errors="ignore")
            self.assertNotIn("#SBATCH --mail-user=", text, f"Hardcoded SLURM mail-user in {path}")
            self.assertNotIn("/home/people/", text, f"Hardcoded HPC home path in {path}")
            self.assertNotIn("@ucd.ie", text, f"Hardcoded institutional email in {path}")

    def test_no_private_analysis_or_manuscript_directories_bundled(self) -> None:
        for forbidden in (
            "manuscript",
            "phylogenomics",
            "workflow/analysis_novel",
            "workflow/analysis_host_range",
            "workflow/analysis_geography",
            "workflow/analysis_surrogate",
        ):
            self.assertFalse((REPO_ROOT / forbidden).exists(), f"Private directory {forbidden} should not be in public repo")


if __name__ == "__main__":
    unittest.main()
