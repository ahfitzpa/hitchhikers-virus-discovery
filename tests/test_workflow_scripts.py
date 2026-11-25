"""Functional and static validation tests for Phase 1 and Phase 2 scripts.

Run with:
    python3 -m unittest discover -s tests -v
"""

from __future__ import annotations

import ast
import py_compile
import re
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
P1_SCRIPTS = REPO_ROOT / "workflow" / "phase1" / "scripts"
P2_SCRIPTS = REPO_ROOT / "workflow" / "phase2" / "scripts"

for _d in (P1_SCRIPTS, P2_SCRIPTS):
    if str(_d) not in sys.path:
        sys.path.insert(0, str(_d))

import family_taxonomy  # noqa: E402
import filter_best_hit  # noqa: E402
import ictv_prepare  # noqa: E402


class TestPythonSyntaxAndDocs(unittest.TestCase):
    def test_all_python_scripts_compile_cleanly(self) -> None:
        py_files = [
            p for p in REPO_ROOT.rglob("*.py")
            if ".git" not in p.parts and "__pycache__" not in p.parts
        ]
        self.assertGreater(len(py_files), 10)
        for p in py_files:
            py_compile.compile(str(p), doraise=True)

    def test_all_workflow_python_scripts_have_module_docstrings(self) -> None:
        for base in (P1_SCRIPTS, P2_SCRIPTS):
            for p in sorted(base.glob("*.py")):
                tree = ast.parse(p.read_text(encoding="utf-8"))
                doc = ast.get_docstring(tree)
                self.assertTrue(
                    doc and len(doc.strip()) >= 15,
                    f"Missing or trivial module docstring in {p.relative_to(REPO_ROOT)}",
                )

    def test_all_referenced_rule_scripts_exist(self) -> None:
        pattern = re.compile(r'["\'](?:\.\./)?(scripts/[A-Za-z0-9_.-]+\.(?:py|sh))["\']')
        for phase in ("phase1", "phase2"):
            phase_dir = REPO_ROOT / "workflow" / phase
            for smk in phase_dir.rglob("*.smk"):
                text = smk.read_text(encoding="utf-8", errors="ignore")
                for rel_script in pattern.findall(text):
                    expected = phase_dir / rel_script
                    self.assertTrue(
                        expected.is_file(),
                        f"{smk.relative_to(REPO_ROOT)} references missing script {rel_script}",
                    )


class TestCoreScriptLogic(unittest.TestCase):
    def test_filter_best_hit_selects_best_evalue_and_filters_fident(self) -> None:
        rows = [
            # Short 100%-identity fragment with worse evalue should not mask the real alignment
            {"query": "q1", "target": "t_short", "pident": 100.0, "alnlen": 39, "evalue": 1e-3, "fident": 1.0},
            {"query": "q1", "target": "t_real", "pident": 82.5, "alnlen": 210, "evalue": 1e-20, "fident": 0.825},
            {"query": "q2", "target": "t_high", "pident": 95.0, "alnlen": 180, "evalue": 1e-30, "fident": 0.95},
        ]
        best = filter_best_hit.best_hit_per_query(rows)
        self.assertEqual(best["q1"]["target"], "t_real")
        kept = filter_best_hit.filter_by_seq_id(best, min_seq_id=0.70, max_seq_id=0.90)
        self.assertEqual(list(kept.keys()), ["q1"])
        self.assertEqual(kept["q1"]["target"], "t_real")

    def test_family_taxonomy_normalise(self) -> None:
        self.assertEqual(family_taxonomy.normalise("Hepeviridae"), "Hepeviridae")
        self.assertEqual(family_taxonomy.normalise("  "), "")

    def test_ictv_safe_name_strips_unsafe_characters(self) -> None:
        cleaned = ictv_prepare.safe_name("Iflavirus sp. (strain A/1)")
        self.assertNotIn(" ", cleaned)
        self.assertNotIn("/", cleaned)

    def test_no_sql_string_interpolation_in_execute_calls(self) -> None:
        for p in (P1_SCRIPTS / "01.centroid_query.py", P1_SCRIPTS / "align_matrix.py"):
            tree = ast.parse(p.read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == "execute":
                    self.assertGreaterEqual(
                        len(node.args),
                        2,
                        f"Expected parameterized cur.execute(sql, params) in {p.name}:{node.lineno}",
                    )
                    self.assertNotIsInstance(
                        node.args[0],
                        (ast.JoinedStr, ast.BinOp),
                        f"SQL string interpolation detected in {p.name}:{node.lineno}",
                    )

    def test_no_unexpanded_pwd_in_sbatch_directives(self) -> None:
        for sh in REPO_ROOT.rglob("*.sh"):
            if ".git" in sh.parts:
                continue
            for lineno, line in enumerate(sh.read_text(encoding="utf-8", errors="ignore").splitlines(), 1):
                if line.startswith("#SBATCH"):
                    self.assertNotIn("$PWD", line, f"Unexpanded $PWD in #SBATCH at {sh.name}:{lineno}")


if __name__ == "__main__":
    unittest.main()
