"""Unit and synthetic-fixture integration tests for Phase 1 & Phase 2 core scripts.

Runs cleanly under both stdlib unittest and pytest:
    python3 -m unittest discover -s tests -v
    pytest -v
"""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
P1_SCRIPTS = REPO_ROOT / "workflow" / "phase1" / "scripts"
P2_SCRIPTS = REPO_ROOT / "workflow" / "phase2" / "scripts"

for _d in (P1_SCRIPTS, P2_SCRIPTS):
    if str(_d) not in sys.path:
        sys.path.insert(0, str(_d))

import collect_stage_counts  # noqa: E402
import extract_rdrp_domains  # noqa: E402
import family_taxonomy  # noqa: E402
import filter_best_hit  # noqa: E402

try:
    import convert_manifest  # noqa: E402
    import merge_taxonomy_calls  # noqa: E402
    import phase1_intake  # noqa: E402
    HAS_PANDAS_BIO = True
except ModuleNotFoundError:
    convert_manifest = None
    merge_taxonomy_calls = None
    phase1_intake = None
    HAS_PANDAS_BIO = False


class TestFilterBestHitCLI(unittest.TestCase):
    def test_end_to_end_cli_filters_and_drops_fident(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            raw = Path(td) / "hits_raw.tsv"
            out = Path(td) / "hits_filtered.tsv"
            raw.write_text(
                "q1\tt_short\t100.0\t39\t1e-3\t1.0\n"
                "q1\tt_real\t82.5\t210\t1e-20\t0.825\n"
                "q2\tt_high\t95.0\t180\t1e-30\t0.95\n",
                encoding="utf-8",
            )
            rc = filter_best_hit.main([
                "--hits-raw", str(raw),
                "--min-seq-id", "0.70",
                "--max-seq-id", "0.90",
                "--out", str(out),
                "--log-level", "ERROR",
            ])
            self.assertEqual(rc, 0)
            lines = out.read_text(encoding="utf-8").strip().splitlines()
            self.assertEqual(len(lines), 1)
            self.assertEqual(lines[0].split("\t"), ["q1", "t_real", "82.5", "210", "1e-20"])

    def test_invalid_seq_id_window_exits_nonzero(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            raw = Path(td) / "empty.tsv"
            out = Path(td) / "out.tsv"
            raw.write_text("", encoding="utf-8")
            with self.assertRaises(SystemExit) as ctx:
                filter_best_hit.main([
                    "--hits-raw", str(raw),
                    "--min-seq-id", "0.90",
                    "--max-seq-id", "0.50",
                    "--out", str(out),
                ])
            self.assertNotEqual(ctx.exception.code, 0)

    def test_malformed_row_raises_value_error(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            bad = Path(td) / "bad.tsv"
            bad.write_text("q1\tt1\t80.0\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                filter_best_hit.read_hits_raw(bad)


class TestFamilyTaxonomyReconcile(unittest.TestCase):
    def test_agree_single_tool_and_unresolved(self) -> None:
        prio = ["genomad", "rdrpcatch", "kraken2"]
        self.assertEqual(family_taxonomy.reconcile({}, prio), ("", "unresolved", ""))
        self.assertEqual(
            family_taxonomy.reconcile({"genomad": "Iflaviridae"}, prio),
            ("Iflaviridae", "single_tool", ""),
        )
        self.assertEqual(
            family_taxonomy.reconcile({"genomad": "Iflaviridae", "kraken2": "Iflaviridae"}, prio),
            ("Iflaviridae", "agree", ""),
        )


class TestCollectStageCounts(unittest.TestCase):
    def test_coverage_counts_missing_and_malformed_rows(self) -> None:
        self.assertEqual(collect_stage_counts.coverage_counts("/nonexistent/path.tsv", 500), (-1, -1))
        with tempfile.TemporaryDirectory() as td:
            cov = Path(td) / "cov.tsv"
            cov.write_text(
                "contig_id\tretained\tcontig_length\n"
                "c1\t1\t1200\n"
                "c2\t1\t300\n"
                "c3\t1\tnot_an_int\n"
                "c4\t0\t5000\n",
                encoding="utf-8",
            )
            retained, long_enough = collect_stage_counts.coverage_counts(str(cov), 500)
            self.assertEqual(retained, 3)
            self.assertEqual(long_enough, 1)


class TestExtractRdRpDomains(unittest.TestCase):
    def test_find_col_and_iter_fasta(self) -> None:
        hdr = ["Contig_name", "RdRp_from(AA)", "RdRp_to(AA)"]
        self.assertEqual(extract_rdrp_domains.find_col(hdr, "RdRp_from"), "RdRp_from(AA)")
        self.assertIsNone(extract_rdrp_domains.find_col(hdr, "nonexistent_col"))

        with tempfile.TemporaryDirectory() as td:
            fa = Path(td) / "test.faa"
            fa.write_text(">seq1 extra\nMAAA\nTTTT\n>seq2\nMGGG\n", encoding="utf-8")
            recs = list(extract_rdrp_domains.iter_fasta(str(fa)))
            self.assertEqual(recs, [("seq1", "MAAATTTT"), ("seq2", "MGGG")])


@unittest.skipUnless(HAS_PANDAS_BIO, "Requires pandas and biopython (installed in CI / tax_core env)")
class TestPhase2PandasScripts(unittest.TestCase):
    def test_merge_taxonomy_calls_collect_and_annotate(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            contigs = base / "contigs.fasta"
            contigs.write_text(">c1\nATGAAATAA\n>c2\nATGCCCTAA\n", encoding="utf-8")

            ktax = base / "ktaxonomy.tsv"
            ktax.write_text("1\t|\t1\t|\tR\t|\t0\t|\troot\n10\t|\t1\t|\tF\t|\t1\t|\tIflaviridae\n", encoding="utf-8")
            k2 = base / "k2.tsv"
            k2.write_text("C\tc1\tIflaviridae (taxid 10)\t9\t10:5\n", encoding="utf-8")

            genomad = base / "genomad.tsv"
            genomad.write_text(
                "seq_name\tvirus_score\tn_hallmarks\ttaxonomy\n"
                "c1\t0.95\t2\tViruses;Riboviria;Iflaviridae\n",
                encoding="utf-8",
            )
            rdrp = base / "rdrp.tsv"
            rdrp.write_text("", encoding="utf-8")

            out_tsv = base / "tool_calls.tsv"
            out_fa = base / "viral.fasta"
            rc = merge_taxonomy_calls.main([
                "--log-level", "ERROR",
                "collect",
                "--contigs", str(contigs),
                "--kraken2", str(k2),
                "--ktaxonomy", str(ktax),
                "--genomad", str(genomad),
                "--rdrpcatch", str(rdrp),
                "--out-tsv", str(out_tsv),
                "--out-fasta", str(out_fa),
            ])
            self.assertEqual(rc, 0)
            self.assertIn("Iflaviridae", out_tsv.read_text(encoding="utf-8"))
            self.assertIn(">c1", out_fa.read_text(encoding="utf-8"))
            self.assertNotIn(">c2", out_fa.read_text(encoding="utf-8"))

    def test_phase1_intake_flag_run_blocks_single_end(self) -> None:
        sugg, flags = phase1_intake.flag_run({"librarylayout": "SINGLE", "platform": "ILLUMINA", "scientificname": "Ips confusus"})
        self.assertEqual(sugg, "N")
        self.assertIn("single-end", flags)

    def test_convert_manifest_infer_layout(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            d = Path(td)
            (d / "SRR1_1.fastq.gz").write_text("", encoding="utf-8")
            self.assertEqual(convert_manifest.infer_layout("SRR1", "Iflaviridae", [d]), "single")
            (d / "SRR1_2.fastq.gz").write_text("", encoding="utf-8")
            self.assertEqual(convert_manifest.infer_layout("SRR1", "Iflaviridae", [d]), "paired")


if __name__ == "__main__":
    unittest.main()
