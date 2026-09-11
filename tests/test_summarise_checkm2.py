"""Tests for paired CheckM2 genetic-code selection and QC provenance."""

from __future__ import annotations

import contextlib
import csv
from decimal import Decimal
import io
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import summarise_checkm2  # noqa: E402


class SummariseCheckM2TestCase(unittest.TestCase):
    """Exercise scientific boundaries and failure cases using native report fields."""

    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def write_report(self, code: int, **overrides: str) -> Path:
        """Write an explicit synthetic report for one translation table."""
        row = {
            "Name": "sample",
            "Completeness": "95" if code == 4 else "80",
            "Contamination": "2" if code == 4 else "1",
            "Coding_Density": "0.9" if code == 4 else "0.8",
            "Average_Gene_Length": "300" if code == 4 else "150",
            "Total_Coding_Sequences": "3000" if code == 4 else "5300",
            "Genome_Size": "3000000",
            "GC_Content": "0.3",
            "Contig_N50": "500000",
            "Translation_Table_Used": str(code),
        }
        row.update(overrides)
        path = self.root / f"gcode{code}.tsv"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t")
            writer.writeheader()
            writer.writerow(row)
        return path

    def run_summary(self, report4: Path, report11: Path) -> dict[str, str]:
        """Run the real summary CLI and read its single output row."""
        output = self.root / "summary.tsv"
        result = summarise_checkm2.main([
            "--accession", "ACC1",
            "--gcode4-report", str(report4),
            "--gcode11-report", str(report11),
            "--output", str(output),
        ])
        self.assertEqual(result, 0)
        with output.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertEqual(reader.fieldnames, list(summarise_checkm2.OUTPUT_COLUMNS))
            rows = list(reader)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["accession"], "ACC1")
        self.assertEqual(rows[0]["Gcode_Rule"], "mean_gene_length_ratio")
        self.assertEqual(rows[0]["Gcode_Length_Ratio_Threshold"], "1.5")
        return rows[0]

    def test_paired_metrics_and_selection_provenance(self) -> None:
        """Retain both reports and make the selected code reproducible."""
        row = self.run_summary(self.write_report(4), self.write_report(11))
        self.assertEqual(row["Completeness_gcode4"], "95")
        self.assertEqual(row["Completeness_gcode11"], "80")
        self.assertEqual(row["Average_Gene_Length_gcode4"], "300")
        self.assertEqual(row["Average_Gene_Length_gcode11"], "150")
        self.assertEqual(row["Gcode_Length_Ratio"], "2")
        self.assertEqual(row["Gcode"], "4")
        self.assertEqual(row["Gcode_Selection_Reason"], "length_ratio_above_threshold")
        self.assertEqual(row["Low_quality"], "false")
        self.assertEqual(row["checkm2_status"], "done")
        self.assertEqual(row["warnings"], "")

    def test_ratio_boundaries_use_unrounded_decimal_measurements(self) -> None:
        """Select code 11 at exactly 1.5, including decimal binary-float traps."""
        cases = [
            ("300", "200", "11"),
            ("299.999999", "200", "11"),
            ("300.000001", "200", "4"),
            ("151.35", "100.9", "11"),
            ("300.000000000000000000001", "200", "4"),
            ("299.999999999999999999999", "200", "11"),
            ("100.9", "100.9", "11"),
            ("150", "300", "11"),
        ]
        for length4, length11, expected in cases:
            with self.subTest(length4=length4, length11=length11):
                row = self.run_summary(
                    self.write_report(4, Average_Gene_Length=length4),
                    self.write_report(11, Average_Gene_Length=length11),
                )
                self.assertEqual(row["Gcode"], expected)
                self.assertEqual(row["Average_Gene_Length_gcode4"], length4)
                self.assertEqual(row["Average_Gene_Length_gcode11"], length11)
                self.assertEqual(
                    Decimal(row["Gcode_Length_Ratio"]) > Decimal("1.5"),
                    expected == "4",
                )
                if expected == "11":
                    self.assertEqual(
                        row["Gcode_Selection_Reason"],
                        "length_ratio_at_or_below_threshold",
                    )

    def test_completeness_does_not_override_length_selection(self) -> None:
        """Even a large completeness advantage cannot override the ratio rule."""
        row = self.run_summary(
            self.write_report(4, Completeness="30"),
            self.write_report(11, Completeness="99"),
        )
        self.assertEqual(row["Gcode"], "4")
        self.assertEqual(row["Low_quality"], "true")
        row = self.run_summary(
            self.write_report(4, Completeness="99", Average_Gene_Length="300"),
            self.write_report(11, Completeness="30", Average_Gene_Length="300"),
        )
        self.assertEqual(row["Gcode"], "11")
        self.assertEqual(row["Low_quality"], "true")

    def test_selected_code_controls_quality_at_exact_boundary(self) -> None:
        """Keep the C - 5T <= 50 criterion with exact decimal arithmetic."""
        for code in (4, 11):
            for completeness, contamination, expected in [
                ("70", "4", "true"),
                ("70.000001", "4", "false"),
                ("50.05", "0.01", "true"),
                ("50.050000000000000000001", "0.01", "false"),
            ]:
                with self.subTest(code=code, completeness=completeness):
                    report4 = self.write_report(4)
                    report11 = self.write_report(
                        11, Average_Gene_Length="300" if code == 11 else "150",
                    )
                    selected = self.write_report(
                        code, Completeness=completeness, Contamination=contamination,
                        Average_Gene_Length="300",
                    )
                    row = self.run_summary(
                        selected if code == 4 else report4,
                        selected if code == 11 else report11,
                    )
                    self.assertEqual(row["Gcode"], str(code))
                    self.assertEqual(row["Low_quality"], expected)

    def test_invalid_report_never_becomes_default_code11(self) -> None:
        """Reject invalid lengths, counts, fractions and translation-table labels."""
        invalid_fields = [
            {"Average_Gene_Length": "0"},
            {"Average_Gene_Length": "-1"},
            {"Average_Gene_Length": "nan"},
            {"Average_Gene_Length": "inf"},
            {"Average_Gene_Length": "NA"},
            {"Completeness": "nan"},
            {"Completeness": "101"},
            {"Contamination": "-1"},
            {"Coding_Density": "1.1"},
            {"Total_Coding_Sequences": "0"},
            {"Total_Coding_Sequences": "1.5"},
            {"Translation_Table_Used": "11"},
            {"Translation_Table_Used": "NA"},
            {"Name": ""},
        ]
        for values in invalid_fields:
            with self.subTest(values=values):
                row = self.run_summary(
                    self.write_report(4, **values), self.write_report(11),
                )
                self.assertEqual(row["Gcode"], "NA")
                self.assertEqual(row["Low_quality"], "NA")
                self.assertEqual(row["Gcode_Length_Ratio"], "NA")
                self.assertEqual(row["Gcode_Selection_Reason"], "invalid_report_pair")
                self.assertEqual(row["checkm2_status"], "failed")
                self.assertIn("checkm2_gcode4_failed", row["warnings"])
                self.assertEqual(row["Completeness_gcode11"], "80")

    def test_missing_empty_and_malformed_reports_preserve_valid_partner(self) -> None:
        """Retain diagnostic metrics from the valid partner of a failed report."""
        report4 = self.write_report(4)
        report11 = self.root / "gcode11.tsv"
        for invalid_text in (
            None, "", "Name\tName\nsample\tsample\n", "Name\tx\nsample\n",
        ):
            with self.subTest(invalid_text=invalid_text):
                if invalid_text is not None:
                    report11.write_text(invalid_text)
                row = self.run_summary(report4, report11)
                self.assertEqual(row["Gcode"], "NA")
                self.assertEqual(row["checkm2_status"], "failed")
                self.assertEqual(row["Completeness_gcode4"], "95")
                self.assertEqual(row["Completeness_gcode11"], "NA")

    def test_mismatched_genomes_and_shared_statistics_fail(self) -> None:
        """Do not compare different genomes even when both reports are valid."""
        for values in (
            {"Name": "other"}, {"Genome_Size": "2999999"}, {"GC_Content": "0.4"},
        ):
            with self.subTest(values=values):
                row = self.run_summary(
                    self.write_report(4), self.write_report(11, **values),
                )
                self.assertEqual(row["Gcode"], "NA")
                self.assertEqual(row["checkm2_status"], "failed")
                self.assertEqual(row["warnings"], "inconsistent_shared_stats")
                self.assertEqual(row["Completeness_gcode4"], "95")
                self.assertEqual(row["Completeness_gcode11"], "80")

    def test_obsolete_rule_option_is_rejected(self) -> None:
        """A removed rule cannot silently run a different scientific method."""
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as error:
                summarise_checkm2.parse_args([
                    "--accession", "ACC1", "--gcode4-report", "four.tsv",
                    "--gcode11-report", "eleven.tsv", "--output", "result.tsv",
                    "--gcode-rule", "strict_delta",
                ])
        self.assertEqual(error.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
