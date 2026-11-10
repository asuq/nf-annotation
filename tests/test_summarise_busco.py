"""Tests for BUSCO summary parsing."""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import summarise_busco


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV file into a header and row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames is not None
        return list(reader.fieldnames or []), list(reader)


def read_tsv_row(path: Path) -> dict[str, str]:
    """Read a single-row TSV file."""
    _, rows = read_tsv(path)
    return rows[0]


class SummariseBuscoTestCase(unittest.TestCase):
    """Cover existing-string, numeric reconstruction, and failure cases."""

    def write_json(self, path: Path, payload: dict[str, object]) -> Path:
        """Write a JSON file and return its path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload), encoding="utf-8")
        return path

    def test_main_preserves_existing_busco_string(self) -> None:
        """Reuse a precomputed one-line summary when it is already present."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            summary = self.write_json(
                tmpdir / "short_summary.json",
                {
                    "results": {
                        "one_line_summary": (
                            "C:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200"
                        )
                    }
                },
            )
            output = tmpdir / "summary.tsv"

            exit_code = summarise_busco.main(
                [
                    "--accession",
                    "ACC1",
                    "--summary",
                    str(summary),
                    "--lineage",
                    "bacillota_odb12",
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            row = read_tsv_row(output)
            self.assertEqual(
                row["BUSCO_bacillota_odb12"],
                "C:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200",
            )
            self.assertEqual(row["busco_status"], "done")
            self.assertEqual(row["warnings"], "")

    def test_main_keeps_both_configured_lineages_in_distinct_columns(self) -> None:
        """Emit stable lineage-specific BUSCO columns for both configured lineages."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            bacillota_summary = self.write_json(
                tmpdir / "bacillota.json",
                {
                    "results": {
                        "one_line_summary": (
                            "C:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200"
                        )
                    }
                },
            )
            mycoplasmatota_summary = self.write_json(
                tmpdir / "mycoplasmatota.json",
                {
                    "results": {
                        "one_line_summary": (
                            "C:97.0%[S:96.0%,D:1.0%],F:2.0%,M:1.0%,n:180"
                        )
                    }
                },
            )
            bacillota_output = tmpdir / "bacillota.tsv"
            mycoplasmatota_output = tmpdir / "mycoplasmatota.tsv"

            self.assertEqual(
                summarise_busco.main(
                    [
                        "--accession",
                        "ACC_BOTH",
                        "--summary",
                        str(bacillota_summary),
                        "--lineage",
                        "bacillota_odb12",
                        "--output",
                        str(bacillota_output),
                    ]
                ),
                0,
            )
            self.assertEqual(
                summarise_busco.main(
                    [
                        "--accession",
                        "ACC_BOTH",
                        "--summary",
                        str(mycoplasmatota_summary),
                        "--lineage",
                        "mycoplasmatota_odb12",
                        "--output",
                        str(mycoplasmatota_output),
                    ]
                ),
                0,
            )

            bacillota_header, bacillota_rows = read_tsv(bacillota_output)
            mycoplasmatota_header, mycoplasmatota_rows = read_tsv(mycoplasmatota_output)

            self.assertEqual(
                bacillota_header,
                [
                    "accession",
                    "lineage",
                    "BUSCO_bacillota_odb12",
                    "busco_status",
                    "warnings",
                ],
            )
            self.assertEqual(
                mycoplasmatota_header,
                [
                    "accession",
                    "lineage",
                    "BUSCO_mycoplasmatota_odb12",
                    "busco_status",
                    "warnings",
                ],
            )
            self.assertNotIn("BUSCO_selected", bacillota_header)
            self.assertNotIn("BUSCO_selected", mycoplasmatota_header)
            self.assertEqual(
                bacillota_rows[0]["BUSCO_bacillota_odb12"],
                "C:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200",
            )
            self.assertEqual(
                mycoplasmatota_rows[0]["BUSCO_mycoplasmatota_odb12"],
                "C:97.0%[S:96.0%,D:1.0%],F:2.0%,M:1.0%,n:180",
            )

    def test_main_builds_busco_string_from_numeric_percentages(self) -> None:
        """Construct a cluster_ani-compatible summary string from numeric JSON values."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            summary = self.write_json(
                tmpdir / "short_summary.json",
                {
                    "results": {
                        "Complete percentage": 96.5,
                        "Single copy percentage": 95.0,
                        "Duplicated percentage": 1.5,
                        "Fragmented percentage": 2.0,
                        "Missing percentage": 1.5,
                        "n_markers": 220,
                    }
                },
            )
            output = tmpdir / "summary.tsv"

            exit_code = summarise_busco.main(
                [
                    "--accession",
                    "ACC2",
                    "--summary",
                    str(summary),
                    "--lineage",
                    "mycoplasmatota_odb12",
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            row = read_tsv_row(output)
            self.assertEqual(
                row["BUSCO_mycoplasmatota_odb12"],
                "C:96.5%[S:95.0%,D:1.5%],F:2.0%,M:1.5%,n:220",
            )
            self.assertEqual(row["busco_status"], "done")

    def test_main_builds_busco_string_from_compact_percentages(self) -> None:
        """Construct a summary string from compact BUSCO C/S/D/F/M percentages."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            summary = self.write_json(
                tmpdir / "short_summary.json",
                {
                    "results": {
                        "C": 98.0,
                        "S": 98.0,
                        "D": 0.0,
                        "F": 1.0,
                        "M": 1.0,
                        "n": 200,
                    }
                },
            )
            output = tmpdir / "summary.tsv"

            exit_code = summarise_busco.main(
                [
                    "--accession",
                    "ACC_COMPACT",
                    "--summary",
                    str(summary),
                    "--lineage",
                    "bacillota_odb12",
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            row = read_tsv_row(output)
            self.assertEqual(
                row["BUSCO_bacillota_odb12"],
                "C:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200",
            )
            self.assertEqual(row["busco_status"], "done")

    def test_main_marks_missing_json_as_failed(self) -> None:
        """Emit an NA row when the BUSCO summary file is unavailable."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            summary = tmpdir / "missing.json"
            output = tmpdir / "summary.tsv"

            exit_code = summarise_busco.main(
                [
                    "--accession",
                    "ACC3",
                    "--summary",
                    str(summary),
                    "--lineage",
                    "bacillota_odb12",
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            row = read_tsv_row(output)
            self.assertEqual(row["BUSCO_bacillota_odb12"], "NA")
            self.assertEqual(row["busco_status"], "failed")
            self.assertEqual(row["warnings"], "busco_summary_failed")

    def test_main_marks_missing_lineage_output_as_failed(self) -> None:
        """Emit an NA row for a missing lineage-specific BUSCO output."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            output = tmpdir / "summary.tsv"

            exit_code = summarise_busco.main(
                [
                    "--accession",
                    "ACC_MISSING",
                    "--summary",
                    str(tmpdir / "missing_mycoplasmatota.json"),
                    "--lineage",
                    "mycoplasmatota_odb12",
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            header, rows = read_tsv(output)
            self.assertEqual(
                header,
                [
                    "accession",
                    "lineage",
                    "BUSCO_mycoplasmatota_odb12",
                    "busco_status",
                    "warnings",
                ],
            )
            self.assertEqual(rows[0]["BUSCO_mycoplasmatota_odb12"], "NA")
            self.assertEqual(rows[0]["busco_status"], "failed")
            self.assertEqual(rows[0]["warnings"], "busco_summary_failed")

    def test_main_builds_busco_string_from_counts(self) -> None:
        """Reconstruct percentages from count-style machine-readable fields."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            summary = self.write_json(
                tmpdir / "short_summary.json",
                {
                    "results": {
                        "Complete": 193,
                        "SingleCopy": 190,
                        "Duplicated": 3,
                        "Fragmented": 4,
                        "Missing": 3,
                        "n_markers_total": 200,
                    }
                },
            )
            output = tmpdir / "summary.tsv"

            summarise_busco.main(
                [
                    "--accession",
                    "ACC4",
                    "--summary",
                    str(summary),
                    "--lineage",
                    "bacillota_odb12",
                    "--output",
                    str(output),
                ]
            )

            row = read_tsv_row(output)
            self.assertEqual(
                row["BUSCO_bacillota_odb12"],
                "C:96.5%[S:95.0%,D:1.5%],F:2.0%,M:1.5%,n:200",
            )

    def test_main_marks_invalid_json_payload_as_failed(self) -> None:
        """Emit a failed row when the summary file is not valid JSON."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            summary = tmpdir / "short_summary.json"
            summary.write_text("{not-json", encoding="utf-8")
            output = tmpdir / "summary.tsv"

            summarise_busco.main(
                [
                    "--accession",
                    "ACC5",
                    "--summary",
                    str(summary),
                    "--lineage",
                    "bacillota_odb12",
                    "--output",
                    str(output),
                ]
            )

            row = read_tsv_row(output)
            self.assertEqual(row["BUSCO_bacillota_odb12"], "NA")
            self.assertEqual(row["busco_status"], "failed")

    def test_main_marks_malformed_machine_readable_summary_as_failed(self) -> None:
        """Reject JSON that contains a malformed BUSCO one-line summary string."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            summary = self.write_json(
                tmpdir / "short_summary.json",
                {
                    "results": {
                        "one_line_summary": "C:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%"
                    }
                },
            )
            output = tmpdir / "summary.tsv"

            exit_code = summarise_busco.main(
                [
                    "--accession",
                    "ACC_BAD",
                    "--summary",
                    str(summary),
                    "--lineage",
                    "bacillota_odb12",
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            row = read_tsv_row(output)
            self.assertEqual(row["BUSCO_bacillota_odb12"], "NA")
            self.assertEqual(row["busco_status"], "failed")
            self.assertEqual(row["warnings"], "busco_summary_failed")

    def test_primary_busco_column_uses_first_configured_lineage(self) -> None:
        """Expose the primary BUSCO column from the first configured lineage."""
        self.assertEqual(
            summarise_busco.primary_busco_column(
                ["mycoplasmatota_odb12", "bacillota_odb12"]
            ),
            "BUSCO_mycoplasmatota_odb12",
        )


if __name__ == "__main__":
    unittest.main()
