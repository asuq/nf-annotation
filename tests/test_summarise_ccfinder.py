"""Tests for CRISPRCasFinder result summarisation."""

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

import summarise_ccfinder


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV file into a header and list of dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames is not None
        return list(reader.fieldnames or []), list(reader)


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    """Read a TSV file into a list of dictionaries."""
    _, rows = read_tsv(path)
    return rows


class SummariseCCFinderTestCase(unittest.TestCase):
    """Cover retained arrays, evidence filtering, and failure handling."""

    def write_json(self, path: Path, payload: dict[str, object]) -> Path:
        """Write a JSON file and return its path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload), encoding="utf-8")
        return path

    def test_main_builds_strain_contig_and_crispr_tables(self) -> None:
        """Ignore evidence-level 1 arrays and aggregate retained arrays correctly."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            result_json = self.write_json(
                tmpdir / "result.json",
                {
                    "Sequences": [
                        {
                            "Id": "contig1",
                            "Length": 1000,
                            "Crisprs": [
                                {
                                    "Name": "c1_a",
                                    "Evidence_Level": 4,
                                    "Start": 1,
                                    "End": 100,
                                    "Spacers": [{}, {}, {}],
                                },
                                {
                                    "Name": "c1_b",
                                    "Evidence_Level": 1,
                                    "Start": 200,
                                    "End": 260,
                                    "Spacers": [{}, {}],
                                },
                            ],
                        },
                        {
                            "Id": "contig2",
                            "Length": 500,
                            "Crisprs": [
                                {
                                    "Name": "c2_a",
                                    "Evidence_Level": 3,
                                    "Start": 50,
                                    "End": 99,
                                    "Spacers": 4,
                                }
                            ],
                        },
                    ]
                },
            )
            outdir = tmpdir / "out"

            exit_code = summarise_ccfinder.main(
                [
                    "--accession",
                    "ACC1",
                    "--result-json",
                    str(result_json),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            strain_header, strain_rows = read_tsv(outdir / "ccfinder_strains.tsv")
            contig_header, contig_rows = read_tsv(outdir / "ccfinder_contigs.tsv")
            crispr_header, crispr_rows = read_tsv(outdir / "ccfinder_crisprs.tsv")

            self.assertEqual(len(strain_rows), 1)
            self.assertEqual(
                strain_header,
                [
                    "accession",
                    "CRISPRS",
                    "SPACERS_SUM",
                    "CRISPR_FRAC",
                    "ccfinder_status",
                    "warnings",
                ],
            )
            self.assertEqual(strain_rows[0]["CRISPRS"], "2")
            self.assertEqual(strain_rows[0]["SPACERS_SUM"], "7")
            self.assertEqual(strain_rows[0]["CRISPR_FRAC"], "0.1")
            self.assertEqual(strain_rows[0]["ccfinder_status"], "done")
            self.assertEqual(strain_rows[0]["warnings"], "")

            self.assertEqual(len(contig_rows), 2)
            self.assertEqual(
                contig_header,
                [
                    "accession",
                    "contig_id",
                    "contig_length",
                    "CRISPRS",
                    "SPACERS_SUM",
                    "CRISPR_FRAC",
                ],
            )
            self.assertEqual(contig_rows[0]["CRISPRS"], "1")
            self.assertEqual(contig_rows[1]["SPACERS_SUM"], "4")

            self.assertEqual(len(crispr_rows), 2)
            self.assertEqual(
                crispr_header,
                [
                    "accession",
                    "contig_id",
                    "crispr_id",
                    "evidence_level",
                    "spacer_count",
                    "start",
                    "end",
                    "crispr_length",
                ],
            )
            self.assertEqual(crispr_rows[0]["crispr_id"], "c1_a")
            self.assertEqual(crispr_rows[1]["crispr_length"], "50")

    def test_main_handles_samples_without_retained_crisprs(self) -> None:
        """Write zero-valued summaries when only evidence-level 1 arrays exist."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            result_json = self.write_json(
                tmpdir / "result.json",
                {
                    "Sequences": [
                        {
                            "Id": "contig1",
                            "Length": 1000,
                            "Crisprs": [
                                {
                                    "Evidence_Level": 1,
                                    "Start": 1,
                                    "End": 100,
                                    "Spacers": 2,
                                }
                            ],
                        }
                    ]
                },
            )
            outdir = tmpdir / "out"

            exit_code = summarise_ccfinder.main(
                [
                    "--accession",
                    "ACC2",
                    "--result-json",
                    str(result_json),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            strain_row = read_tsv_rows(outdir / "ccfinder_strains.tsv")[0]
            contig_rows = read_tsv_rows(outdir / "ccfinder_contigs.tsv")
            crispr_rows = read_tsv_rows(outdir / "ccfinder_crisprs.tsv")
            self.assertEqual(strain_row["CRISPRS"], "0")
            self.assertEqual(strain_row["SPACERS_SUM"], "0")
            self.assertEqual(strain_row["CRISPR_FRAC"], "0")
            self.assertEqual(strain_row["ccfinder_status"], "done")
            self.assertEqual(contig_rows[0]["CRISPRS"], "0")
            self.assertEqual(contig_rows[0]["SPACERS_SUM"], "0")
            self.assertEqual(contig_rows[0]["CRISPR_FRAC"], "0")
            self.assertEqual(crispr_rows, [])

    def test_main_handles_samples_with_no_crispr_arrays(self) -> None:
        """Write zero-valued summaries when the sample has no CRISPR arrays at all."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            result_json = self.write_json(
                tmpdir / "result.json",
                {
                    "Sequences": [
                        {
                            "Id": "contig1",
                            "Length": 750,
                            "Crisprs": [],
                        },
                        {
                            "Id": "contig2",
                            "Length": 250,
                        },
                    ]
                },
            )
            outdir = tmpdir / "out"

            exit_code = summarise_ccfinder.main(
                [
                    "--accession",
                    "ACC_EMPTY",
                    "--result-json",
                    str(result_json),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            strain_row = read_tsv_rows(outdir / "ccfinder_strains.tsv")[0]
            contig_rows = read_tsv_rows(outdir / "ccfinder_contigs.tsv")
            crispr_rows = read_tsv_rows(outdir / "ccfinder_crisprs.tsv")

            self.assertEqual(strain_row["accession"], "ACC_EMPTY")
            self.assertEqual(strain_row["CRISPRS"], "0")
            self.assertEqual(strain_row["SPACERS_SUM"], "0")
            self.assertEqual(strain_row["CRISPR_FRAC"], "0")
            self.assertEqual(strain_row["ccfinder_status"], "done")
            self.assertEqual(strain_row["warnings"], "")
            self.assertEqual(len(contig_rows), 2)
            self.assertEqual(contig_rows[0]["CRISPRS"], "0")
            self.assertEqual(contig_rows[0]["CRISPR_FRAC"], "0")
            self.assertEqual(contig_rows[1]["CRISPRS"], "0")
            self.assertEqual(contig_rows[1]["CRISPR_FRAC"], "0")
            self.assertEqual(crispr_rows, [])

    def test_main_accepts_json_with_null_optional_fields(self) -> None:
        """Ignore upstream optional fields that were normalised to JSON null."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            result_json = self.write_json(
                tmpdir / "result.json",
                {
                    "Sequences": [
                        {
                            "Id": "contig_nulls",
                            "Length": 1000,
                            "Crisprs": [
                                {
                                    "Name": "c_nulls",
                                    "Evidence_Level": 2,
                                    "Start": 10,
                                    "End": 100,
                                    "Spacers": 5,
                                    "Conservation_DRs": None,
                                    "Conservation_Spacers": 0,
                                }
                            ],
                        }
                    ]
                },
            )
            outdir = tmpdir / "out"

            exit_code = summarise_ccfinder.main(
                [
                    "--accession",
                    "ACC_NULLS",
                    "--result-json",
                    str(result_json),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            strain_row = read_tsv_rows(outdir / "ccfinder_strains.tsv")[0]
            crispr_row = read_tsv_rows(outdir / "ccfinder_crisprs.tsv")[0]
            self.assertEqual(strain_row["CRISPRS"], "1")
            self.assertEqual(strain_row["ccfinder_status"], "done")
            self.assertEqual(crispr_row["crispr_id"], "c_nulls")
            self.assertEqual(crispr_row["evidence_level"], "2")

    def test_build_strain_row_emits_stable_master_table_values(self) -> None:
        """Build stable strain-level master-table fields from retained CRISPR records."""
        crispr_records = [
            summarise_ccfinder.CrisprRecord(
                accession="ACC_STABLE",
                contig_id="contig1",
                crispr_id="c1",
                evidence_level=4,
                spacer_count=3,
                start=1,
                end=101,
            ),
            summarise_ccfinder.CrisprRecord(
                accession="ACC_STABLE",
                contig_id="contig2",
                crispr_id="c2",
                evidence_level=3,
                spacer_count=4,
                start=50,
                end=99,
            ),
        ]
        contig_summaries = [
            summarise_ccfinder.ContigSummary(
                accession="ACC_STABLE",
                contig_id="contig1",
                contig_length=1000,
                crisprs=1,
                spacers_sum=3,
                crispr_frac=0.101,
            ),
            summarise_ccfinder.ContigSummary(
                accession="ACC_STABLE",
                contig_id="contig2",
                contig_length=500,
                crisprs=1,
                spacers_sum=4,
                crispr_frac=0.1,
            ),
        ]

        row = summarise_ccfinder.build_strain_row(
            accession="ACC_STABLE",
            crispr_records=crispr_records,
            contig_summaries=contig_summaries,
            status="done",
            warnings=[],
        )

        self.assertEqual(
            row,
            {
                "accession": "ACC_STABLE",
                "CRISPRS": "2",
                "SPACERS_SUM": "7",
                "CRISPR_FRAC": "0.100667",
                "ccfinder_status": "done",
                "warnings": "",
            },
        )

    def test_main_marks_invalid_json_as_failed(self) -> None:
        """Emit NA-like strain output and empty detail tables on failure."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            outdir = tmpdir / "out"
            result_json = tmpdir / "missing.json"

            exit_code = summarise_ccfinder.main(
                [
                    "--accession",
                    "ACC3",
                    "--result-json",
                    str(result_json),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            strain_row = read_tsv_rows(outdir / "ccfinder_strains.tsv")[0]
            self.assertEqual(strain_row["CRISPRS"], "NA")
            self.assertEqual(strain_row["ccfinder_status"], "failed")
            self.assertEqual(strain_row["warnings"], "ccfinder_summary_failed")
            self.assertEqual(read_tsv_rows(outdir / "ccfinder_contigs.tsv"), [])
            self.assertEqual(read_tsv_rows(outdir / "ccfinder_crisprs.tsv"), [])

    def test_main_marks_malformed_json_as_failed(self) -> None:
        """Treat malformed JSON as a failed CRISPR summary."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            outdir = tmpdir / "out"
            result_json = tmpdir / "result.json"
            result_json.write_text("{\n", encoding="utf-8")

            exit_code = summarise_ccfinder.main(
                [
                    "--accession",
                    "ACC_BAD_JSON",
                    "--result-json",
                    str(result_json),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            strain_row = read_tsv_rows(outdir / "ccfinder_strains.tsv")[0]
            self.assertEqual(strain_row["CRISPRS"], "NA")
            self.assertEqual(strain_row["SPACERS_SUM"], "NA")
            self.assertEqual(strain_row["CRISPR_FRAC"], "NA")
            self.assertEqual(strain_row["ccfinder_status"], "failed")
            self.assertEqual(strain_row["warnings"], "ccfinder_summary_failed")
            self.assertEqual(read_tsv_rows(outdir / "ccfinder_contigs.tsv"), [])
            self.assertEqual(read_tsv_rows(outdir / "ccfinder_crisprs.tsv"), [])

    def test_main_generates_fallback_ids_and_invalid_entry_warning(self) -> None:
        """Handle missing IDs and skip malformed retained CRISPR entries with a warning."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            result_json = self.write_json(
                tmpdir / "result.json",
                {
                    "Sequences": [
                        {
                            "Length": 400,
                            "Crisprs": [
                                {
                                    "Evidence_Level": 4,
                                    "Start": 10,
                                    "End": 50,
                                    "Spacers": 2,
                                },
                                {
                                    "Evidence_Level": 4,
                                    "Start": 80,
                                    "End": 60,
                                    "Spacers": 1,
                                },
                            ],
                        }
                    ]
                },
            )
            outdir = tmpdir / "out"

            summarise_ccfinder.main(
                [
                    "--accession",
                    "ACC4",
                    "--result-json",
                    str(result_json),
                    "--outdir",
                    str(outdir),
                ]
            )

            strain_row = read_tsv_rows(outdir / "ccfinder_strains.tsv")[0]
            contig_row = read_tsv_rows(outdir / "ccfinder_contigs.tsv")[0]
            crispr_row = read_tsv_rows(outdir / "ccfinder_crisprs.tsv")[0]

            self.assertEqual(strain_row["CRISPRS"], "1")
            self.assertEqual(strain_row["warnings"], "invalid_crispr_entry")
            self.assertEqual(contig_row["contig_id"], "contig_1")
            self.assertEqual(crispr_row["crispr_id"], "contig_1_crispr_1")


if __name__ == "__main__":
    unittest.main()
