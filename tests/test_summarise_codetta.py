"""Tests for Codetta output summarisation and NCBI table matching."""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

from Bio.Data import CodonTable

ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import summarise_codetta


def read_tsv_row(path: Path) -> dict[str, str]:
    """Read one single-row TSV output."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    return rows[0]


class SummariseCodettaTestCase(unittest.TestCase):
    """Cover Codetta parsing, fallback handling, and table matching."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write one UTF-8 text file."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def write_inference_table(self, directory: Path, genetic_code: str) -> None:
        """Write one Codetta-style inference table from a 64-character code."""
        lines = [
            "# Codon inferences",
            "# codon   inference  std code  diff?    N aligned  N used",
        ]
        for codon, amino_acid in zip(
            summarise_codetta.CODON_ORDER, genetic_code, strict=True
        ):
            lines.append(f"{codon}\t{amino_acid}\tX\t.\t1\t1")
        lines.append("#")
        self.write_text_file(
            directory / "codetta_inference.txt", "\n".join(lines) + "\n"
        )

    def run_cli(self, *, accession: str, codetta_dir: Path) -> dict[str, str]:
        """Run the Codetta summariser and return the emitted row."""
        output = codetta_dir / "codetta_summary.tsv"
        exit_code = summarise_codetta.main(
            [
                "--accession",
                accession,
                "--codetta-dir",
                str(codetta_dir),
                "--output",
                str(output),
            ]
        )
        self.assertEqual(exit_code, 0)
        return read_tsv_row(output)

    def test_main_reconstructs_table_11_equivalent_candidates_from_inference_rows(
        self,
    ) -> None:
        """Return every compatible table when table 11 is not uniquely identifiable."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            codetta_dir = Path(tmpdir_name)
            genetic_code = summarise_codetta.build_ncbi_table_string(
                CodonTable.generic_by_id[11]
            )
            self.write_inference_table(codetta_dir, genetic_code)

            row = self.run_cli(accession="ACC11", codetta_dir=codetta_dir)

            self.assertEqual(row["accession"], "ACC11")
            self.assertEqual(row["Codetta_Genetic_Code"], genetic_code)
            self.assertEqual(row["Codetta_NCBI_Table_Candidates"], "1;11")
            self.assertEqual(row["codetta_status"], "done")
            self.assertEqual(row["warnings"], "codetta_multiple_ncbi_tables")

    def test_main_reconstructs_table_4_from_inference_rows(self) -> None:
        """Parse the 64-codon table and recover a unique NCBI table 4 match."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            codetta_dir = Path(tmpdir_name)
            genetic_code = summarise_codetta.build_ncbi_table_string(
                CodonTable.generic_by_id[4]
            )
            self.write_inference_table(codetta_dir, genetic_code)

            row = self.run_cli(accession="ACC4", codetta_dir=codetta_dir)

            self.assertEqual(row["Codetta_NCBI_Table_Candidates"], "4")
            self.assertEqual(row["warnings"], "")

    def test_main_reconstructs_one_nonstandard_table_from_inference_rows(self) -> None:
        """Keep support for one non-4-or-11 NCBI candidate match."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            codetta_dir = Path(tmpdir_name)
            genetic_code = summarise_codetta.build_ncbi_table_string(
                CodonTable.generic_by_id[2]
            )
            self.write_inference_table(codetta_dir, genetic_code)

            row = self.run_cli(accession="ACC2", codetta_dir=codetta_dir)

            self.assertEqual(row["Codetta_NCBI_Table_Candidates"], "2")
            self.assertEqual(row["warnings"], "")

    def test_main_falls_back_to_the_log_line_when_needed(self) -> None:
        """Use the terminal genetic-code line when the inference table is unusable."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            codetta_dir = Path(tmpdir_name)
            genetic_code = summarise_codetta.build_ncbi_table_string(
                CodonTable.generic_by_id[11]
            )
            self.write_text_file(codetta_dir / "codetta_inference.txt", "# broken\n")
            self.write_text_file(
                codetta_dir / "codetta.log", f"Genetic code: {genetic_code}\n"
            )

            row = self.run_cli(accession="ACC_LOG", codetta_dir=codetta_dir)

            self.assertEqual(row["Codetta_Genetic_Code"], genetic_code)
            self.assertEqual(row["Codetta_NCBI_Table_Candidates"], "1;11")
            self.assertEqual(row["warnings"], "codetta_multiple_ncbi_tables")

    def test_main_reports_multiple_matching_candidate_tables(self) -> None:
        """Return all compatible NCBI tables in stable ascending order."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            codetta_dir = Path(tmpdir_name)
            table_4 = summarise_codetta.build_ncbi_table_string(
                CodonTable.generic_by_id[4]
            )
            table_11 = summarise_codetta.build_ncbi_table_string(
                CodonTable.generic_by_id[11]
            )
            genetic_code = "".join(
                observed if observed == expected else "?"
                for observed, expected in zip(table_4, table_11, strict=True)
            )
            self.write_inference_table(codetta_dir, genetic_code)

            expected_candidates = summarise_codetta.match_ncbi_tables(genetic_code)
            self.assertGreater(len(expected_candidates), 1)

            row = self.run_cli(accession="ACC_MULTI", codetta_dir=codetta_dir)

            self.assertEqual(
                row["Codetta_NCBI_Table_Candidates"],
                ";".join(str(table_id) for table_id in expected_candidates),
            )
            self.assertEqual(row["warnings"], "codetta_multiple_ncbi_tables")

    def test_main_reports_unassigned_when_no_table_matches(self) -> None:
        """Return unassigned when the observed code fits no NCBI table."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            codetta_dir = Path(tmpdir_name)
            table_11 = summarise_codetta.build_ncbi_table_string(
                CodonTable.generic_by_id[11]
            )
            genetic_code = "!" + table_11[1:]
            self.write_inference_table(codetta_dir, genetic_code)

            row = self.run_cli(accession="ACC_NONE", codetta_dir=codetta_dir)

            self.assertEqual(row["Codetta_NCBI_Table_Candidates"], "unassigned")
            self.assertEqual(row["warnings"], "codetta_unassigned")

    def test_main_writes_failure_row_when_outputs_are_missing(self) -> None:
        """Return a stable failure row when Codetta did not emit usable files."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            codetta_dir = Path(tmpdir_name)

            row = self.run_cli(accession="ACC_FAIL", codetta_dir=codetta_dir)

            self.assertEqual(row["Codetta_Genetic_Code"], "NA")
            self.assertEqual(row["Codetta_NCBI_Table_Candidates"], "NA")
            self.assertEqual(row["codetta_status"], "failed")
            self.assertEqual(row["warnings"], "codetta_failed")


if __name__ == "__main__":
    unittest.main()
