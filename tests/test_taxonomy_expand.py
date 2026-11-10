"""Tests for taxonomy expansion from a pinned NCBI taxdump."""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import taxonomy_expand


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a TSV file into row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class TaxonomyExpandTestCase(unittest.TestCase):
    """Cover lineage expansion and invalid Tax_ID handling."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write text to a test file and return the path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def test_expand_requested_tax_ids_from_metadata(self) -> None:
        """Expand valid metadata Tax_ID values into lineage columns."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                (
                    "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id\n"
                    "ACC1\tfalse\tNA\t/tmp/acc1.fna\tACC1\n"
                    "ACC2\tfalse\tNA\t/tmp/acc2.fna\tACC2"
                )
                + "\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\nACC1\t8\nACC2\t999999" + "\n",
            )
            taxdump = tmpdir / "taxdump"
            self.write_text_file(
                taxdump / "nodes.dmp",
                (
                    "1\t|\t1\t|\tno rank\t|\n"
                    "2\t|\t1\t|\tsuperkingdom\t|\n"
                    "3\t|\t2\t|\tphylum\t|\n"
                    "4\t|\t3\t|\tclass\t|\n"
                    "5\t|\t4\t|\torder\t|\n"
                    "6\t|\t5\t|\tfamily\t|\n"
                    "7\t|\t6\t|\tgenus\t|\n"
                    "8\t|\t7\t|\tspecies\t|"
                )
                + "\n",
            )
            self.write_text_file(
                taxdump / "names.dmp",
                (
                    "1\t|\troot\t|\t\t|\tscientific name\t|\n"
                    "2\t|\tBacteria\t|\t\t|\tscientific name\t|\n"
                    "3\t|\tFirmicutes\t|\t\t|\tscientific name\t|\n"
                    "4\t|\tMollicutes\t|\t\t|\tscientific name\t|\n"
                    "5\t|\tMycoplasmatales\t|\t\t|\tscientific name\t|\n"
                    "6\t|\tMycoplasmataceae\t|\t\t|\tscientific name\t|\n"
                    "7\t|\tMycoplasma\t|\t\t|\tscientific name\t|\n"
                    "8\t|\tMycoplasma testus\t|\t\t|\tscientific name\t|"
                )
                + "\n",
            )
            output = tmpdir / "taxonomy_expanded.tsv"

            exit_code = taxonomy_expand.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--taxdump",
                    str(taxdump),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            rows = read_tsv(output)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["Tax_ID"], "8")
            self.assertEqual(rows[0]["superkingdom"], "Bacteria")
            self.assertEqual(rows[0]["species"], "Mycoplasma testus")

    def test_missing_tax_id_column_writes_header_only(self) -> None:
        """Handle metadata tables without Tax_ID by emitting an empty output table."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id\nACC1\tfalse\tNA\t/tmp/acc1.fna\tACC1\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tOrganism_Name\nACC1\tExample\n",
            )
            taxdump = tmpdir / "taxdump"
            self.write_text_file(taxdump / "nodes.dmp", "1\t|\t1\t|\tno rank\t|\n")
            self.write_text_file(
                taxdump / "names.dmp",
                "1\t|\troot\t|\t\t|\tscientific name\t|\n",
            )
            output = tmpdir / "taxonomy_expanded.tsv"

            exit_code = taxonomy_expand.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--taxdump",
                    str(taxdump),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertEqual(read_tsv(output), [])


if __name__ == "__main__":
    unittest.main()
