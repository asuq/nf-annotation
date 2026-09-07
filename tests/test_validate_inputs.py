"""Tests for the input-validation CLI."""

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

import master_table_contract  # noqa: E402
import validate_inputs  # noqa: E402


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a TSV file into dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class ValidateInputsTestCase(unittest.TestCase):
    """Cover the strict sample-manifest and metadata validation contract."""

    def test_empty_manifest_is_rejected(self) -> None:
        with self.assertRaisesRegex(validate_inputs.ValidationError, "at least one"):
            validate_inputs.validate_samples(
                validate_inputs.REQUIRED_SAMPLE_COLUMNS, [], {}
            )

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write a UTF-8 text file and return its path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def test_main_writes_expected_outputs_and_collision_safe_ids(self) -> None:
        """Validate manifest success, metadata fallback warnings, and collision suffixes."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            fasta_one = self.write_text_file(tmpdir / "genomes" / "a.fna", ">a\nACGT\n")
            fasta_two = self.write_text_file(tmpdir / "genomes" / "b.fna", ">b\nTGCA\n")
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "\n".join(
                    [
                        "accession,is_new,assembly_level,genome_fasta,organism_name",
                        f"ACC-1,false,,{fasta_one},Known one",
                        f"ACC 1,true,Scaffold,{fasta_two},Candidate two",
                    ]
                )
                + "\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.csv",
                "\n".join(
                    [
                        "accession,Tax_ID,Organism_Name",
                        "ACC-1,123,Known one",
                    ]
                )
                + "\n",
            )
            outdir = tmpdir / "validated"

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            validated_rows = read_tsv(outdir / "validated_samples.tsv")
            accession_map_rows = read_tsv(outdir / "accession_map.tsv")
            warning_rows = read_tsv(outdir / "validation_warnings.tsv")
            sample_status_rows = read_tsv(outdir / "sample_status.tsv")

            self.assertEqual(len(validated_rows), 2)
            self.assertEqual(validated_rows[0]["is_new"], "false")
            self.assertEqual(validated_rows[1]["is_new"], "true")
            self.assertEqual(validated_rows[0]["assembly_level"], "NA")
            self.assertEqual(validated_rows[1]["assembly_level"], "Scaffold")
            self.assertEqual(validated_rows[0]["genome_fasta"], str(fasta_one.resolve()))
            self.assertEqual(validated_rows[1]["genome_fasta"], str(fasta_two.resolve()))

            internal_ids = {row["accession"]: row["internal_id"] for row in accession_map_rows}
            self.assertEqual(len(set(internal_ids.values())), 2)
            self.assertTrue(internal_ids["ACC-1"].startswith("ACC_1_"))
            self.assertTrue(internal_ids["ACC 1"].startswith("ACC_1_"))

            warning_codes = {
                (row["accession"], row["warning_code"]) for row in warning_rows
            }
            self.assertIn(
                ("ACC 1", "missing_metadata_for_new_sample"),
                warning_codes,
            )
            self.assertIn(
                ("ACC-1", "internal_id_collision_resolved"),
                warning_codes,
            )
            self.assertIn(
                ("ACC 1", "internal_id_collision_resolved"),
                warning_codes,
            )
            self.assertEqual(len(sample_status_rows), 2)
            sample_status_by_accession = {
                row["accession"]: row for row in sample_status_rows
            }
            self.assertEqual(
                sample_status_by_accession["ACC-1"]["validation_status"],
                "done",
            )
            self.assertEqual(
                sample_status_by_accession["ACC-1"]["internal_id"],
                internal_ids["ACC-1"],
            )
            self.assertEqual(
                sample_status_by_accession["ACC-1"]["warnings"],
                "internal_id_collision_resolved",
            )
            self.assertEqual(
                sample_status_by_accession["ACC 1"]["warnings"],
                "internal_id_collision_resolved;missing_metadata_for_new_sample",
            )
            self.assertEqual(
                sample_status_by_accession["ACC 1"]["taxonomy_status"],
                "na",
            )
            self.assertEqual(
                sample_status_by_accession["ACC 1"]["gcode"],
                "NA",
            )
            self.assertEqual(
                sample_status_by_accession["ACC 1"]["ani_included"],
                "na",
            )

    def test_main_canonicalises_sample_assembly_levels(self) -> None:
        """Canonicalise sample-manifest assembly levels before staging outputs."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            fasta_one = self.write_text_file(tmpdir / "genomes" / "one.fna", ">a\nACGT\n")
            fasta_two = self.write_text_file(tmpdir / "genomes" / "two.fna", ">b\nTGCA\n")
            fasta_three = self.write_text_file(tmpdir / "genomes" / "three.fna", ">c\nGATC\n")
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "\n".join(
                    [
                        "accession,is_new,assembly_level,genome_fasta",
                        f"ACC1,true,Scaffolds,{fasta_one}",
                        f"ACC2,true,cOnTiGs,{fasta_two}",
                        f"ACC3,false,CHROMOSOMES,{fasta_three}",
                    ]
                )
                + "\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.csv",
                "\n".join(
                    [
                        "accession,Tax_ID",
                        "ACC3,3",
                    ]
                )
                + "\n",
            )
            outdir = tmpdir / "validated"

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            validated_rows = {
                row["accession"]: row for row in read_tsv(outdir / "validated_samples.tsv")
            }
            accession_map_rows = {
                row["accession"]: row for row in read_tsv(outdir / "accession_map.tsv")
            }
            self.assertEqual(validated_rows["ACC1"]["assembly_level"], "Scaffold")
            self.assertEqual(validated_rows["ACC2"]["assembly_level"], "Contig")
            self.assertEqual(validated_rows["ACC3"]["assembly_level"], "Chromosome")
            self.assertEqual(accession_map_rows["ACC1"]["assembly_level"], "Scaffold")
            self.assertEqual(accession_map_rows["ACC2"]["assembly_level"], "Contig")
            self.assertEqual(accession_map_rows["ACC3"]["assembly_level"], "Chromosome")

    def test_main_uses_staged_sample_status_columns_asset(self) -> None:
        """Allow one valid sample-status contract file to be provided explicitly."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            fasta_path = self.write_text_file(tmpdir / "genomes" / "a.fna", ">a\nACGT\n")
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "\n".join(
                    [
                        "accession,is_new,assembly_level,genome_fasta",
                        f"ACC1,false,,{fasta_path}",
                    ]
                )
                + "\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.csv",
                "\n".join(
                    [
                        "accession,Tax_ID,Organism_Name",
                        "ACC1,123,Known one",
                    ]
                )
                + "\n",
            )
            sample_status_columns = self.write_text_file(
                tmpdir / "sample_status_columns.txt",
                "\n".join(
                    master_table_contract.build_sample_status_columns(["custom_odb12"])
                )
                + "\n",
            )
            outdir = tmpdir / "validated"

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--sample-status-columns",
                    str(sample_status_columns),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            sample_status_rows = read_tsv(outdir / "sample_status.tsv")
            self.assertIn("busco_custom_odb12_status", sample_status_rows[0])
            self.assertEqual(sample_status_rows[0]["validation_status"], "done")

    def test_main_seeds_custom_busco_status_columns_from_lineages(self) -> None:
        """Seed BUSCO status columns directly from configured lineage names."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            fasta_path = self.write_text_file(tmpdir / "genomes" / "a.fna", ">a\nACGT\n")
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "accession,is_new,assembly_level,genome_fasta\n"
                f"ACC1,false,,{fasta_path}\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.csv",
                "accession,Tax_ID,Organism_Name\nACC1,123,Known one\n",
            )
            outdir = tmpdir / "validated"

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--busco-lineage",
                    "custom_odb12",
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            sample_status_rows = read_tsv(outdir / "sample_status.tsv")
            self.assertEqual(sample_status_rows[0]["busco_custom_odb12_status"], "na")
            self.assertEqual(sample_status_rows[0]["gcode"], "NA")
            self.assertEqual(sample_status_rows[0]["ani_included"], "na")

    def test_main_rejects_missing_required_sample_headers(self) -> None:
        """Fail when the sample manifest omits a required locked header."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            fasta_path = self.write_text_file(tmpdir / "genome.fna", ">a\nACGT\n")
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "\n".join(
                    [
                        "accession,is_new,genome_fasta",
                        f"ACC1,false,{fasta_path}",
                    ]
                )
                + "\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.csv",
                "accession,Tax_ID\nACC1,1\n",
            )

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--outdir",
                    str(tmpdir / "validated"),
                ]
            )

            self.assertEqual(exit_code, 1)

    def test_main_rejects_existing_samples_missing_metadata(self) -> None:
        """Fail when an existing sample has no matching metadata row."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            fasta_path = self.write_text_file(tmpdir / "genome.fna", ">a\nACGT\n")
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "\n".join(
                    [
                        "accession,is_new,assembly_level,genome_fasta",
                        f"ACC1,false,,{fasta_path}",
                    ]
                )
                + "\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\nOTHER\t1\n",
            )
            outdir = tmpdir / "validated"

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse((outdir / "validated_samples.tsv").exists())

    def test_main_rejects_accessions_with_path_separators(self) -> None:
        """Fail when an accession is unsafe for published output folders."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            fasta_path = self.write_text_file(tmpdir / "genome.fna", ">a\nACGT\n")
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "\n".join(
                    [
                        "accession,is_new,assembly_level,genome_fasta",
                        f"BAD/ACC,false,,{fasta_path}",
                    ]
                )
                + "\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.csv",
                "accession,Tax_ID\nBAD/ACC,1\n",
            )
            outdir = tmpdir / "validated"

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse((outdir / "accession_map.tsv").exists())

    def test_main_rejects_duplicate_sample_accessions(self) -> None:
        """Fail when the sample manifest repeats the same accession."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            fasta_one = self.write_text_file(tmpdir / "a.fna", ">a\nACGT\n")
            fasta_two = self.write_text_file(tmpdir / "b.fna", ">b\nTGCA\n")
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "\n".join(
                    [
                        "accession,is_new,assembly_level,genome_fasta",
                        f"ACC1,false,,{fasta_one}",
                        f"ACC1,true,Scaffold,{fasta_two}",
                    ]
                )
                + "\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.csv",
                "accession,Tax_ID\nACC1,1\n",
            )

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--outdir",
                    str(tmpdir / "validated"),
                ]
            )

            self.assertEqual(exit_code, 1)

    def test_main_rejects_invalid_boolean_values(self) -> None:
        """Fail when `is_new` cannot be normalised to true or false."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            fasta_path = self.write_text_file(tmpdir / "genome.fna", ">a\nACGT\n")
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "\n".join(
                    [
                        "accession,is_new,assembly_level,genome_fasta",
                        f"ACC1,maybe,,{fasta_path}",
                    ]
                )
                + "\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.csv",
                "accession,Tax_ID\nACC1,1\n",
            )

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--outdir",
                    str(tmpdir / "validated"),
                ]
            )

            self.assertEqual(exit_code, 1)

    def test_main_allows_existing_samples_without_assembly_level(self) -> None:
        """Allow blank assembly levels when `is_new=false` and metadata is present."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            fasta_path = self.write_text_file(tmpdir / "genome.fna", ">a\nACGT\n")
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "\n".join(
                    [
                        "accession,is_new,assembly_level,genome_fasta",
                        f"ACC1,false,,{fasta_path}",
                    ]
                )
                + "\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.csv",
                "accession,Tax_ID\nACC1,1\n",
            )
            outdir = tmpdir / "validated"

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            validated_rows = read_tsv(outdir / "validated_samples.tsv")
            self.assertEqual(validated_rows[0]["assembly_level"], "NA")

    def test_main_rejects_uppercase_sample_headers(self) -> None:
        """Fail when required sample-manifest headers are not lower-case."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            fasta_path = self.write_text_file(tmpdir / "genome.fna", ">a\nACGT\n")
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "\n".join(
                    [
                        "accession,is_new,Assembly_Level,genome_fasta",
                        f"ACC1,false,,{fasta_path}",
                    ]
                )
                + "\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.csv",
                "accession,Tax_ID\nACC1,1\n",
            )

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--outdir",
                    str(tmpdir / "validated"),
                ]
            )

            self.assertEqual(exit_code, 1)

    def test_main_rejects_new_samples_without_assembly_level(self) -> None:
        """Fail when `is_new=true` rows omit the required assembly_level."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            fasta_path = self.write_text_file(tmpdir / "genome.fna", ">a\nACGT\n")
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "\n".join(
                    [
                        "accession,is_new,assembly_level,genome_fasta",
                        f"ACC1,true,,{fasta_path}",
                    ]
                )
                + "\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.csv",
                "accession,Tax_ID\nACC1,1\n",
            )

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--outdir",
                    str(tmpdir / "validated"),
                ]
            )

            self.assertEqual(exit_code, 1)

    def test_main_rejects_invalid_new_sample_assembly_level(self) -> None:
        """Fail when a new sample uses a non-parseable assembly level."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            fasta_path = self.write_text_file(tmpdir / "genome.fna", ">a\nACGT\n")
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "\n".join(
                    [
                        "accession,is_new,assembly_level,genome_fasta",
                        f"ACC1,true,metagenome,{fasta_path}",
                    ]
                )
                + "\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.csv",
                "accession,Tax_ID\nACC1,1\n",
            )
            outdir = tmpdir / "validated"

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse((outdir / "validated_samples.tsv").exists())

    def test_main_rejects_missing_genome_fasta_paths(self) -> None:
        """Fail when a manifest row points to a genome FASTA that does not exist."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            missing_fasta_path = tmpdir / "missing.fna"
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "\n".join(
                    [
                        "accession,is_new,assembly_level,genome_fasta",
                        f"ACC1,false,,{missing_fasta_path}",
                    ]
                )
                + "\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.csv",
                "accession,Tax_ID\nACC1,1\n",
            )

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--outdir",
                    str(tmpdir / "validated"),
                ]
            )

            self.assertEqual(exit_code, 1)

    def test_main_can_defer_missing_genome_fasta_checks(self) -> None:
        """Allow workflow-level staging to own host-path existence checks."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            missing_fasta_path = tmpdir / "missing.fna"
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "\n".join(
                    [
                        "accession,is_new,assembly_level,genome_fasta",
                        f"ACC1,false,,{missing_fasta_path}",
                    ]
                )
                + "\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.csv",
                "accession,Tax_ID\nACC1,1\n",
            )
            outdir = tmpdir / "validated"

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--defer-genome-fasta-check",
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            validated_rows = read_tsv(outdir / "validated_samples.tsv")
            self.assertEqual(validated_rows[0]["genome_fasta"], str(missing_fasta_path.resolve()))

    def test_main_rejects_accessions_that_sanitise_to_empty_ids(self) -> None:
        """Fail when sanitisation removes every character from the accession."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            fasta_path = self.write_text_file(tmpdir / "genome.fna", ">a\nACGT\n")
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "\n".join(
                    [
                        "accession,is_new,assembly_level,genome_fasta",
                        f"!!!,false,,{fasta_path}",
                    ]
                )
                + "\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.csv",
                "accession,Tax_ID\n!!!,1\n",
            )

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--outdir",
                    str(tmpdir / "validated"),
                ]
            )

            self.assertEqual(exit_code, 1)

    def test_main_rejects_duplicate_metadata_accessions(self) -> None:
        """Fail when the metadata table contains ambiguous duplicate keys."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            fasta_path = self.write_text_file(tmpdir / "genome.fna", ">a\nACGT\n")
            sample_csv = self.write_text_file(
                tmpdir / "samples.csv",
                "\n".join(
                    [
                        "accession,is_new,assembly_level,genome_fasta",
                        f"ACC1,false,,{fasta_path}",
                    ]
                )
                + "\n",
            )
            metadata_csv = self.write_text_file(
                tmpdir / "metadata.csv",
                "\n".join(
                    [
                        "accession,Tax_ID",
                        "ACC1,1",
                        "ACC1,2",
                    ]
                )
                + "\n",
            )

            exit_code = validate_inputs.main(
                [
                    "--sample-csv",
                    str(sample_csv),
                    "--metadata",
                    str(metadata_csv),
                    "--outdir",
                    str(tmpdir / "validated"),
                ]
            )

            self.assertEqual(exit_code, 1)

    def test_sanitise_accession_replaces_non_ascii_characters(self) -> None:
        """Keep only ASCII alphanumerics and underscores in internal IDs."""
        self.assertEqual(validate_inputs.sanitise_accession("A-B C"), "A_B_C")
        self.assertEqual(validate_inputs.sanitise_accession("A__B"), "A_B")
        self.assertEqual(validate_inputs.sanitise_accession("__A--B__"), "A_B")
        self.assertEqual(validate_inputs.sanitise_accession("Acc\u00e9ssion"), "Acc_ssion")

    def test_add_collision_suffixes_is_deterministic(self) -> None:
        """Assign the same collision suffixes regardless of input ordering."""
        forward = validate_inputs.add_collision_suffixes(["ACC-1", "ACC 1"])
        reverse = validate_inputs.add_collision_suffixes(["ACC 1", "ACC-1"])

        self.assertEqual(forward, reverse)
        self.assertEqual(set(forward), {"ACC-1", "ACC 1"})
        self.assertTrue(forward["ACC-1"].startswith("ACC_1_"))
        self.assertTrue(forward["ACC 1"].startswith("ACC_1_"))
        self.assertNotEqual(forward["ACC-1"], forward["ACC 1"])


if __name__ == "__main__":
    unittest.main()
