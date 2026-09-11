"""Tests for the HPC matrix validation helper."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "validate_hpc_matrix.py"
sys.path.insert(0, str(ROOT / "bin"))
SPEC = importlib.util.spec_from_file_location("validate_hpc_matrix", SCRIPT)
assert SPEC is not None
assert SPEC.loader is not None
validate_hpc_matrix = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(validate_hpc_matrix)


def write_text(path: Path, text: str) -> None:
    """Write one text file, creating parent directories first."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


class ValidateHpcMatrixTestCase(unittest.TestCase):
    """Exercise the HPC matrix validation helper."""

    def run_helper(self, *args: str) -> subprocess.CompletedProcess[str]:
        """Run the helper with text capture enabled."""
        return subprocess.run(
            ["python3", str(SCRIPT), *args],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_dbprep_validation_accepts_minimal_prepared_tree(self) -> None:
        """Validate a minimal successful dbprep case."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            results_dir = root / "results"
            write_text(
                results_dir / "runtime_database_report.tsv",
                (
                    "component\tstatus\tsource\tdestination\tdetails\n"
                    f"taxdump\tprepared\tnative_download\t{root / 'db' / 'taxdump'}\tfiles=2\n"
                ),
            )
            write_text(
                results_dir / "nextflow_args.txt",
                f"--taxdump {root / 'db' / 'taxdump'}\n",
            )
            for name in ("trace.tsv", "report.html", "timeline.html", "dag.html"):
                write_text(results_dir / "pipeline_info" / name, "ok\n")

            taxdump_dir = root / "db" / "taxdump"
            write_text(taxdump_dir / "names.dmp", "names\n")
            write_text(taxdump_dir / "nodes.dmp", "nodes\n")
            write_text(taxdump_dir / ".nf_myco_ready.json", "{}\n")

            result = self.run_helper(
                "dbprep",
                "--results-dir",
                str(results_dir),
                "--taxdump",
                str(taxdump_dir),
                "--expected-component",
                "taxdump",
                "--expected-status",
                "prepared",
                "--expected-arg=--taxdump",
            )

            self.assertEqual(result.returncode, 0, msg=result.stderr)

    def test_dbprep_validation_accepts_codetta_component(self) -> None:
        """Validate one minimal successful Codetta dbprep case."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            results_dir = root / "results"
            codetta_dir = root / "db" / "codetta"
            write_text(
                results_dir / "runtime_database_report.tsv",
                (
                    "component\tstatus\tsource\tdestination\tdetails\n"
                    f"codetta\tprepared\tnative_download\t{codetta_dir}\tfiles=5\n"
                ),
            )
            write_text(
                results_dir / "nextflow_args.txt",
                f"--codetta_db {codetta_dir}\n",
            )
            for name in ("trace.tsv", "report.html", "timeline.html", "dag.html"):
                write_text(results_dir / "pipeline_info" / name, "ok\n")

            for suffix in ("", ".h3f", ".h3i", ".h3m", ".h3p"):
                write_text(codetta_dir / f"Pfam-A_enone.hmm{suffix}", "stub\n")
            write_text(codetta_dir / ".nf_myco_ready.json", "{}\n")

            result = self.run_helper(
                "dbprep",
                "--results-dir",
                str(results_dir),
                "--codetta-db",
                str(codetta_dir),
                "--expected-component",
                "codetta",
                "--expected-status",
                "prepared",
                "--expected-arg=--codetta_db",
            )

            self.assertEqual(result.returncode, 0, msg=result.stderr)

    def test_dbprep_validation_rejects_codetta_without_ready_marker(self) -> None:
        """Reject one Codetta dbprep result when the destination lacks the marker."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            results_dir = root / "results"
            codetta_dir = root / "db" / "codetta"
            write_text(
                results_dir / "runtime_database_report.tsv",
                (
                    "component\tstatus\tsource\tdestination\tdetails\n"
                    f"codetta\tprepared\tnative_download\t{codetta_dir}\tfiles=5\n"
                ),
            )
            write_text(
                results_dir / "nextflow_args.txt",
                f"--codetta_db {codetta_dir}\n",
            )
            for name in ("trace.tsv", "report.html", "timeline.html", "dag.html"):
                write_text(results_dir / "pipeline_info" / name, "ok\n")

            for suffix in ("", ".h3f", ".h3i", ".h3m", ".h3p"):
                write_text(codetta_dir / f"Pfam-A_enone.hmm{suffix}", "stub\n")

            result = self.run_helper(
                "dbprep",
                "--results-dir",
                str(results_dir),
                "--codetta-db",
                str(codetta_dir),
                "--expected-component",
                "codetta",
                "--expected-status",
                "prepared",
                "--expected-arg=--codetta_db",
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Missing ready marker for codetta", result.stderr)

    def test_medium_run_validation_accepts_allowed_phyla(self) -> None:
        """Validate a minimal successful medium real-data run."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            outdir = root / "out"
            write_text(
                outdir / "tables" / "master_table.tsv",
                (
                    "Accession\tphylum\n"
                    "ACC1\tMycoplasmatota\n"
                    "ACC2\tBacillota\n"
                ),
            )
            write_text(
                outdir / "tables" / "sample_status.tsv",
                "accession\tvalidation_status\nACC1\tdone\nACC2\tdone\n",
            )
            write_text(
                outdir / "tables" / "tool_and_db_versions.tsv",
                (
                    "component\tkind\tversion\timage_or_path\tnotes\n"
                    f"eggnog_mapper\ttool\t2.1.13\tNA\tNA\n"
                    f"taxdump\tdatabase\t20240914\t{(root / 'db' / 'taxdump').resolve()}\tNA\n"
                ),
            )
            for name in ("trace.tsv", "report.html", "timeline.html", "dag.html"):
                write_text(outdir / "pipeline_info" / name, "ok\n")

            result = self.run_helper(
                "medium-run",
                "--outdir",
                str(outdir),
                "--db-root",
                str(root / "db"),
                "--sample-count",
                "2",
                "--allowed-phylum",
                "Mycoplasmatota",
                "--allowed-phylum",
                "Bacillota",
            )

            self.assertEqual(result.returncode, 0, msg=result.stderr)

    def test_medium_run_validation_rejects_unresolved_gcode(self) -> None:
        """Unresolved genetic code is a failed comparison under the ratio rule."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            outdir = root / "out"
            write_text(
                outdir / "tables" / "master_table.tsv",
                (
                    "Accession\tphylum\n"
                    "ACC1\tMycoplasmatota\n"
                    "ACC2\tBacillota\n"
                ),
            )
            write_text(
                outdir / "tables" / "sample_status.tsv",
                (
                    "accession\tvalidation_status\tgcode_status\tcodetta_status\t"
                    "prokka_status\tccfinder_status\tpadloc_status\teggnog_status\twarnings\n"
                    "ACC1\tdone\tdone\tdone\tdone\tdone\tdone\tdone\tcodetta_multiple_ncbi_tables\n"
                    "ACC2\tdone\tfailed\tdone\tskipped\tskipped\tskipped\tskipped\t"
                    "gcode_na;codetta_multiple_ncbi_tables\n"
                ),
            )
            write_text(
                outdir / "tables" / "tool_and_db_versions.tsv",
                (
                    "component\tkind\tversion\timage_or_path\tnotes\n"
                    "eggnog_mapper\ttool\t2.1.13\tNA\tNA\n"
                    f"taxdump\tdatabase\t20240914\t{(root / 'db' / 'taxdump').resolve()}\tNA\n"
                ),
            )
            for name in ("trace.tsv", "report.html", "timeline.html", "dag.html"):
                write_text(outdir / "pipeline_info" / name, "ok\n")

            result = self.run_helper(
                "medium-run",
                "--outdir",
                str(outdir),
                "--db-root",
                str(root / "db"),
                "--sample-count",
                "2",
                "--allowed-phylum",
                "Mycoplasmatota",
                "--allowed-phylum",
                "Bacillota",
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("unexpected medium sample-status failures", result.stderr)

    def test_medium_run_validation_accepts_isolated_secondary_busco_failure(self) -> None:
        """Allow one medium run where only the non-primary BUSCO lineage failed."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            outdir = root / "out"
            write_text(
                outdir / "tables" / "master_table.tsv",
                (
                    "Accession\tphylum\n"
                    "ACC1\tMycoplasmatota\n"
                    "ACC2\tBacillota\n"
                ),
            )
            write_text(
                outdir / "tables" / "sample_status.tsv",
                (
                    "accession\tvalidation_status\tgcode_status\t"
                    "busco_bacillota_odb12_status\tbusco_mycoplasmatota_odb12_status\t"
                    "codetta_status\tprokka_status\tccfinder_status\tpadloc_status\t"
                    "eggnog_status\twarnings\n"
                    "ACC1\tdone\tdone\tdone\tdone\tdone\tdone\tdone\tdone\tdone\t\n"
                    "ACC2\tdone\tdone\tdone\tfailed\tdone\tdone\tdone\tdone\tdone\t"
                    "busco_summary_failed\n"
                ),
            )
            write_text(
                outdir / "tables" / "tool_and_db_versions.tsv",
                (
                    "component\tkind\tversion\timage_or_path\tnotes\n"
                    "eggnog_mapper\ttool\t2.1.13\tNA\tNA\n"
                    f"taxdump\tdatabase\t20240914\t{(root / 'db' / 'taxdump').resolve()}\tNA\n"
                ),
            )
            for name in ("trace.tsv", "report.html", "timeline.html", "dag.html"):
                write_text(outdir / "pipeline_info" / name, "ok\n")

            result = self.run_helper(
                "medium-run",
                "--outdir",
                str(outdir),
                "--db-root",
                str(root / "db"),
                "--sample-count",
                "2",
                "--allowed-phylum",
                "Mycoplasmatota",
                "--allowed-phylum",
                "Bacillota",
            )

            self.assertEqual(result.returncode, 0, msg=result.stderr)

    def test_versions_table_validation_accepts_symlinked_db_paths(self) -> None:
        """Allow reported DB paths that use one symlinked alias of the HPC root."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            real_db_root = root / "real_db"
            real_db_root.mkdir()
            alias_db_root = root / "flash_db"
            alias_db_root.symlink_to(real_db_root, target_is_directory=True)

            versions_path = root / "tool_and_db_versions.tsv"
            write_text(
                versions_path,
                (
                    "component\tkind\tversion\timage_or_path\tnotes\n"
                    "eggnog_mapper\ttool\t2.1.13\tNA\tNA\n"
                    f"taxdump\tdatabase\t20240914\t{alias_db_root / 'ncbi_taxdump_20240914'}\tPinned taxdump\n"
                    f"checkm2_db\tdatabase\tNA\t{alias_db_root / 'checkm2' / 'CheckM2_database'}\tCheckM2 database\n"
                    f"codetta_db\tdatabase\tNA\t{alias_db_root / 'codetta' / 'Pfam-A_enone'}\tCodetta profile database\n"
                    f"busco_datasets\tdatabase\tbacillota_odb12;mycoplasmatota_odb12\t{alias_db_root / 'busco'}\tConfigured BUSCO lineages in order\n"
                    f"bacillota_odb12\tdatabase\tNA\t{alias_db_root / 'busco' / 'bacillota_odb12'}\tBUSCO lineage dataset\n"
                    f"mycoplasmatota_odb12\tdatabase\tNA\t{alias_db_root / 'busco' / 'mycoplasmatota_odb12'}\tBUSCO lineage dataset\n"
                    f"eggnog_db\tdatabase\tNA\t{alias_db_root / 'Eggnog_db' / 'Eggnog_Diamond_db'}\teggNOG database\n"
                ),
            )

            validate_hpc_matrix.assert_versions_table_clean(versions_path, real_db_root)

    def test_medium_phylum_helper_allows_missing_metadata_case_na(self) -> None:
        """Allow `phylum=NA` only for the deliberate missing-metadata case."""
        master_rows = {
            "MYCO_NEW_1": {"phylum": "NA"},
            "ACC1": {"phylum": "Mycoplasmatota"},
            "ACC2": {"phylum": "Bacillota"},
        }

        validate_hpc_matrix.assert_medium_phyla_allowed(
            master_rows,
            allowed_phyla={"Mycoplasmatota", "Bacillota"},
            missing_metadata_case_accessions={"MYCO_NEW_1"},
        )

    def test_medium_phylum_helper_rejects_unexpected_na(self) -> None:
        """Reject `phylum=NA` for accessions outside the missing-metadata role."""
        master_rows = {
            "ACC_BAD": {"phylum": "NA"},
            "ACC1": {"phylum": "Mycoplasmatota"},
            "ACC2": {"phylum": "Bacillota"},
        }

        with self.assertRaises(validate_hpc_matrix.ValidateHpcMatrixError) as context:
            validate_hpc_matrix.assert_medium_phyla_allowed(
                master_rows,
                allowed_phyla={"Mycoplasmatota", "Bacillota"},
                missing_metadata_case_accessions=set(),
            )

        self.assertIn("Observed unexpected phylum=NA rows in medium run: ACC_BAD", str(context.exception))

    def test_medium_status_helper_rejects_second_failed_status(self) -> None:
        """Reject one medium row when any second failed status appears."""
        with tempfile.TemporaryDirectory() as tmpdir:
            status_path = Path(tmpdir) / "sample_status.tsv"
            write_text(
                status_path,
                (
                    "accession\tgcode_status\tprokka_status\twarnings\n"
                    "ACC1\tfailed\tfailed\tgcode_na\n"
                ),
            )

            with self.assertRaises(validate_hpc_matrix.ValidateHpcMatrixError) as context:
                validate_hpc_matrix.assert_medium_statuses_allowed(status_path)

            self.assertIn(
                "Found unexpected medium sample-status failures: ACC1 (gcode_status, prokka_status)",
                str(context.exception),
            )

    def test_medium_status_helper_rejects_missing_gcode_na_warning(self) -> None:
        """Reject one medium row when the gcode failure lacks the expected warning token."""
        with tempfile.TemporaryDirectory() as tmpdir:
            status_path = Path(tmpdir) / "sample_status.tsv"
            write_text(
                status_path,
                (
                    "accession\tgcode_status\twarnings\n"
                    "ACC1\tfailed\tcodetta_multiple_ncbi_tables\n"
                ),
            )

            with self.assertRaises(validate_hpc_matrix.ValidateHpcMatrixError) as context:
                validate_hpc_matrix.assert_medium_statuses_allowed(status_path)

            self.assertIn(
                "Found unexpected medium sample-status failures: ACC1 (gcode_status)",
                str(context.exception),
            )

    def test_medium_status_helper_rejects_primary_busco_failure(self) -> None:
        """Reject one medium row when the primary BUSCO lineage failed."""
        with tempfile.TemporaryDirectory() as tmpdir:
            status_path = Path(tmpdir) / "sample_status.tsv"
            write_text(
                status_path,
                (
                    "accession\tbusco_bacillota_odb12_status\twarnings\n"
                    "ACC1\tfailed\tbusco_summary_failed\n"
                ),
            )

            with self.assertRaises(validate_hpc_matrix.ValidateHpcMatrixError) as context:
                validate_hpc_matrix.assert_medium_statuses_allowed(status_path)

            self.assertIn(
                "Found unexpected medium sample-status failures: ACC1 (busco_bacillota_odb12_status)",
                str(context.exception),
            )

    def test_medium_status_helper_rejects_secondary_busco_failure_without_expected_warning(self) -> None:
        """Reject one medium row when a secondary BUSCO failure lacks the expected warning."""
        with tempfile.TemporaryDirectory() as tmpdir:
            status_path = Path(tmpdir) / "sample_status.tsv"
            write_text(
                status_path,
                (
                    "accession\tbusco_mycoplasmatota_odb12_status\twarnings\n"
                    "ACC1\tfailed\tmissing_busco_mycoplasmatota_odb12_summary\n"
                ),
            )

            with self.assertRaises(validate_hpc_matrix.ValidateHpcMatrixError) as context:
                validate_hpc_matrix.assert_medium_statuses_allowed(status_path)

            self.assertIn(
                "Found unexpected medium sample-status failures: ACC1 (busco_mycoplasmatota_odb12_status)",
                str(context.exception),
            )

    def test_strict_status_helper_rejects_any_failed_row(self) -> None:
        """Keep the tracked-run helper strict for any failed status row."""
        with tempfile.TemporaryDirectory() as tmpdir:
            status_path = Path(tmpdir) / "sample_status.tsv"
            write_text(
                status_path,
                (
                    "accession\tgcode_status\twarnings\n"
                    "ACC1\tfailed\tgcode_na\n"
                ),
            )

            with self.assertRaises(validate_hpc_matrix.ValidateHpcMatrixError) as context:
                validate_hpc_matrix.assert_no_failed_statuses(status_path)

            self.assertIn("Found failed sample-status rows: ACC1", str(context.exception))


if __name__ == "__main__":
    unittest.main()
