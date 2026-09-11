"""Tests for the OIST HPC matrix wrapper."""

from __future__ import annotations

import os
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "run_oist_hpc_matrix.sh"


class RunOistHpcMatrixScriptTestCase(unittest.TestCase):
    """Exercise the OIST HPC matrix wrapper."""

    def make_ready_hpc_root(self) -> Path:
        """Create one temporary HPC root with ready markers for all DBs."""
        root = Path(tempfile.mkdtemp(prefix="oist_hpc_matrix_"))
        ready_roots = (
            root / "db" / "ncbi_taxdump_20240914",
            root / "db" / "checkm2" / "CheckM2_database",
            root / "db" / "codetta" / "Pfam-A_enone",
            root / "db" / "busco",
            root / "db" / "Eggnog_db" / "Eggnog_Diamond_db",
        )
        for ready_root in ready_roots:
            ready_root.mkdir(parents=True, exist_ok=True)
            (ready_root / ".nf_myco_ready.json").write_text("{}", encoding="ascii")
        return root

    def run_wrapper(self, *args: str) -> subprocess.CompletedProcess[str]:
        """Run the wrapper with text capture enabled."""
        return subprocess.run(
            ["bash", str(SCRIPT), *args],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_wrapper_help_lists_supported_modes(self) -> None:
        """Print wrapper help for the HPC matrix runner."""
        result = self.run_wrapper("--help")

        self.assertEqual(result.returncode, 0)
        self.assertIn(
            "<prepare|medium-prepare|db-full|db-reuse|db-matrix|p1|p2|all>",
            result.stdout,
        )
        self.assertIn("Run the OIST HPC validation campaign", result.stdout)
        self.assertIn("--hpc-root PATH", result.stdout)
        self.assertIn("--gcode-rule RULE", result.stdout)
        self.assertIn(
            "Host python3 must be >= 3.12 for the harness and matrix validators.",
            result.stdout,
        )
        self.assertNotIn("--medium-candidates-tsv", result.stdout)

    def test_prepare_dry_run_prints_cohort_command(self) -> None:
        """Print the tracked-cohort preparation command."""
        result = self.run_wrapper(
            "--dry-run",
            "--hpc-root",
            "/tmp/hpc",
            "prepare",
        )

        self.assertEqual(result.returncode, 0)
        self.assertIn(
            "INFO: [prepare] Running tracked 9-sample cohort prepare",
            result.stdout,
        )
        self.assertIn(
            "bin/run_pipeline_test.sh prepare --work-root /tmp/hpc/acceptance",
            result.stdout,
        )

    def test_medium_prepare_dry_run_prints_medium_prepare_command(self) -> None:
        """Print the fixed medium-cohort preparation command."""
        result = self.run_wrapper(
            "--dry-run",
            "--hpc-root",
            "/tmp/hpc",
            "medium-prepare",
        )

        self.assertEqual(result.returncode, 0)
        self.assertIn(
            "INFO: [medium-prepare] Running fixed medium Mycoplasmatota/Bacillota cohort prepare",
            result.stdout,
        )
        self.assertIn(
            "python3 bin/run_acceptance_tests.py prepare --work-root /tmp/hpc/medium",
            result.stdout,
        )
        self.assertIn(
            f"--source-catalog {ROOT}/assets/tables/medium/source_catalog.tsv",
            result.stdout,
        )
        self.assertIn(
            f"--cohort-plan {ROOT}/assets/tables/medium/cohort_plan.tsv",
            result.stdout,
        )

    def test_db_full_dry_run_uses_wrapper_and_validator(self) -> None:
        """Validate a fresh DB root as one prepared first-pass run."""
        result = self.run_wrapper(
            "--dry-run",
            "--hpc-root",
            "/tmp/hpc",
            "db-full",
        )

        self.assertEqual(result.returncode, 0)
        self.assertIn(
            "INFO: [db-full] Running full runtime DB prep gate",
            result.stdout,
        )
        self.assertIn("bin/run_pipeline_test.sh dbprep-slurm", result.stdout)
        self.assertIn("--dbprep-profile oist", result.stdout)
        self.assertIn("--codetta-db /tmp/hpc/db/codetta/Pfam-A_enone", result.stdout)
        self.assertIn("python3 bin/validate_hpc_matrix.py dbprep", result.stdout)
        self.assertIn("--expected-component codetta", result.stdout)
        self.assertIn("--expected-arg=--codetta_db", result.stdout)
        self.assertIn("--expected-status prepared", result.stdout)
        self.assertNotIn("--padloc-db", result.stdout)

    def test_db_full_dry_run_uses_present_for_ready_root(self) -> None:
        """Validate an already-ready DB root as a reuse pass."""
        hpc_root = self.make_ready_hpc_root()
        result = self.run_wrapper(
            "--dry-run",
            "--hpc-root",
            str(hpc_root),
            "db-full",
        )

        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertIn("--expected-status present", result.stdout)

    def test_db_full_dry_run_rejects_partially_ready_root(self) -> None:
        """Fail fast when only some canonical DB roots are ready."""
        hpc_root = Path(tempfile.mkdtemp(prefix="oist_hpc_matrix_partial_"))
        ready_root = hpc_root / "db" / "ncbi_taxdump_20240914"
        ready_root.mkdir(parents=True, exist_ok=True)
        (ready_root / ".nf_myco_ready.json").write_text("{}", encoding="ascii")

        result = self.run_wrapper(
            "--dry-run",
            "--hpc-root",
            str(hpc_root),
            "db-full",
        )

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Inconsistent runtime DB root", result.stderr)

    def test_p1_dry_run_uses_profile_defaults(self) -> None:
        """Avoid explicit max-resource overrides in the real tracked run."""
        result = self.run_wrapper(
            "--dry-run",
            "--hpc-root",
            "/tmp/hpc",
            "p1",
        )

        self.assertEqual(result.returncode, 0)
        self.assertIn(
            "INFO: [p1] Running tracked 9-sample edge-case pipeline test",
            result.stdout,
        )
        self.assertIn("nextflow run . -profile oist", result.stdout)
        self.assertIn("--sample_csv /tmp/hpc/acceptance/generated/sample_sheet.csv", result.stdout)
        self.assertIn("--codetta_db /tmp/hpc/db/codetta/Pfam-A_enone", result.stdout)
        self.assertIn("--eggnog_db /tmp/hpc/db/Eggnog_db/Eggnog_Diamond_db", result.stdout)
        self.assertNotIn("--max_cpus", result.stdout)
        self.assertNotIn("--max_memory", result.stdout)
        self.assertNotIn("--max_time", result.stdout)
        self.assertNotIn("--padloc_db", result.stdout)

    def test_p1_dry_run_forwards_gcode_rule_override(self) -> None:
        """Pass an explicit gcode rule override through to the pipeline run."""
        result = self.run_wrapper(
            "--dry-run",
            "--hpc-root",
            "/tmp/hpc",
            "--gcode-rule",
            "delta_then_11",
            "p1",
        )

        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertIn("--gcode_rule delta_then_11", result.stdout)

    def test_p2_auto_prepares_medium_inputs(self) -> None:
        """Refresh the fixed medium cohort inputs before each p2 run."""
        result = self.run_wrapper(
            "--dry-run",
            "--hpc-root",
            "/tmp/hpc",
            "p2",
        )

        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertIn(
            "INFO: [p2] Running medium Mycoplasmatota/Bacillota pipeline test",
            result.stdout,
        )
        self.assertIn(
            "python3 bin/run_acceptance_tests.py prepare --work-root /tmp/hpc/medium",
            result.stdout,
        )
        self.assertIn(
            "--sample_csv /tmp/hpc/medium/generated/sample_sheet.csv",
            result.stdout,
        )
        self.assertIn(
            "--metadata /tmp/hpc/medium/generated/metadata.tsv",
            result.stdout,
        )
        self.assertIn(
            "--metadata-tsv /tmp/hpc/medium/generated/metadata.tsv",
            result.stdout,
        )
        self.assertIn(
            f"--cohort-plan {ROOT}/assets/tables/medium/cohort_plan.tsv",
            result.stdout,
        )
        self.assertIn(
            f"--source-catalog {ROOT}/assets/tables/medium/source_catalog.tsv",
            result.stdout,
        )

    def test_fixed_medium_inputs_are_not_reused_implicitly(self) -> None:
        """Refresh the tracked medium cohort instead of trusting stale generated files."""
        script_text = SCRIPT.read_text(encoding="utf-8")

        self.assertIn('MEDIUM_SAMPLE_CSV="${MEDIUM_ROOT}/generated/sample_sheet.csv"', script_text)
        self.assertIn('MEDIUM_METADATA="${MEDIUM_ROOT}/generated/metadata.tsv"', script_text)
        self.assertIn("run_medium_prepare", script_text)
        self.assertNotIn('! -f "${MEDIUM_SAMPLE_CSV}"', script_text)
        self.assertNotIn('! -f "${MEDIUM_METADATA}"', script_text)

    def test_p2_override_inputs_keep_looser_validator_contract(self) -> None:
        """Do not force fixed-cohort validation arguments for custom medium inputs."""
        result = self.run_wrapper(
            "--dry-run",
            "--hpc-root",
            "/tmp/hpc",
            "p2",
            "--medium-sample-csv",
            "/tmp/medium_sample_sheet.csv",
            "--medium-metadata",
            "/tmp/medium_metadata.tsv",
        )

        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertIn("--sample_csv /tmp/medium_sample_sheet.csv", result.stdout)
        self.assertIn("--metadata /tmp/medium_metadata.tsv", result.stdout)
        self.assertNotIn("--metadata-tsv /tmp/medium_metadata.tsv", result.stdout)
        self.assertNotIn("--cohort-plan", result.stdout)
        self.assertNotIn("--source-catalog", result.stdout)

    def test_all_auto_prepares_medium_inputs(self) -> None:
        """Run the medium prepare step automatically during the full campaign."""
        result = self.run_wrapper(
            "--dry-run",
            "--hpc-root",
            "/tmp/hpc",
            "all",
        )

        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertIn(
            "python3 bin/run_acceptance_tests.py prepare --work-root /tmp/hpc/medium",
            result.stdout,
        )

    def test_db_matrix_uses_taxdump_helper_failure_text(self) -> None:
        """Keep taxdump negative-case expectations aligned with the helper."""
        script_text = SCRIPT.read_text(encoding="utf-8")

        self.assertIn(
            "No local source was supplied for taxdump, and remote download is disabled.",
            script_text,
        )
        self.assertIn(
            "Destination must be a directory for taxdump",
            script_text,
        )

    def test_db_matrix_uses_codetta_missing_marker_failure_text(self) -> None:
        """Keep the Codetta marker-less case aligned with the helper contract."""
        script_text = SCRIPT.read_text(encoding="utf-8")

        self.assertIn(
            'fail_log="${valid_root}/codetta_missing_marker.log"',
            script_text,
        )
        self.assertIn(
            "Destination is valid but incomplete for codetta",
            script_text,
        )

    def test_db_matrix_dry_run_resets_disposable_case_root(self) -> None:
        """Reset only the disposable db_cases tree before seeding matrix cases."""
        result = self.run_wrapper(
            "--dry-run",
            "--hpc-root",
            "/tmp/hpc",
            "db-matrix",
        )
        resolved_case_root = str(Path("/tmp/hpc/db_cases").resolve())

        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertIn(f"rm -rf {resolved_case_root}", result.stdout)
        self.assertIn(f"mkdir -p {resolved_case_root}", result.stdout)
        self.assertIn(
            "INFO: [db-matrix] Running db3_valid_without_marker/checkm2",
            result.stdout,
        )
        self.assertIn(
            "INFO: [db-matrix] Running db3_valid_without_marker/codetta",
            result.stdout,
        )
        self.assertIn(
            "--codetta_db /tmp/hpc/db_cases/db3_valid_without_marker/codetta",
            result.stdout,
        )
        self.assertNotIn(
            "--results-dir /tmp/hpc/db_cases/db3_valid_without_marker/results_codetta",
            result.stdout,
        )
        self.assertIn(
            "INFO: [db-matrix] Running db5_taxdump_missing_no_download",
            result.stdout,
        )
        self.assertIn(
            "INFO: [db-matrix] Running db5_codetta_missing_no_download",
            result.stdout,
        )
        self.assertIn(
            "INFO: [db-matrix] Running db10_checkm2_nested_layout",
            result.stdout,
        )

    def test_db_matrix_reset_targets_db_cases_only(self) -> None:
        """Guard the reset so it cannot target the wider HPC root."""
        script_text = SCRIPT.read_text(encoding="utf-8")

        self.assertIn("preflight_host_python()", script_text)
        self.assertIn('readonly MINIMUM_HOST_PYTHON="3.12"', script_text)
        self.assertIn('resolved_case_root="$(python3 -c', script_text)
        self.assertIn('if [[ "${resolved_case_root}" == "/" ]]; then', script_text)
        self.assertIn('if [[ "${resolved_case_root}" != */db_cases ]]; then', script_text)

    def test_wrapper_rejects_old_host_python_before_running_matrix_steps(self) -> None:
        """Fail fast with a clear message when python3 on PATH is too old."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fake_python = Path(tmpdir) / "python3"
            fake_python.write_text(
                "#!/usr/bin/env bash\n"
                "if [[ \"$1\" == \"--version\" ]]; then\n"
                "  printf 'Python 3.6.8\\n'\n"
                "  exit 0\n"
                "fi\n"
                "if [[ \"$1\" == \"-c\" ]]; then\n"
                "  exit 1\n"
                "fi\n"
                "exit 99\n",
                encoding="ascii",
            )
            fake_python.chmod(0o755)

            env = os.environ.copy()
            env["PATH"] = f"{tmpdir}{os.pathsep}{env['PATH']}"
            result = subprocess.run(
                ["bash", str(SCRIPT), "--hpc-root", "/tmp/hpc", "prepare"],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
                env=env,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Host python3 must be >= 3.12", result.stderr)
            self.assertIn("Python 3.6.8", result.stderr)

    def test_all_accepts_medium_inputs_after_mode(self) -> None:
        """Accept medium inputs even when they appear after the mode token."""
        result = self.run_wrapper(
            "--dry-run",
            "--hpc-root",
            "/tmp/hpc",
            "all",
            "--medium-sample-csv",
            "/tmp/medium_sample_sheet.csv",
            "--medium-metadata",
            "/tmp/medium_metadata.tsv",
        )

        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertIn(
            "--sample_csv /tmp/medium_sample_sheet.csv",
            result.stdout,
        )
        self.assertIn(
            "--metadata /tmp/medium_metadata.tsv",
            result.stdout,
        )

    def test_wrapper_has_valid_bash_syntax(self) -> None:
        """Pass Bash syntax validation."""
        result = subprocess.run(
            ["bash", "-n", str(SCRIPT)],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )

        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")
        self.assertEqual(result.stderr, "")


if __name__ == "__main__":
    unittest.main()
