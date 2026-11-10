"""Tests for the manual pipeline test wrapper."""

from __future__ import annotations

import os
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "run_pipeline_test.sh"


class RunPipelineTestScriptTestCase(unittest.TestCase):
    """Exercise the Bash wrapper around the acceptance harness."""

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
        """Print wrapper help before a mode is chosen."""
        result = self.run_wrapper("--help")

        self.assertEqual(result.returncode, 0)
        self.assertIn(
            "bin/run_pipeline_test.sh [--dry-run] <prepare|unit|stub|local|slurm|dbprep-slurm|all>",
            result.stdout,
        )
        self.assertIn(
            "Manual wrapper around bin/run_acceptance_tests.py.", result.stdout
        )
        self.assertIn(
            "Host python3 must be >= 3.12 for the acceptance harness and validators.",
            result.stdout,
        )
        self.assertIn(
            "dbprep-slurm and all depend on the current runtime-db helper image and Codetta-aware helper CLIs.",
            result.stdout,
        )

    def test_wrapper_dry_run_stub_prints_delegated_command(self) -> None:
        """Print the delegated stub command without executing it."""
        result = self.run_wrapper("--dry-run", "stub")

        self.assertEqual(result.returncode, 0)
        self.assertEqual(
            result.stdout.strip(),
            "python3 bin/run_acceptance_tests.py stub",
        )
        self.assertEqual(result.stderr, "")

    def test_wrapper_dry_run_local_preserves_forwarded_arguments(self) -> None:
        """Forward real-run arguments unchanged during dry-run output."""
        result = self.run_wrapper(
            "--dry-run",
            "local",
            "--taxdump",
            "/tmp/taxdump",
            "--checkm2-db",
            "/tmp/checkm2",
            "--codetta-db",
            "/tmp/codetta",
            "--busco-db",
            "/tmp/busco",
            "--eggnog-db",
            "/tmp/eggnog",
        )

        self.assertEqual(result.returncode, 0)
        self.assertEqual(
            result.stdout.strip(),
            "python3 bin/run_acceptance_tests.py local --taxdump /tmp/taxdump "
            "--checkm2-db /tmp/checkm2 --codetta-db /tmp/codetta --busco-db /tmp/busco "
            "--eggnog-db /tmp/eggnog",
        )
        self.assertEqual(result.stderr, "")

    def test_wrapper_rejects_invalid_mode(self) -> None:
        """Fail with a clear mode-validation error."""
        result = self.run_wrapper("bogus")

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Error: unsupported mode bogus.", result.stderr)
        self.assertIn(
            "Expected one of: prepare unit stub local slurm dbprep-slurm all",
            result.stderr,
        )

    def test_wrapper_forwards_mode_help_to_acceptance_harness(self) -> None:
        """Delegate per-mode help text to the Python harness."""
        result = self.run_wrapper("local", "--help")

        self.assertEqual(result.returncode, 0)
        self.assertIn("Run the real-data local acceptance cohort.", result.stdout)
        self.assertIn("--checkm2-db CHECKM2_DB", result.stdout)
        self.assertIn("--codetta-db CODETTA_DB", result.stdout)

    def test_wrapper_slurm_help_shows_singularity_runtime_flags(self) -> None:
        """Expose renamed Singularity runtime flags through delegated help."""
        result = self.run_wrapper("slurm", "--help")

        self.assertEqual(result.returncode, 0)
        self.assertIn("--singularity-cache-dir SINGULARITY_CACHE_DIR", result.stdout)
        self.assertIn(
            "--singularity-run-options SINGULARITY_RUN_OPTIONS", result.stdout
        )
        self.assertIn("--slurm-qos SLURM_QOS", result.stdout)
        self.assertNotIn("--slurm-account", result.stdout)

    def test_wrapper_dry_run_slurm_preserves_singularity_arguments(self) -> None:
        """Forward renamed Singularity arguments unchanged during dry-run output."""
        result = self.run_wrapper(
            "--dry-run",
            "slurm",
            "--taxdump",
            "/tmp/taxdump",
            "--checkm2-db",
            "/tmp/checkm2",
            "--busco-db",
            "/tmp/busco",
            "--eggnog-db",
            "/tmp/eggnog",
            "--singularity-cache-dir",
            "/tmp/singularity-cache",
            "--singularity-run-options",
            "bind=/db",
            "--slurm-qos",
            "2h",
            "--slurm-cluster-options=--qos=2h",
        )

        self.assertEqual(result.returncode, 0)
        self.assertEqual(
            result.stdout.strip(),
            "python3 bin/run_acceptance_tests.py slurm --taxdump /tmp/taxdump "
            "--checkm2-db /tmp/checkm2 --busco-db /tmp/busco "
            "--eggnog-db /tmp/eggnog "
            "--singularity-cache-dir /tmp/singularity-cache "
            "--singularity-run-options bind=/db "
            "--slurm-qos 2h --slurm-cluster-options=--qos=2h",
        )
        self.assertEqual(result.stderr, "")

    def test_wrapper_dbprep_slurm_help_shows_prep_specific_flags(self) -> None:
        """Expose the dedicated database-prep SLURM mode through delegated help."""
        result = self.run_wrapper("dbprep-slurm", "--help")

        self.assertEqual(result.returncode, 0)
        self.assertIn("This mode runs prepare_databases.nf on SLURM", result.stdout)
        self.assertIn("--dbprep-profile DBPREP_PROFILE", result.stdout)
        self.assertIn("--codetta-db CODETTA_DB", result.stdout)
        self.assertIn("--busco-db BUSCO_DB", result.stdout)
        self.assertIn("--force-runtime-database-rebuild", result.stdout)

    def test_wrapper_dbprep_modes_guard_against_the_stale_helper_tag(self) -> None:
        """Keep the dbprep preflight tied to the refreshed helper image contract."""
        script_text = SCRIPT.read_text(encoding="utf-8")

        self.assertIn("preflight_runtime_db_helper_image()", script_text)
        self.assertIn("preflight_host_python()", script_text)
        self.assertIn('readonly MINIMUM_HOST_PYTHON="3.12"', script_text)
        self.assertIn(
            'readonly NEXTFLOW_CONFIG="${REPO_ROOT}/nextflow.config"', script_text
        )
        self.assertIn("params.runtime_db_helper_container", script_text)
        self.assertIn("quay.io/asuq1617/nf-myco_db:0.3", script_text)
        self.assertIn("quay.io/asuq1617/nf-myco_db:0.2", script_text)
        self.assertIn("dbprep-slurm | all)", script_text)

    def test_wrapper_rejects_old_host_python_before_running_the_harness(self) -> None:
        """Fail fast with a clear message when python3 on PATH is too old."""
        with tempfile.TemporaryDirectory() as tmpdir:
            fake_python = Path(tmpdir) / "python3"
            fake_python.write_text(
                "#!/usr/bin/env bash\n"
                'if [[ "$1" == "--version" ]]; then\n'
                "  printf 'Python 3.6.8\\n'\n"
                "  exit 0\n"
                "fi\n"
                'if [[ "$1" == "-c" ]]; then\n'
                "  exit 1\n"
                "fi\n"
                "exit 99\n",
                encoding="ascii",
            )
            fake_python.chmod(0o755)

            env = os.environ.copy()
            env["PATH"] = f"{tmpdir}{os.pathsep}{env['PATH']}"
            result = subprocess.run(
                ["bash", str(SCRIPT), "prepare"],
                cwd=ROOT,
                text=True,
                capture_output=True,
                check=False,
                env=env,
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("host python3 must be >= 3.12", result.stderr)
            self.assertIn("Python 3.6.8", result.stderr)

    def test_wrapper_dry_run_dbprep_slurm_preserves_runtime_arguments(self) -> None:
        """Forward dbprep-slurm arguments unchanged during dry-run output."""
        result = self.run_wrapper(
            "--dry-run",
            "dbprep-slurm",
            "--dbprep-profile",
            "oist",
            "--taxdump",
            "/tmp/taxdump",
            "--checkm2-db",
            "/tmp/checkm2",
            "--codetta-db",
            "/tmp/codetta",
            "--busco-db",
            "/tmp/busco",
            "--eggnog-db",
            "/tmp/eggnog",
            "--force-runtime-database-rebuild",
            "--slurm-queue",
            "short",
            "--slurm-qos",
            "2h",
            "--slurm-cluster-options=--constraint=ssd",
            "--singularity-cache-dir",
            "/tmp/singularity-cache",
        )

        self.assertEqual(result.returncode, 0)
        self.assertEqual(
            result.stdout.strip(),
            "python3 bin/run_acceptance_tests.py dbprep-slurm --dbprep-profile oist "
            "--taxdump /tmp/taxdump --checkm2-db /tmp/checkm2 --codetta-db /tmp/codetta --busco-db /tmp/busco "
            "--eggnog-db /tmp/eggnog "
            "--force-runtime-database-rebuild --slurm-queue short "
            "--slurm-qos 2h --slurm-cluster-options=--constraint=ssd "
            "--singularity-cache-dir /tmp/singularity-cache",
        )
        self.assertEqual(result.stderr, "")

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
