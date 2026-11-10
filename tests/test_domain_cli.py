"""Reject ambiguous domain configuration before any workflow task starts."""

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")


@unittest.skipUnless(NEXTFLOW, "Nextflow is required")
class DomainCliTests(unittest.TestCase):
    def test_missing_invalid_and_removed_options_fail_before_tasks(self):
        cases = [
            (
                "main.nf",
                ["--domain", "bacteria", "--pipeline_stage", "ani"],
                "has been removed",
            ),
            (
                "main.nf",
                ["--domain", "bacteria", "--ani_publication_mode", "reference"],
                "has been removed",
            ),
            (
                "main.nf",
                ["--domain", "bacteria", "--ani_import_batches", "8"],
                "has been removed",
            ),
            (
                "main.nf",
                ["--domain", "bacteria", "--ani_import_max_forks", "2"],
                "has been removed",
            ),
            (
                "main.nf",
                ["--domain", "bacteria", "--ani_score_profile", "mag"],
                "has been removed",
            ),
            (
                "main.nf",
                ["-profile", "test", "--domain", "bacteria", "--aai_threshold", "0.91"],
                "--aai_threshold",
            ),
            (
                "reannotate.nf",
                ["--source_fingerprint_workers", "9"],
                "must be an integer from 1 to 8",
            ),
            (
                "main.nf",
                [
                    "--domain",
                    "bacteria",
                    "--pipeline_stage",
                    "clustering",
                    "--publication_mode",
                    "reference",
                ],
                "--qc_publication_root is required",
            ),
            (
                "main.nf",
                ["--domain", "bacteria", "--cohort_validation_workers", "9"],
                "must be an integer from 1 to 8",
            ),
            (
                "main.nf",
                ["--domain", "bacteria", "--source_fingerprint_workers", "9"],
                "must be an integer from 1 to 8",
            ),
            (
                "main.nf",
                ["--domain", "bacteria", "--source_fingerprint_workers", "0"],
                "must be an integer from 1 to 8",
            ),
            (
                "main.nf",
                ["--domain", "bacteria", "--publication_mode", "unknown"],
                "must be self_contained or reference",
            ),
            (
                "main.nf",
                ["--domain", "bacteria", "--publication_mode", "reference"],
                "requires --pipeline_stage clustering",
            ),
            ("main.nf", [], "--domain is required"),
            ("main.nf", ["--domain", "auto"], "--domain is required"),
            ("main.nf", ["--pipeline_stage", "qc"], "--domain is required"),
            ("main.nf", ["--pipeline_stage", "clustering"], "--domain is required"),
        ]
        for entrypoint in ("main.nf", "prepare_databases.nf", "reannotate.nf"):
            for name in ("domain_column", "prokka_kingdom", "barrnap_kingdom"):
                cases.append(
                    (
                        entrypoint,
                        ["--" + name, "unused"],
                        "--" + name + " is unsupported",
                    )
                )
        for entrypoint, arguments, expected in cases:
            with (
                self.subTest(entrypoint=entrypoint, arguments=arguments),
                tempfile.TemporaryDirectory(dir="/tmp") as name,
            ):
                result = subprocess.run(
                    [NEXTFLOW, "run", str(ROOT / entrypoint), *arguments],
                    cwd=name,
                    capture_output=True,
                    text=True,
                    check=False,
                    timeout=90,
                    env=dict(
                        os.environ,
                        NXF_ANSI_LOG="false",
                        NXF_DISABLE_CHECK_LATEST="true",
                    ),
                )
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(expected, result.stdout + result.stderr)
                self.assertNotIn("Submitted process", result.stdout + result.stderr)

    def test_non_analysis_entrypoints_do_not_require_domain(self):
        # Reach the existing input validation without scheduling a task.
        for entrypoint, expected in (
            ("reannotate.nf", "--annotation_from requires"),
            ("prepare_databases.nf", "destination"),
        ):
            with (
                self.subTest(entrypoint=entrypoint),
                tempfile.TemporaryDirectory(dir="/tmp") as name,
            ):
                result = subprocess.run(
                    [NEXTFLOW, "run", str(ROOT / entrypoint)],
                    cwd=name,
                    capture_output=True,
                    text=True,
                    check=False,
                    timeout=90,
                    env=dict(
                        os.environ,
                        NXF_ANSI_LOG="false",
                        NXF_DISABLE_CHECK_LATEST="true",
                    ),
                )
                log = result.stdout + result.stderr
                self.assertNotEqual(result.returncode, 0)
                self.assertNotIn("--domain is required", log)
                self.assertIn(expected, log)
                self.assertNotIn("Submitted process", log)
