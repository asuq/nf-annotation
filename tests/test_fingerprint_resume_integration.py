"""Upgrade v0.4.5 bookkeeping without repeating a completed ANI calculation."""

import json
import shutil
import subprocess
import unittest

import test_cohort_update_integration as fixture


@unittest.skipUnless(
    fixture.NEXTFLOW and shutil.which("git"), "Nextflow and git are required"
)
class FingerprintResumeIntegrationTests(unittest.TestCase):
    def test_v045_upgrade_and_worker_changes_keep_completed_ani_cached(self):
        paths = (
            "bin/reuse_source_validation.py",
            "modules/local/reuse_source_validation.nf",
            "modules/local/prepare_cohort_update.nf",
            "subworkflows/local/cohort_ani.nf",
            "conf/base.config",
            "main.nf",
            "nextflow.config",
        )
        baseline = {}
        for path in paths:
            result = subprocess.run(
                ["git", "show", f"6e60ca05c963a1d5ab71042dd4680e277f0df245:{path}"],
                cwd=fixture.ROOT,
                capture_output=True,
                check=False,
            )
            if result.returncode:
                self.skipTest("The v0.4.5 commit is required for upgrade qualification")
            baseline[path] = result.stdout
        fixture.CohortUpdateIntegrationTestCase.setUpClass()
        self.addCleanup(fixture.CohortUpdateIntegrationTestCase.tearDownClass)
        case = fixture.CohortUpdateIntegrationTestCase(methodName="runTest")
        current = {path: (case.project / path).read_bytes() for path in paths}
        source = case.run_pipeline(
            "parallel-source", ["A", "B"], extra=["--pipeline_stage", "qc"]
        )
        for path, content in baseline.items():
            (case.project / path).write_bytes(content)
        args = ["--pipeline_stage", "ani", "--ani_16s_policy", "ignore"]
        output = case.run_pipeline(
            "parallel-view", ["A", "B"], source=source, extra=args
        )
        before = (output / "cohort/fastani/fastani.matrix").read_bytes()
        (case.root / "baseline-trace.tsv").write_bytes(
            (output / "pipeline_info/trace.tsv").read_bytes()
        )
        for path, content in current.items():
            (case.project / path).write_bytes(content)
        for requested in (4, 1):
            case.run_pipeline(
                "parallel-view",
                ["A", "B"],
                source=source,
                resume=True,
                extra=[
                    *args,
                    "--source_fingerprint_workers",
                    str(requested),
                    "--max_cpus",
                    "2",
                ],
            )
            trace = fixture.read_rows(output / "pipeline_info/trace.tsv")
            for process in ("FASTANI", "BUILD_FASTANI_INPUTS"):
                rows = [
                    row
                    for row in trace
                    if row["name"].split(" (")[0].endswith(":" + process)
                ]
                self.assertEqual([row["status"] for row in rows], ["CACHED"], process)
            self.assertEqual(
                before, (output / "cohort/fastani/fastani.matrix").read_bytes()
            )
            metrics = json.loads(
                (
                    output
                    / "pipeline_info/source_fingerprinting/fingerprint_metrics.json"
                ).read_text()
            )
            self.assertEqual(metrics["phase"], "complete")
            self.assertEqual(metrics["workers"], min(requested, 2))
            self.assertEqual(metrics["completed_roots"], 1)
            self.assertGreater(metrics["bytes_hashed"], 0)
