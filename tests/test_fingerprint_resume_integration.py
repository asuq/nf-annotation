"""Worker changes reuse completed v0.4.8 computation and genome reporting."""

import json
import shutil
import unittest

import test_cohort_update_integration as fixture


@unittest.skipUnless(
    fixture.NEXTFLOW and shutil.which("git"), "Nextflow and git are required"
)
class FingerprintResumeIntegrationTests(unittest.TestCase):
    def test_worker_changes_keep_completed_ani_and_reporting_cached(self):
        fixture.CohortUpdateIntegrationTestCase.setUpClass()
        self.addCleanup(fixture.CohortUpdateIntegrationTestCase.tearDownClass)
        case = fixture.CohortUpdateIntegrationTestCase(methodName="runTest")
        source = case.run_pipeline(
            "parallel-source", ["A", "B"], extra=["--pipeline_stage", "qc"]
        )
        args = ["--pipeline_stage", "ani", "--ani_16s_policy", "ignore"]
        output = case.run_pipeline(
            "parallel-view", ["A", "B"], source=source, extra=args
        )
        before = (output / "cohort/fastani/fastani.matrix").read_bytes()
        (case.root / "baseline-trace.tsv").write_bytes(
            (output / "pipeline_info/trace.tsv").read_bytes()
        )
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
            for process in ("FASTANI", "BUILD_FASTANI_INPUTS", "REPORT_GENOMES"):
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
