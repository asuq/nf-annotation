"""Explicit review data survives normal publication, resume and cohort reuse."""

import json
import sys
import unittest
from pathlib import Path

import test_cohort_update_integration as fixture

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from reviewed_codes import fasta_identity


@unittest.skipUnless(fixture.NEXTFLOW, "Nextflow is required")
class ReviewedCodeIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fixture.CohortUpdateIntegrationTestCase.setUpClass()
        cls.case = fixture.CohortUpdateIntegrationTestCase(methodName="runTest")

    @classmethod
    def tearDownClass(cls):
        fixture.CohortUpdateIntegrationTestCase.tearDownClass()

    def test_review_resume_and_reuse_use_the_same_decision(self):
        manifest = self.case.inputs / "reviews.json"
        review = {
            "accession": "A",
            "fasta_identity_md5": fasta_identity(self.case.genomes["A"]),
            "genetic_code": 4,
            "review_id": "synthetic",
            "reason": "integration control",
        }
        manifest.write_text(
            json.dumps({"schema": "reviewed-genetic-codes-v1", "reviews": [review]})
        )
        options = ["--pipeline_stage", "qc", "--reviewed_codes", str(manifest)]
        source = self.case.run_pipeline("review-source", ["A"], extra=options)
        master = fixture.read_rows(source / "tables/master_table.tsv")[0]
        self.assertEqual(json.loads(master["Gcode_Review"]), review)
        self.assertEqual(master["Gcode_Automatic"], "4")
        self.case.run_pipeline("review-source", ["A"], extra=options, resume=True)
        for row in self.case.heavy_tasks(source):
            self.assertEqual(row["status"], "CACHED")
        imported = self.case.run_pipeline(
            "review-reuse", ["A"], source=source, extra=options
        )
        self.assertEqual(
            fixture.read_rows(imported / "tables/master_table.tsv")[0]["Gcode_Review"],
            master["Gcode_Review"],
        )
        self.assertFalse(self.case.heavy_tasks(imported))
        self.case.run_pipeline(
            "review-missing",
            ["A"],
            source=source,
            extra=["--pipeline_stage", "qc"],
            expected_code=1,
        )
        review["reason"] = "changed review identity"
        manifest.write_text(
            json.dumps({"schema": "reviewed-genetic-codes-v1", "reviews": [review]})
        )
        self.case.run_pipeline(
            "review-changed", ["A"], source=source, extra=options, expected_code=1
        )
