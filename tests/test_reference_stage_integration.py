"""Reference views, chained ANI reuse and independent reannotation through Nextflow."""

import os
import shutil
import subprocess
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import test_cohort_update_integration as fixture
from annotation_common import read_json
from annotation_source import validate_source


@unittest.skipUnless(fixture.NEXTFLOW, "Nextflow is required")
class ReferenceStageIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fixture.CohortUpdateIntegrationTestCase.setUpClass()
        cls.case = fixture.CohortUpdateIntegrationTestCase(methodName="runTest")
        cls.root = cls.case.root
        shutil.copy2(fixture.ROOT / "reannotate.nf", cls.case.project / "reannotate.nf")

    @classmethod
    def tearDownClass(cls):
        fixture.CohortUpdateIntegrationTestCase.tearDownClass()

    def ani(self, name, source, mode):
        return self.case.run_pipeline(
            name,
            ["A", "LOW", "NOGCODE"],
            source=source,
            extra=[
                "--pipeline_stage",
                "ani",
                "--ani_16s_policy",
                "ignore",
                "--ani_publication_mode",
                mode,
            ],
        )

    def test_reference_views_reuse_and_materialise_without_changing_science(self):
        source = self.case.run_pipeline(
            "ref-source", ["A", "LOW", "NOGCODE"], extra=["--pipeline_stage", "qc"]
        )
        portable = self.ani("copy-view", source, "self_contained")
        reference = self.ani("reference-view", source, "reference")
        self.assertEqual(
            read_json(reference / "annotation_results.json")["pipeline_contract"],
            "nf-annotation-ani-reference-v1",
        )
        validate_source(reference)
        matched = 0
        for path in (portable / "tables").rglob("*.tsv"):
            if path.name == "tool_and_db_versions.tsv":
                continue
            self.assertEqual(
                path.read_bytes(),
                (reference / path.relative_to(portable)).read_bytes(),
                path.name,
            )
            matched += 1
        self.assertGreaterEqual(matched, 20)
        audits = [
            read_json(path)
            for path in (reference / "pipeline_info/ani_import").glob("*.json")
        ]
        self.assertTrue(
            all(
                row["files_copied"] == row["files_cloned"] == 0
                for audit in audits
                for row in audit["samples"]
            )
        )
        chained = self.ani("chained-view", reference, "reference")
        providers = read_json(chained / "ani_references.json")["sources"]
        self.assertEqual([row["root"] for row in providers], [str(source.resolve())])
        shutil.rmtree(reference)
        validate_source(chained)
        materialised = self.ani("materialised-view", chained, "self_contained")
        reannotation = self.root / "reannotation"
        launch = self.root / "reannotation-launch"
        launch.mkdir()
        result = subprocess.run(
            [
                fixture.NEXTFLOW,
                "run",
                str(self.case.project / "reannotate.nf"),
                "-profile",
                "test",
                "--annotation_from",
                str(chained),
                "--outdir",
                str(reannotation),
                "-work-dir",
                str(self.root / "reannotation-work"),
            ],
            cwd=launch,
            check=False,
            capture_output=True,
            text=True,
            timeout=180,
            env=dict(os.environ, NXF_ANSI_LOG="false", NXF_DISABLE_CHECK_LATEST="true"),
        )
        (launch / "captured.log").write_text(result.stdout + result.stderr)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        validate_source(reannotation)
        shutil.rmtree(chained)
        shutil.rmtree(source)
        for output in (portable, materialised, reannotation):
            validate_source(output)
            self.assertFalse(
                any(path.is_symlink() for path in (output / "samples").rglob("*"))
            )
