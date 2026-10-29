"""Frozen v0.4.4 QC remains untouched through relocated v2 reuse and reannotation."""

import io
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import test_cohort_update_integration as fixture
from annotation_common import digest, read_json
from annotation_source import validate_source
from fingerprinting.source_fingerprint import fingerprint

BASELINE = "6b8c0b886f7a49c83b08f3cefbc58d44c392136c"


@unittest.skipUnless(
    fixture.NEXTFLOW and shutil.which("git"), "Nextflow and git are required"
)
class V044RelocationIntegrationTests(unittest.TestCase):
    def assert_no_qc(self, output):
        names = {
            row["name"].split(" (")[0].split(":")[-1]
            for row in fixture.read_rows(output / "pipeline_info/trace.tsv")
        }
        self.assertFalse(
            names
            & {
                "BARRNAP",
                "CHECKM2_GCODE4",
                "CHECKM2_GCODE11",
                "BUSCO",
                "PROKKA",
                "CODETTA",
                "CCFINDER",
                "ANNOTATION_SEARCH",
            }
        )

    def test_frozen_qc_move_reuse_materialise_and_reannotate(self):
        original = fixture.ROOT
        archive = subprocess.run(
            ["git", "archive", BASELINE], cwd=original, capture_output=True, check=False
        )
        if archive.returncode:
            self.skipTest("Frozen v0.4.4 commit required")
        temporary = tempfile.TemporaryDirectory(dir="/tmp")
        self.addCleanup(temporary.cleanup)
        baseline = Path(temporary.name) / "baseline"
        baseline.mkdir()
        with tarfile.open(fileobj=io.BytesIO(archive.stdout)) as bundle:
            bundle.extractall(baseline, filter="data")
        shutil.copytree(
            original / "external/nf-helper/conf", baseline / "external/nf-helper/conf"
        )
        fixture.ROOT = baseline
        try:
            fixture.CohortUpdateIntegrationTestCase.setUpClass()
        finally:
            fixture.ROOT = original
        old = fixture.CohortUpdateIntegrationTestCase(methodName="runTest")
        old_root = old.root
        self.addCleanup(lambda: shutil.rmtree(old_root, ignore_errors=True))
        accessions = ["A", "B", "LOW", "NOGCODE"]
        qc = old.run_pipeline("v044-qc", accessions, extra=["--pipeline_stage", "qc"])
        expected = {
            p.relative_to(qc).as_posix(): digest(p)
            for p in qc.rglob("*")
            if p.is_file()
        }
        qc_sha = fingerprint(qc)["sha256"]
        fixture.CohortUpdateIntegrationTestCase.setUpClass()
        self.addCleanup(fixture.CohortUpdateIntegrationTestCase.tearDownClass)
        case = fixture.CohortUpdateIntegrationTestCase(methodName="runTest")
        collection = case.root / "collection"
        boundary = collection / "02-qc"
        boundary.mkdir(parents=True)
        provider = boundary / "batch-a"
        shutil.move(qc, provider)
        shutil.rmtree(old_root)
        for accession in accessions:
            case.genomes[accession] = (
                provider / f"samples/{accession}/staged/{accession}.fasta"
            )
        case.input_fingerprints = {name: digest(path) for name, path in case.genomes.items()}
        view = collection / "03-ani/view"
        args = ["--pipeline_stage", "ani", "--ani_16s_policy", "ignore"]
        case.run_pipeline(
            "relative-view",
            accessions,
            source=provider,
            output=view,
            extra=[
                *args,
                "--ani_publication_mode",
                "reference",
                "--qc_publication_root",
                str(boundary),
            ],
        )
        self.assert_no_qc(view)
        before = fingerprint(view, workers=1)["sha256"]
        document = (view / "ani_references.json").read_bytes()
        self.assertEqual(
            expected,
            {
                p.relative_to(provider).as_posix(): digest(p)
                for p in provider.rglob("*")
                if p.is_file()
            },
        )
        # Same-location resume changes only allocation, and must keep completed ANI cached.
        case.run_pipeline(
            "relative-view",
            accessions,
            source=provider,
            output=view,
            resume=True,
            extra=[
                *args,
                "--ani_publication_mode",
                "reference",
                "--qc_publication_root",
                str(boundary),
                "--cohort_validation_workers",
                "1",
                "--source_fingerprint_workers",
                "1",
            ],
        )
        trace = fixture.read_rows(view / "pipeline_info/trace.tsv")
        for process in ("FASTANI", "BUILD_FASTANI_INPUTS"):
            self.assertEqual(
                [
                    r["status"]
                    for r in trace
                    if r["name"].split(" (")[0].endswith(":" + process)
                ],
                ["CACHED"],
            )
        shutil.rmtree(case.root / "work-relative-view")
        shutil.rmtree(case.root / "launch-relative-view")
        moved = Path(temporary.name) / "another parent" / "collection"
        moved.parent.mkdir()
        shutil.copytree(collection, moved, symlinks=True)
        shutil.rmtree(collection)
        provider = moved / "02-qc/batch-a"
        view = moved / "03-ani/view"
        validate_source(view)
        self.assertEqual(before, fingerprint(view, workers=4)["sha256"])
        self.assertEqual(qc_sha, fingerprint(provider)["sha256"])
        self.assertEqual(document, (view / "ani_references.json").read_bytes())
        for accession in accessions:
            case.genomes[accession] = (
                provider / f"samples/{accession}/staged/{accession}.fasta"
            )
        successor = moved / "03-ani/successor"
        case.run_pipeline(
            "successor",
            accessions,
            source=view,
            output=successor,
            extra=[
                *args,
                "--ani_publication_mode",
                "reference",
                "--qc_publication_root",
                str(moved / "02-qc"),
            ],
        )
        self.assert_no_qc(successor)
        shutil.rmtree(view)
        validate_source(successor)
        materialised = case.run_pipeline(
            "materialised", accessions, source=successor, extra=args
        )
        self.assert_no_qc(materialised)
        shutil.copy2(original / "reannotate.nf", case.project / "reannotate.nf")
        output = case.root / "reannotated"
        launch = case.root / "reannotation-launch"
        launch.mkdir()
        result = subprocess.run(
            [
                fixture.NEXTFLOW,
                "run",
                str(case.project / "reannotate.nf"),
                "-profile",
                "test",
                "--annotation_from",
                str(successor),
                "--outdir",
                str(output),
                "-work-dir",
                str(case.root / "annotation-work"),
            ],
            cwd=launch,
            capture_output=True,
            text=True,
            check=False,
            timeout=180,
            env=dict(os.environ, NXF_DISABLE_CHECK_LATEST="true", NXF_ANSI_LOG="false"),
        )
        (launch / "captured.log").write_text(result.stdout + result.stderr)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assert_no_qc(output)
        self.assertEqual(
            expected,
            {
                p.relative_to(provider).as_posix(): digest(p)
                for p in provider.rglob("*")
                if p.is_file()
            },
        )
        shutil.rmtree(moved)
        for independent in (materialised, output):
            validate_source(independent)
            self.assertFalse(
                any(p.is_symlink() for p in (independent / "samples").rglob("*"))
            )
            self.assertTrue(
                read_json(independent / "annotation_results.json")["complete"]
            )
