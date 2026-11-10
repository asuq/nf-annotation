"""Both clustering methods share QC while AAI eligibility and native reuse stay independent."""

import shutil
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import test_cohort_update_integration as fixture
from annotation_common import read_json, digest


@unittest.skipUnless(fixture.NEXTFLOW, "Nextflow required")
class AaiStageIntegrationTests(unittest.TestCase):
    def test_marker_preparation_starts_before_the_busco_cohort_finishes(self):
        fixture.CohortUpdateIntegrationTestCase.setUpClass()
        self.addCleanup(fixture.CohortUpdateIntegrationTestCase.tearDownClass)
        case = fixture.CohortUpdateIntegrationTestCase(methodName="runTest")
        marker = case.root / "aai-started"
        caller = case.project / "tests/cohort_update_fixture.py"
        source = caller.read_text()
        gate = (
            "    if args.tool == 'busco' and args.accession == 'B':\n"
            "        import time\n"
            "        deadline = time.monotonic() + 45\n"
            f"        while not Path({str(marker)!r}).exists():\n"
            "            if time.monotonic() > deadline: raise AssertionError('AAI waited for BUSCO')\n"
            "            time.sleep(0.05)\n"
        )
        caller.write_text(
            source.replace(
                "def make_stub(args: argparse.Namespace) -> None:\n",
                "def make_stub(args: argparse.Namespace) -> None:\n" + gate,
            )
        )
        driver = case.project / "tests/fastaai_workflow_fixture.py"
        driver.write_text(
            driver.read_text().replace(
                "    prepare_aai.prepare(args)",
                "    prepare_aai.prepare(args)\n    Path("
                + repr(str(marker))
                + ").write_text('ready')",
            )
        )
        result = case.run_pipeline("ready-aai", ["A", "B"])
        self.assertTrue(marker.is_file())
        trace = fixture.read_rows(result / "pipeline_info/trace.tsv")
        self.assertEqual(sum(":PREPARE_AAI_BATCH" in row["name"] for row in trace), 2)
        self.assertTrue((result / "cohort/aai_clusters/aai_summary.tsv").is_file())

    def test_16s_independence_and_threshold_score_changes_reuse_native_work(self):
        fixture.CohortUpdateIntegrationTestCase.setUpClass()
        self.addCleanup(fixture.CohortUpdateIntegrationTestCase.tearDownClass)
        case = fixture.CohortUpdateIntegrationTestCase(methodName="runTest")
        with (case.project / "nextflow.config").open("a") as handle:
            handle.write(
                "\ntrace.fields = 'task_id,hash,name,status,exit,cpus,memory'\ntrace.raw = true\n"
            )
        path = case.inputs / "PARTIAL.fasta"
        path.write_bytes(case.genomes["B"].read_bytes())
        case.genomes["PARTIAL"] = path
        case.input_fingerprints["PARTIAL"] = digest(path)
        with case.metadata.open("a") as handle:
            handle.write("PARTIAL\t8\tSynthetic partial\tScaffold\tNA\n")
        samples = ["A", "PARTIAL"]
        qc = case.run_pipeline("aai-qc", samples, extra=["--pipeline_stage", "qc"])
        names = {
            r["name"].split(" (")[0].split(":")[-1]
            for r in fixture.read_rows(qc / "pipeline_info/trace.tsv")
        }
        self.assertFalse(
            names & {"PREPARE_AAI_BATCH", "FASTAAI_COMPARE", "CLUSTER_AAI"}
        )
        args = ["--pipeline_stage", "clustering"]
        output = case.run_pipeline("aai-view", samples, source=qc, extra=args)
        master = {
            r["accession"]: r
            for r in fixture.read_rows(output / "tables/master_table.tsv")
        }
        self.assertEqual(master["PARTIAL"]["Cluster_ID"], "NA")
        self.assertNotEqual(master["PARTIAL"]["AAI_Cluster_ID"], "NA")
        self.assertEqual(
            master["A"]["AAI_Cluster_ID"], master["PARTIAL"]["AAI_Cluster_ID"]
        )
        self.assertTrue((output / "cohort/fastaai/Bacteria/jaccard.npy").is_file())
        trace = fixture.read_rows(output / "pipeline_info/trace.tsv")
        clusters = [
            r for r in trace if r["name"].split(" (")[0].endswith(":CLUSTER_AAI")
        ]
        self.assertEqual(len(clusters), 1)
        self.assertEqual(int(clusters[0]["cpus"]), 1)
        self.assertEqual(int(clusters[0]["memory"]), 3 * 1024**3)

        for extra in (
            ["--aai_threshold", ".85"],
            ["--aai_threshold", ".85", "--score_profile", "mag"],
        ):
            case.run_pipeline(
                "aai-view", samples, source=qc, extra=[*args, *extra], resume=True
            )
            trace = fixture.read_rows(output / "pipeline_info/trace.tsv")
            for process in ("FASTANI", "FASTAAI_COMPARE"):
                self.assertEqual(
                    [
                        r["status"]
                        for r in trace
                        if r["name"].split(" (")[0].endswith(":" + process)
                    ],
                    ["CACHED"],
                    process,
                )
            metrics = [
                read_json(p)
                for p in (output / "pipeline_info/fastaai").glob("batch_*.json")
            ]
            self.assertEqual(sum(m["extracted"] for m in metrics), 0)
            reps = fixture.read_rows(
                output / "cohort/aai_clusters/aai_representatives.tsv"
            )
            self.assertEqual(len(reps), 2)
        # Reuse the published profiles in an independently materialised successor.
        successor = case.run_pipeline(
            "aai-successor", samples, source=output, extra=args
        )
        metrics = [
            read_json(p)
            for p in (successor / "pipeline_info/fastaai").glob("batch_*.json")
        ]
        self.assertEqual(sum(m["extracted"] for m in metrics), 0)
        shutil.rmtree(output)
        for p in (successor / "cohort/fastaai/profiles").glob("*/profile.json"):
            self.assertFalse(p.is_symlink())
            self.assertTrue((p.parent / "markers.faa").is_file())


if __name__ == "__main__":
    unittest.main()
