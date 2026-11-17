"""A short genome retains native prediction evidence when excluded from analyses."""

import hashlib
import unittest

import test_cohort_update_integration as fixture


@unittest.skipUnless(fixture.NEXTFLOW, "Nextflow required")
class ShortAssemblyIntegrationTests(unittest.TestCase):
    def test_post_prediction_length_boundary(self):
        fixture.CohortUpdateIntegrationTestCase.setUpClass()
        self.addCleanup(fixture.CohortUpdateIntegrationTestCase.tearDownClass)
        case = fixture.CohortUpdateIntegrationTestCase()
        for name in ("A", "B", "FAILED"):
            path = case.genomes[name]
            seq = path.read_text().splitlines()[1][:19999]
            # At the cutoff the last base belongs to a contig Prokka filters.
            tail = ">small\nC\n" if name == "B" else ""
            path.write_text(f">contig_{name}\n{seq}\n{tail}")
            case.input_fingerprints[name] = hashlib.sha256(
                path.read_bytes()
            ).hexdigest()
        root = case.run_pipeline(
            "short-qc", ["A", "B", "FAILED"], extra=["--pipeline_stage", "qc"]
        )
        rows = fixture.read_rows(root / "tables/master_table.tsv")
        self.assertEqual([row["accession"] for row in rows], ["A", "B", "FAILED"])
        self.assertEqual(
            [row["Assembly_Bases"] for row in rows], ["19999", "20000", "19999"]
        )
        self.assertEqual(
            [row["Assembly_Exclusion_Reason"] for row in rows],
            ["assembly_below_20000_bp", "NA", "assembly_below_20000_bp"],
        )
        self.assertEqual(rows[0]["CheckM2_Low_quality"], "false")
        self.assertEqual(rows[0]["Low_quality"], "true")
        statuses = fixture.read_rows(root / "tables/sample_status.tsv")
        self.assertEqual(
            [row["prokka_status"] for row in statuses], ["done", "done", "failed"]
        )
        intact = fixture.read_rows(root / "cohort/16s/all_best_16S_manifest.tsv")
        self.assertEqual([r["accession"] for r in intact], ["B"])
        clustered = case.run_pipeline(
            "short-clusters",
            ["A", "B", "FAILED"],
            extra=["--pipeline_stage", "clustering"],
            source=root,
        )
        rows = fixture.read_rows(clustered / "tables/master_table.tsv")
        self.assertEqual([r["Cluster_ID"] == "NA" for r in rows], [True, False, True])
        self.assertEqual(
            [r["AAI_Cluster_ID"] == "NA" for r in rows], [True, False, True]
        )
