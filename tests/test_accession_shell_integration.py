"""Original identifiers remain literal data through native wrapper/report boundaries."""

import csv
import hashlib
import unittest
import test_cohort_update_integration as fixture


@unittest.skipUnless(fixture.NEXTFLOW, "Nextflow required")
class AccessionShellIntegrationTests(unittest.TestCase):
    def test_original_accessions_are_quoted_and_preserved(self):
        fixture.CohortUpdateIntegrationTestCase.setUpClass()
        self.addCleanup(fixture.CohortUpdateIntegrationTestCase.tearDownClass)
        case = fixture.CohortUpdateIntegrationTestCase(methodName="runTest")
        names = ["O'Brien", "A$(touch injected)", "-sample"]
        with case.metadata.open("a", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            for index, name in enumerate(names):
                genome = case.inputs / f"literal{index}.fasta"
                genome.write_bytes(case.genomes["A"].read_bytes())
                case.genomes[name] = genome
                case.input_fingerprints[name] = hashlib.sha256(
                    genome.read_bytes()
                ).hexdigest()
                writer.writerow([name, "8", "Synthetic bacterium", "Scaffold", "NA"])
        result = case.run_pipeline(
            "literal-identifiers", names, extra=["--pipeline_stage", "qc"]
        )
        masters = fixture.read_rows(result / "tables/master_table.tsv")
        self.assertEqual({row["accession"] for row in masters}, set(names))
        self.assertTrue(all(row["Gcode"] == "4" for row in masters))
        self.assertFalse(list(case.root.rglob("injected")))
        rows = fixture.read_rows(result / "tables/accession_map.tsv")
        for row in rows:
            self.assertRegex(row["internal_id"], r"^[A-Za-z0-9_]+$")
