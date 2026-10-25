"""Mixed-domain QC publications and cumulative ANI through real Nextflow."""

import json
import shutil
import unittest

import test_cohort_update_integration as fixture


@unittest.skipUnless(fixture.NEXTFLOW, "Nextflow is required")
class DomainStageIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        fixture.CohortUpdateIntegrationTestCase.setUpClass()
        cls.case = fixture.CohortUpdateIntegrationTestCase(methodName="runTest")
        cls.root = cls.case.root
        cls.config = cls.root / "domains.config"
        cls.config.write_text(
            "params { busco_lineages = []; domain = 'metadata'; busco_domain_lineages = [Bacteria: ['bacteria_odb12'], Archaea: ['archaea_odb12']] }\n"
        )
        lines = cls.case.metadata.read_text().splitlines()
        cls.case.metadata.write_text(
            lines[0]
            + "\tDomain\n"
            + "".join(
                line + ("\tArchaea\n" if line.split("\t")[0] == "B" else "\tBacteria\n")
                for line in lines[1:]
            )
        )
        for lineage in ("bacteria_odb12", "archaea_odb12"):
            target = cls.case.project / "assets/fixtures/stub/busco" / lineage
            target.mkdir()
            (target / "dataset.cfg").write_text(f"name={lineage}\n")
        path = cls.case.project / "tests/cohort_update_fixture.py"
        text = path.read_text().replace(
            '"300" if args.translation_table == "4" else "150"',
            '"300" if args.translation_table == "4" and accession != "B" else "150"',
        )
        text = text.replace(
            "translate(table=4)", 'translate(table=11 if accession == "B" else 4)'
        ).replace(
            '"transl_table": ["4"]',
            '"transl_table": ["11" if accession == "B" else "4"]',
        )
        text = (
            text.replace(
                '"C": 98.0,', '"C": 92.0 if args.lineage == "archaea_odb12" else 98.0,'
            )
            .replace(
                '"S": 98.0,', '"S": 92.0 if args.lineage == "archaea_odb12" else 98.0,'
            )
            .replace(
                '"M": 1.0,', '"M": 7.0 if args.lineage == "archaea_odb12" else 1.0,'
            )
        )
        path.write_text(text)

    @classmethod
    def tearDownClass(cls):
        fixture.CohortUpdateIntegrationTestCase.tearDownClass()

    def run_stage(self, name, accessions, stage, **kwargs):
        return self.case.run_pipeline(
            name,
            accessions,
            domain=None,
            extra=[
                "-c",
                str(self.config),
                "--pipeline_stage",
                stage,
                "--ani_16s_policy",
                "ignore",
                "--annotation_tools",
                "eggnog,cogclassifier,pfam,kofam,padloc",
            ],
            **kwargs,
        )

    def test_qc_then_multi_source_ani_preserves_domains_codes_and_evidence(self):
        shared_config = self.config.read_text()
        self.config.write_text(
            "params { domain = 'metadata'; busco_domain_lineages = [Bacteria: ['bacteria_odb12']] }\n"
        )
        a = self.run_stage("domain-qc-a", ["A", "LOW", "NOGCODE"], "qc")
        self.config.write_text(
            "params { domain = 'metadata'; busco_domain_lineages = [Archaea: ['archaea_odb12']] }\n"
        )
        b = self.run_stage("domain-qc-b", ["B"], "qc")
        self.config.write_text(shared_config)
        for root in (a, b):
            names = {
                r["name"].split(" (")[0].split(":")[-1]
                for r in fixture.read_rows(root / "pipeline_info/trace.tsv")
            }
            self.assertFalse(
                names
                & {
                    "FASTANI",
                    "CLUSTER_ANI",
                    "SELECT_ANI_REPRESENTATIVES",
                    "CCFINDER",
                    "ANNOTATION_SEARCH",
                }
            )
            self.assertTrue(
                json.loads((root / "annotation_results.json").read_text())["complete"]
            )
            self.assertTrue(
                all(
                    r["ani_included"] == "na"
                    for r in fixture.read_rows(root / "tables/sample_status.tsv")
                )
            )
        self.assertEqual(
            {p.name for p in (a / "samples/A/busco").iterdir()}, {"bacteria_odb12"}
        )
        self.assertEqual(
            {p.name for p in (b / "samples/B/busco").iterdir()}, {"archaea_odb12"}
        )
        shutil.rmtree(self.root / "work-domain-qc-a")
        shutil.rmtree(self.root / "work-domain-qc-b")
        output = self.run_stage(
            "domain-ani", ["A", "B", "LOW", "NOGCODE"], "ani", source=[a, b]
        )
        names = {
            r["name"].split(" (")[0].split(":")[-1]
            for r in fixture.read_rows(output / "pipeline_info/trace.tsv")
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
        self.assertIn("FASTANI", names)
        master = {
            r["accession"]: r
            for r in fixture.read_rows(output / "tables/master_table.tsv")
        }
        status = {
            r["accession"]: r
            for r in fixture.read_rows(output / "tables/sample_status.tsv")
        }
        self.assertEqual((master["A"]["Gcode"], master["B"]["Gcode"]), ("4", "11"))
        self.assertEqual(master["B"]["Primary_BUSCO_Lineage"], "archaea_odb12")
        self.assertEqual(
            master["B"]["Primary_BUSCO_Score"], master["B"]["BUSCO_archaea_odb12"]
        )
        self.assertTrue(master["B"]["Primary_BUSCO_Score"].startswith("C:92.0%"))
        self.assertEqual(master["B"]["BUSCO_bacteria_odb12"], "NA")
        self.assertEqual(status["B"]["busco_bacteria_odb12_status"], "na")
        self.assertEqual(status["B"]["ani_included"], "true")
        self.assertEqual(status["LOW"]["ani_included"], "false")
        self.assertIn("low_quality", status["LOW"]["ani_exclusion_reason"])
        self.assertEqual(status["NOGCODE"]["prokka_status"], "skipped")
        self.assertEqual(status["NOGCODE"]["ani_included"], "false")
        ani = {
            r["accession"]: r
            for r in fixture.read_rows(output / "cohort/fastani/ani_metadata.tsv")
        }
        self.assertEqual(ani["B"]["BUSCO_selected"], master["B"]["BUSCO_archaea_odb12"])
        self.assertEqual(ani["B"]["primary_busco_lineage"], "archaea_odb12")
        for accession, source in (("A", a), ("B", b)):
            self.assertEqual(
                (source / f"samples/{accession}/prokka/prokka.faa").read_bytes(),
                (output / f"samples/{accession}/prokka/prokka.faa").read_bytes(),
            )
        self.run_stage(
            "domain-ani-missing", ["A", "D"], "ani", source=a, expected_code=1
        )
        self.assertIn(
            "completed QC publications for every requested sample",
            (self.root / "launch-domain-ani-missing/captured.log").read_text(),
        )
        original_config = self.config.read_text()
        self.config.write_text(
            original_config.replace("['archaea_odb12']", "['other_archaeal_odb12']")
        )
        self.run_stage("domain-incompatible", ["B"], "ani", source=b, expected_code=1)
        self.config.write_text(original_config)

    def test_single_domain_modes_produce_matching_provenance(self):
        original = self.config.read_text()
        try:
            for mode, accession, lineage in (
                ("bacteria", "A", "bacteria_odb12"),
                ("archaea", "B", "archaea_odb12"),
            ):
                self.config.write_text(
                    original.replace("domain = 'metadata'", f"domain = '{mode}'")
                )
                output = self.run_stage("single-" + mode, [accession], "qc")
                row = fixture.read_rows(output / "tables/master_table.tsv")[0]
                self.assertEqual(row["Annotation_Domain"], mode.capitalize())
                self.assertEqual(row["Annotation_Domain_Source"], "domain")
                self.assertEqual(row["Primary_BUSCO_Lineage"], lineage)
        finally:
            self.config.write_text(original)

    def test_missing_and_conflicting_metadata_fail_before_qc(self):
        original = self.case.metadata.read_text()
        try:
            self.case.metadata.write_text(original.replace("\tArchaea\n", "\tNA\n"))
            self.run_stage("domain-missing", ["B"], "qc", expected_code=1)
            log = (self.root / "launch-domain-missing/captured.log").read_text()
            self.assertIn("Sample 'B': Domain must be Bacteria or Archaea", log)
            self.assertNotIn("PER_SAMPLE_QC:", log)
        finally:
            self.case.metadata.write_text(original)
        original_config = self.config.read_text()
        try:
            self.config.write_text(original_config + "params.domain = 'bacteria'\n")
            self.run_stage("domain-conflict", ["B"], "qc", expected_code=1)
            self.assertIn(
                "Domain metadata conflicts",
                (self.root / "launch-domain-conflict/captured.log").read_text(),
            )
        finally:
            self.config.write_text(original_config)

    def test_legacy_source_cannot_acquire_an_unrecorded_domain_policy(self):
        legacy = self.case.run_pipeline("legacy-source", ["A"])
        original = self.config.read_text()
        try:
            self.config.write_text(
                "params { domain = 'metadata'; busco_domain_lineages = [Bacteria: ['bacillota_odb12', 'mycoplasmatota_odb12']] }\n"
            )
            self.run_stage(
                "legacy-unsafe", ["A"], "ani", source=legacy, expected_code=1
            )
            self.assertIn(
                "Domain/BUSCO policy changed or is unrecorded",
                (self.root / "launch-legacy-unsafe/captured.log").read_text(),
            )
        finally:
            self.config.write_text(original)

    def test_changed_policy_code_is_revalidated_on_resume(self):
        self.run_stage("policy-resume", ["A"], "qc")
        manifest = self.case.inputs / "policy-resume.csv"
        old_method = self.case.manifest
        helper = self.case.project / "bin/sample_policy.py"
        original = helper.read_text()
        try:
            self.case.manifest = lambda name, accessions: manifest
            # Change only the imported helper; inputs and the invoking CLI stay fixed.
            helper.write_text(
                original
                + "\n_original_resolve = resolve_sample\ndef resolve_sample(*args, **kwargs):\n    raise ValueError('changed policy control')\n"
            )
            self.run_stage("policy-resume", ["A"], "qc", resume=True, expected_code=1)
            self.assertIn(
                "changed policy control",
                (self.root / "launch-policy-resume/captured.log").read_text(),
            )
        finally:
            self.case.manifest = old_method
            helper.write_text(original)


if __name__ == "__main__":
    unittest.main()
