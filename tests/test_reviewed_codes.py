"""Reviewed code choices cannot escape their genome identity or paired evidence."""

import gzip
import json
import shutil
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import reviewed_codes as reviews
import summarise_checkm2 as qc
import test_summarise_checkm2 as fixture


class ReviewedCodeTests(unittest.TestCase):
    root: Path
    setUp = fixture.SummariseCheckM2TestCase.setUp
    write_report = fixture.SummariseCheckM2TestCase.write_report

    def test_empty_summary_header_matches_live_contract(self):
        header = (
            Path(__file__).resolve().parents[1]
            / "assets/tables/headers/checkm2_summary.tsv"
        )
        self.assertEqual(
            header.read_text().rstrip("\n").split("\t"), list(qc.OUTPUT_COLUMNS)
        )

    def entry(self):
        self.genome = self.root / "genome.fasta"
        self.genome.write_text(">one description\nATGTGA\n>two\nAcgTN\n")
        return {
            "accession": "ACC1",
            "fasta_identity_md5": reviews.fasta_identity(self.genome),
            "genetic_code": 11,
            "review_id": "synthetic-review",
            "reason": "explicit test decision",
        }

    def summary(self, entries, **kwargs):
        return qc.build_output_row(
            "ACC1",
            qc.parse_report(self.write_report(4)),
            qc.parse_report(self.write_report(11, Completeness="30")),
            warnings=[],
            genome=self.genome,
            genome_id="sample",
            review_candidates=entries,
            **kwargs,
        )

    def test_effective_call_retains_automatic_call_and_recomputes_quality(self):
        entry = self.entry()
        row = self.summary([entry])
        self.assertEqual(
            (row["Gcode_Automatic"], row["Gcode"], row["Low_quality"]),
            ("4", "11", "true"),
        )
        self.assertEqual(row["Gcode_Automatic_Reason"], "length_ratio_above_threshold")
        self.assertEqual(json.loads(row["Gcode_Review"]), entry)
        self.assertEqual(row["Gcode_Selection_Reason"], "reviewed_code_manifest")
        unreviewed = self.summary([])
        self.assertEqual((unreviewed["Gcode"], unreviewed["Gcode_Review"]), ("4", "NA"))
        self.assertEqual(self.summary([entry | {"genetic_code": 4}])["Gcode"], "4")

    def test_fasta_identity_preserves_order_ids_case_and_dna_not_wrapping(self):
        entry = self.entry()
        self.genome.write_text(">one changed description\nATG\nTGA\n>two\nAcgTN\n")
        self.assertEqual(
            reviews.fasta_identity(self.genome), entry["fasta_identity_md5"]
        )
        zipped = self.root / "compressed"
        zipped.write_bytes(gzip.compress(self.genome.read_bytes()))
        self.assertEqual(reviews.fasta_identity(zipped), entry["fasta_identity_md5"])
        for data in (
            ">one\nATGTGA\n>two\nACGTN\n",
            ">two\nAcgTN\n>one\nATGTGA\n",
            ">changed\nATGTGA\n>two\nAcgTN\n",
        ):
            self.genome.write_text(data)
            with self.assertRaisesRegex(ValueError, "identity"):
                reviews.select_review("ACC1", self.genome, [entry])

    def test_invalid_pair_foreign_review_or_missing_identity_cannot_override(self):
        entry = self.entry()
        report = qc.parse_report(self.write_report(4))
        for report11 in (
            None,
            report,
            qc.parse_report(self.write_report(11, Name="foreign")),
        ):
            with self.assertRaisesRegex(ValueError, "paired"):
                qc.build_output_row(
                    "ACC1",
                    report,
                    report11,
                    warnings=[],
                    genome=self.genome,
                    review_candidates=[entry],
                )
        for candidate in (
            entry | {"accession": "FOREIGN"},
            entry | {"fasta_identity_md5": "0" * 32},
        ):
            with self.assertRaises(ValueError):
                self.summary([candidate])

    def test_manifest_rejects_ambiguous_or_malformed_choices(self):
        entry = self.entry()
        path = self.root / "reviews.json"
        for entries in (
            [entry, entry],
            [entry | {"genetic_code": True}],
            [entry | {"genetic_code": 5}],
            [entry | {"reason": ""}],
        ):
            path.write_text(
                json.dumps({"schema": "reviewed-genetic-codes-v1", "reviews": entries})
            )
            with self.assertRaises(ValueError):
                reviews.load_manifest(path)
        path.write_text(
            json.dumps({"schema": "reviewed-genetic-codes-v1", "reviews": [entry]})
        )
        self.assertEqual(reviews.load_manifest(path), {"ACC1": [entry]})

    def test_publication_replays_choice_and_rejects_changed_evidence(self):
        entry = self.entry()
        master = self.summary([entry])
        sample = {
            "accession": "ACC1",
            "internal_id": "sample",
            "review_candidates_base64": reviews.encode_candidates([entry]),
        }
        (self.root / "staged").mkdir()
        shutil.copy2(self.genome, self.root / "staged/sample.fasta")
        for code in (4, 11):
            directory = self.root / f"checkm2_gcode{code}"
            directory.mkdir()
            shutil.copy2(
                self.root / f"gcode{code}.tsv", directory / "quality_report.tsv"
            )
        reviews.validate_published_decision(sample, self.root, master)
        for changed in (
            master | {"Gcode_Automatic": "11"},
            master | {"Gcode": "4"},
            master | {"Low_quality": "false"},
        ):
            with self.assertRaisesRegex(ValueError, "disagrees"):
                reviews.validate_published_decision(sample, self.root, changed)
        with self.assertRaises(ValueError):
            reviews.validate_published_decision(
                sample | {"review_candidates_base64": reviews.encode_candidates([])},
                self.root,
                master,
            )
        (self.root / "staged/sample.fasta").write_text(">different\nATGTGA\n")
        with self.assertRaisesRegex(ValueError, "identity"):
            reviews.validate_published_decision(sample, self.root, master)
