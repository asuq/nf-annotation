"""Shared seed conservation and separate native annotation of each proteome."""

from __future__ import annotations

import copy
import csv
import json
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import normalize_eggnog as egg
import test_annotation_batch_tasks as fixture
from annotation_common import AnnotationError
from annotation_result import inventory
from eggnog_batches import validate_batch
from eggnog_native_fixture import native_annotations


class EggnogBatchNormalizationTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixture.AnnotationBatchTaskTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.raw, self.resource = self.fixture.raw, self.fixture.resource
        self.batch, self.proteins = validate_batch(self.fixture.batchdir)
        self.genes = {
            protein["accession"]: protein["gene_id"] for protein in self.proteins
        }
        self.queries = {
            protein["accession"]: protein["tool_id"] for protein in self.proteins
        }

    def normalize(self, proteins=None):
        return egg.normalize_batch(
            self.raw,
            self.proteins if proteins is None else proteins,
            self.resource,
            batch=self.batch,
        )

    def annotations(self, index):
        return (
            self.raw
            / "annotations"
            / f"member{index:08d}"
            / "eggnog.emapper.annotations"
        )

    def annotation_rows(self, index):
        lines = [
            line
            for line in self.annotations(index).read_text().splitlines()
            if not line.startswith("##")
        ]
        lines[0] = lines[0].removeprefix("#")
        return list(csv.DictReader(lines, delimiter="\t"))

    def test_native_tables_are_parsed_per_proteome_with_one_ontology_load(self):
        before = inventory(self.raw)
        with patch.object(
            egg, "_normalize_annotations", wraps=egg._normalize_annotations
        ) as parser, patch.object(
            egg, "ontology_terms", wraps=egg.ontology_terms
        ) as ontology:
            members = self.normalize()
        self.assertEqual(parser.call_count, 3)
        self.assertEqual(ontology.call_count, 1)
        self.assertEqual(list(members), ["A", "B", "C"])
        self.assertEqual(sum(member.reported_hits for member in members.values()), 2)
        self.assertEqual(members["A"].features["ko"], {self.genes["A"]: {"K00001"}})
        self.assertEqual(members["B"].features["ko"], {self.genes["B"]: {"K00002"}})
        self.assertEqual(
            {error["field"] for error in members["B"].errors}, {"go", "ec"}
        )
        self.assertEqual(members["A"].definitions, {"go": {}})
        self.assertEqual(members["B"].definitions, {})
        for accession, member in members.items():
            for _, rows in member.tables.values():
                self.assertTrue(all(row["accession"] == accession for row in rows))
        self.assertEqual(inventory(self.raw), before)

    def test_name_confidence_is_retained_even_when_primary_matrices_agree(self):
        # The native qualification exposed a medium/low Preferred_name change
        # that KO/GO/EC/Pfam-only comparisons did not detect.
        for index, confidence in enumerate(("m", "l")):
            rows = self.annotation_rows(index)
            rows[0]["Preferred_name"] = "hflB"
            rows[0]["annotation_confidence"] = (
                confidence + rows[0]["annotation_confidence"][1:]
            )
            native_annotations(self.annotations(index), rows, 1)
        members = self.normalize()
        self.assertEqual(
            members["A"].features["Preferred_name"], {self.genes["A"]: {"hflB"}}
        )
        self.assertNotIn("Preferred_name", members["B"].features)
        self.assertIn(self.genes["B"], members["B"].accepted)
        row = members["B"].tables["eggnog_annotations.tsv"][1][0]
        self.assertEqual(
            json.loads(row["excluded_fields"])["Preferred_name"], "low confidence"
        )

    def test_no_hit_member_has_header_only_tables_and_no_borrowed_evidence(self):
        empty = self.normalize()["C"]
        self.assertEqual(empty.reported_hits, 0)
        self.assertEqual(empty.mapped | empty.accepted, set())
        self.assertEqual(empty.features, {})
        self.assertEqual(empty.errors, [])
        self.assertEqual(empty.definitions, {})
        output = self.fixture.root / "no_hits"
        empty.publish(output)
        for name, columns in egg.TABLE_COLUMNS.items():
            self.assertEqual((output / name).read_text(), "\t".join(columns) + "\n")

    def test_go_definitions_remain_specific_to_each_native_proteome(self):
        for index, namespace, term in (
            (0, "gos_mf", "GO:0000001"),
            (1, "gos_bp", "GO:0000002"),
        ):
            rows = self.annotation_rows(index)
            rows[0]["GOs"] = term
            confidence = list(rows[0]["annotation_confidence"])
            confidence[egg.FIELDS.index("GOs")] = "h"
            rows[0]["annotation_confidence"] = "".join(confidence)
            native_annotations(self.annotations(index), rows, 1)
            namespace_row = dict.fromkeys(egg.GO_HEADER, "-")
            namespace_row.update(query=rows[0]["query"])
            namespace_row[namespace] = term
            namespace_row[namespace + "_confidence"] = "high"
            sidecar = self.annotations(index).with_name(
                "eggnog.emapper.annotations.go_namespaces.tsv"
            )
            sidecar.write_text(
                "\t".join(egg.GO_HEADER)
                + "\n"
                + "\t".join(namespace_row[name] for name in egg.GO_HEADER)
                + "\n"
            )
        members = self.normalize()
        self.assertEqual(members["A"].definitions, {"go": {"GO:0000001": "Control 1"}})
        self.assertEqual(members["B"].definitions, {"go": {"GO:0000002": "Control 2"}})
        self.assertEqual(members["C"].definitions, {})

    def test_seed_without_annotation_remains_mapped_and_unaccepted(self):
        native_annotations(self.annotations(1), [], 1)
        sidecar = self.annotations(1).with_name(
            "eggnog.emapper.annotations.go_namespaces.tsv"
        )
        sidecar.write_text("\t".join(egg.GO_HEADER) + "\n")
        member = self.normalize()["B"]
        self.assertEqual(member.reported_hits, 1)
        self.assertEqual(member.mapped, {self.genes["B"]})
        self.assertEqual(member.accepted, set())
        self.assertEqual(member.features, {})
        self.assertEqual(len(member.tables["eggnog_seed_hits.tsv"][1]), 1)
        self.assertEqual(member.tables["eggnog_annotations.tsv"][1], [])

    def test_original_protein_ids_cannot_shadow_canonical_queries(self):
        proteins = copy.deepcopy(self.proteins)
        proteins[0]["protein_id"] = self.queries["B"]
        members = self.normalize(proteins)
        self.assertEqual(members["A"].features["ko"], {self.genes["A"]: {"K00001"}})
        self.assertEqual(members["B"].features["ko"], {self.genes["B"]: {"K00002"}})

    def test_foreign_annotation_queries_and_unknown_seed_aliases_fail(self):
        rows = self.annotation_rows(0)
        rows[0]["query"] = self.queries["B"]
        native_annotations(self.annotations(0), rows, 1)
        with self.assertRaisesRegex(AnnotationError, "unknown protein ID"):
            self.normalize()
        seeds = self.raw / "search/eggnog.emapper.seed_orthologs"
        original = seeds.read_text()
        for value in ("unknown", "gene_1"):
            seeds.write_text(original.replace(self.queries["A"], value))
            with self.assertRaisesRegex(AnnotationError, "unknown protein ID"):
                self.normalize()

    def test_duplicate_or_empty_canonical_input_ids_fail(self):
        for key in ("tool_id", "gene_id"):
            proteins = copy.deepcopy(self.proteins)
            proteins[1][key] = proteins[0][key]
            with self.assertRaisesRegex(AnnotationError, "Duplicate canonical"):
                self.normalize(proteins)
        proteins = copy.deepcopy(self.proteins)
        proteins[0]["tool_id"] = ""
        with self.assertRaisesRegex(AnnotationError, "Invalid canonical"):
            self.normalize(proteins)
        with self.assertRaisesRegex(AnnotationError, "no declared proteins"):
            self.normalize([])


if __name__ == "__main__":
    unittest.main()
