"""Conservation and sample isolation when projecting one native eggNOG batch."""

from __future__ import annotations

import copy
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import normalize_eggnog as egg
import test_annotation_normalizers as fixture
from annotation_common import AnnotationError


class EggnogBatchNormalizationTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixture.NormalizerTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.tearDown)
        self.raw, self.resource = self.fixture.raw, self.fixture.resource
        self.proteins = [
            dict(
                accession=accession,
                gene_id=f"{accession}::gene{index}",
                protein_id=f"gene{index}",
                tool_id=f"p{accession}{index}",
                length="200",
            )
            for accession, indices in (("A", (1, 2)), ("B", (1, 2)), ("C", (1,)))
            for index in indices
        ]
        self.rows = [
            self.annotation("pB1", GOs="GO:0000003", EC="ec:bad", KEGG_ko="K00003"),
            self.annotation("pA2", GOs="GO:0000002", KEGG_ko="K00002"),
            self.annotation(
                "pA1", GOs="GO:0000001", EC="ec:1.2.3.4", PFAMs="Control_Pfam"
            ),
        ]
        self.go_rows = [
            self.namespaces("pB1", gos_mf="GO:0000003"),
            self.namespaces("pA2", gos_bp="GO:0000002"),
            self.namespaces("pA1", gos_mf="GO:0000001"),
        ]
        self.fixture.write_eggnog(self.rows, self.go_rows)
        # B2 has a validated native seed but no functional annotation row.
        path = self.raw / "eggnog.emapper.seed_orthologs"
        text = path.read_text().replace("## 3 queries scanned\n", "")
        path.write_text(
            text
            + "pB2\t1\t0\t100\t1\t100\t1\t100\t90\t50\t50\n"
            + "## 4 queries scanned\n"
        )

    def annotation(self, identifier, **fields):
        row = dict.fromkeys(egg.HEADER, "-")
        row.update(
            query=identifier,
            seed_ortholog="1234",
            evalue="0",
            score="100",
            COG_category="J",
            **fields,
        )
        row["annotation_confidence"] = "".join(
            "h" if name in fields else "-" for name in egg.FIELDS
        )
        return row

    def namespaces(self, identifier, **fields):
        row = dict.fromkeys(egg.GO_HEADER, "-")
        row.update(query=identifier, **fields)
        for name in fields:
            row[name + "_confidence"] = "high"
        return row

    def test_one_parse_conserves_every_evidence_row_and_assignment(self):
        complete = egg.normalize(self.raw, self.proteins, self.resource)
        with patch.object(egg, "normalize", wraps=egg.normalize) as parser:
            members = egg.normalize_batch(self.raw, self.proteins, self.resource)
        self.assertEqual(parser.call_count, 1)
        self.assertEqual(list(members), ["A", "B", "C"])
        self.assertEqual(
            sum(member.reported_hits for member in members.values()),
            complete.reported_hits,
        )
        self.assertEqual(
            set().union(*(member.mapped for member in members.values())),
            complete.mapped,
        )
        self.assertEqual(
            set().union(*(member.accepted for member in members.values())),
            complete.accepted,
        )
        for name, (columns, rows) in complete.tables.items():
            self.assertEqual(
                sum(len(member.tables[name][1]) for member in members.values()),
                len(rows),
            )
            for accession, member in members.items():
                self.assertEqual(
                    member.tables[name],
                    (columns, [row for row in rows if row["accession"] == accession]),
                )
        self.assertEqual(
            {
                name: {
                    gene: values
                    for member in members.values()
                    for gene, values in member.features.get(name, {}).items()
                }
                for name in complete.features
            },
            complete.features,
        )
        self.assertEqual(members["A"].errors, [])
        self.assertEqual(members["C"].errors, [])
        self.assertEqual(members["B"].errors, complete.errors)
        self.assertEqual(
            {error["field"] for error in members["B"].errors}, {"ec", "go"}
        )
        self.assertEqual(
            members["A"].definitions["go"],
            {"GO:0000001": "Control 1", "GO:0000002": "Control 2"},
        )
        self.assertEqual(members["B"].definitions, {})
        self.assertEqual(members["C"].definitions, {})
        self.assertEqual(members["B"].reported_hits, 2)
        self.assertIn("B::gene2", members["B"].mapped)
        self.assertNotIn("B::gene2", members["B"].accepted)
        self.assertEqual(members["B"].features["ko"], {"B::gene1": {"K00003"}})

    def test_no_hit_member_has_header_only_tables_and_no_borrowed_evidence(self):
        members = egg.normalize_batch(self.raw, self.proteins, self.resource)
        empty = members["C"]
        self.assertEqual(empty.reported_hits, 0)
        self.assertEqual(empty.mapped | empty.accepted, set())
        self.assertEqual(empty.features, {})
        self.assertEqual(empty.errors, [])
        output = self.fixture.root / "no_hits"
        empty.publish(output)
        for name, columns in egg.TABLE_COLUMNS.items():
            self.assertEqual((output / name).read_text(), "\t".join(columns) + "\n")
        self.assertEqual(
            {path.name for path in output.iterdir()}, set(egg.TABLE_COLUMNS)
        )

    def test_original_protein_ids_never_shadow_another_canonical_query(self):
        proteins = copy.deepcopy(self.proteins)
        proteins[0]["protein_id"] = "pB1"
        members = egg.normalize_batch(self.raw, proteins, self.resource)
        self.assertEqual(members["B"].features["ko"], {"B::gene1": {"K00003"}})
        self.assertEqual(members["A"].features["ko"], {"A::gene2": {"K00002"}})

    def test_unknown_native_queries_and_original_aliases_fail(self):
        path = self.raw / "eggnog.emapper.seed_orthologs"
        original = path.read_text()
        for value in ("unknown", "gene1"):
            path.write_text(original.replace("pA1", value))
            with self.assertRaisesRegex(AnnotationError, "unknown protein ID"):
                egg.normalize_batch(self.raw, self.proteins, self.resource)

    def test_duplicate_or_empty_canonical_input_ids_fail(self):
        for key in ("tool_id", "gene_id"):
            proteins = copy.deepcopy(self.proteins)
            proteins[1][key] = proteins[0][key]
            with self.assertRaisesRegex(AnnotationError, "Duplicate canonical"):
                egg.normalize_batch(self.raw, proteins, self.resource)
        proteins = copy.deepcopy(self.proteins)
        proteins[0]["tool_id"] = ""
        with self.assertRaisesRegex(AnnotationError, "Invalid canonical"):
            egg.normalize(self.raw, proteins, self.resource)
        with self.assertRaisesRegex(AnnotationError, "no declared proteins"):
            egg.normalize_batch(self.raw, [], self.resource)


if __name__ == "__main__":
    unittest.main()
