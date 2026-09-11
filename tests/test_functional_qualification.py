"""Protect independent qualification checks against misleading coordinate/count joins."""

import copy
import unittest
from collections import Counter
from decimal import Decimal

from qualify_functional_cohort import (
    AnnotationError,
    category_counts,
    crispr_coordinates,
    matrix_counts,
)


class FunctionalQualificationTests(unittest.TestCase):
    def payload(self):
        return {
            "Sequences": [
                {
                    "Id": "contig",
                    "Length": 100,
                    "Crisprs": [
                        {"Evidence_Level": 1},
                        {
                            "Name": "array1",
                            "Evidence_Level": 2,
                            "Start": 11,
                            "End": 40,
                            "Spacers": 2,
                        },
                    ],
                }
            ],
        }

    def test_native_contig_alias_and_established_evidence_policy(self):
        rows = crispr_coordinates(self.payload(), {"contig.2": 100})
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["source_contig_id"], "contig.2")
        self.assertEqual(
            (rows[0]["start"], rows[0]["end"], rows[0]["spacer_count"]), (11, 40, 2)
        )

    def test_colliding_version_stripped_contig_names_are_rejected(self):
        with self.assertRaisesRegex(AnnotationError, "Ambiguous"):
            crispr_coordinates(self.payload(), {"contig.1": 100, "contig.2": 100})

    def test_native_coordinates_require_complete_unique_bounded_contigs(self):
        original = self.payload()
        variants = []
        wrong_length = copy.deepcopy(original)
        wrong_length["Sequences"][0]["Length"] = 99
        variants.append(wrong_length)
        outside = copy.deepcopy(original)
        outside["Sequences"][0]["Crisprs"][1]["End"] = 101
        variants.append(outside)
        duplicate = copy.deepcopy(original)
        duplicate["Sequences"].append(copy.deepcopy(duplicate["Sequences"][0]))
        variants.append(duplicate)
        invalid_level = copy.deepcopy(original)
        invalid_level["Sequences"][0]["Crisprs"][1]["Evidence_Level"] = 5
        variants.append(invalid_level)
        variants.append({"Sequences": []})
        for payload in variants:
            with self.subTest(payload=payload), self.assertRaises(AnnotationError):
                crispr_coordinates(payload, {"contig.2": 100})

    def test_matrix_counts_preserve_gene_units_and_missing_fields(self):
        rows = [
            {"features": '["K00001","K00002"]'},
            {"features": '["K00001"]'},
            {"features": "NA"},
        ]
        counts, missing = matrix_counts(rows, "features")
        self.assertEqual(counts, Counter(K00001=2, K00002=1))
        self.assertTrue(missing)
        self.assertEqual(
            matrix_counts([{"features": "[]"}], "features"), (Counter(), False)
        )
        with self.assertRaisesRegex(AnnotationError, "Duplicate"):
            matrix_counts([{"features": '["K00001","K00001"]'}], "features")

    def test_category_weights_conserve_one_unit_per_assigned_gene(self):
        rows = [
            {"cogclassifier_categories": '["R","S"]'},
            {"cogclassifier_categories": '["R","S","J"]'},
            {"cogclassifier_categories": "[]"},
        ]
        counts, weights, missing = category_counts(rows)
        self.assertFalse(missing)
        self.assertEqual(counts, Counter(R=2, S=2, J=1))
        self.assertLess(abs(sum(weights.values()) - 2), Decimal("1e-25"))
        self.assertLess(abs(weights["R"] - Decimal(5) / 6), Decimal("1e-25"))
        self.assertTrue(category_counts(rows + [{"cogclassifier_categories": "NA"}])[2])
        with self.assertRaisesRegex(AnnotationError, "Unknown"):
            category_counts([{"cogclassifier_categories": '["?"]'}])
