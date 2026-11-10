"""Genus thresholds, missing evidence and weighted quality-only representatives."""

import sys
import tempfile
import unittest
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from aai_common import aai_estimate, jaccard_cutoff
from cluster_aai import cluster_domain


def record(name, **overrides):
    row = dict(
        accession=name,
        domain="Bacteria",
        status="ready",
        reason="NA",
        marker_count=70,
        profile=name,
        quality=dict(
            path="proteins.faa",
            gcode="11",
            checkm2_completeness="95",
            checkm2_contamination="1",
            assembly_level="Complete Genome",
            n50="1000000",
            scaffolds="1",
            genome_size="1000000",
            busco="C:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200",
            organism_name=name,
        ),
    )
    row.update(overrides)
    return row


class AaiClusteringTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory(dir="/tmp")
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)

    def run_case(self, records, values, threshold=0.65):
        from aai_evidence import from_pairs

        path = self.root / "evidence"
        names = sorted(r["accession"] for r in records if r["status"] == "ready")
        from_pairs(path, names, values)
        return cluster_domain(records, path, threshold, "default", "Bacteria")

    def test_ani_scores_match_frozen_v048_for_every_profile(self):
        import json
        from ani_common import Genome
        from ani_representatives import select_representative_for_indices

        fixture = json.loads(
            (
                Path(__file__).resolve().parents[1]
                / "assets/testdata/fastaai/ani_v048_scores.json"
            ).read_text()
        )
        names = [row["Accession"] for row in fixture["metadata"]]
        for profile, expected in fixture["expected"].items():
            meta = {row["Accession"]: Genome(**row) for row in fixture["metadata"]}
            _, representative, _ = select_representative_for_indices(
                list(range(len(names))),
                names,
                meta,
                np.asarray(fixture["matrix"]),
                profile,
            )
            self.assertEqual(representative, expected["representative"])
            self.assertEqual(
                {a: f"{g.Score:.6f}" for a, g in meta.items()}, expected["scores"]
            )

    def test_calibration_and_bounds(self):
        for value in (0.31, 0.4, 0.65, 0.89):
            self.assertAlmostEqual(
                float(aai_estimate(jaccard_cutoff(value))), value * 100, places=10
            )
        for value in (float("nan"), float("inf"), 0.29, 0.91):
            with self.assertRaises(ValueError):
                jaccard_cutoff(value)
        self.assertEqual(aai_estimate(1), ">90%")
        self.assertEqual(aai_estimate(0), "<30%")

    def test_complete_linkage_does_not_join_a_chain(self):
        rows, reps, members, _ = self.run_case(
            [record("A"), record("B"), record("C")],
            {
                ("A", "B"): jaccard_cutoff(0.7),
                ("B", "C"): jaccard_cutoff(0.69),
                ("A", "C"): jaccard_cutoff(0.6),
            },
        )
        groups = {r["Accession"]: r["AAI_Cluster_ID"] for r in rows}
        self.assertEqual(groups["A"], groups["B"])
        self.assertNotEqual(groups["A"], groups["C"])
        self.assertEqual(len(reps), 2)

    def test_exact_cutoff_and_below(self):
        cutoff = jaccard_cutoff(0.65)
        for j, expected in ((cutoff, 1), (cutoff - 1e-8, 2), (cutoff + 1e-8, 1)):
            _, reps, _, _ = self.run_case([record("A"), record("B")], {("A", "B"): j})
            self.assertEqual(len(reps), expected)

    def test_calibrated_endpoints_are_inclusive(self):
        for threshold in (0.30, 0.90):
            _, reps, _, _ = self.run_case(
                [record("A"), record("B")],
                {("A", "B"): jaccard_cutoff(threshold)},
                threshold=threshold,
            )
            self.assertEqual(len(reps), 1)

    def test_missing_pair_prevents_a_complete_linkage_merge(self):
        rows, reps, _, _ = self.run_case(
            [record("A"), record("B"), record("C")],
            {("A", "B"): 0.8, ("B", "C"): 0.9},
        )
        self.assertEqual(len(reps), 2)
        groups = {r["Accession"]: r["AAI_Cluster_ID"] for r in rows}
        self.assertNotEqual(groups["A"], groups["C"])

    def test_domains_keep_distinct_cluster_ids(self):
        from aai_evidence import from_pairs

        identifiers = []
        for domain in ("Bacteria", "Archaea"):
            path = self.root / domain
            from_pairs(path, ["A", "B"], {("A", "B"): 1.0})
            rows, reps, _, _ = cluster_domain(
                [record("A", domain=domain), record("B", domain=domain)],
                path,
                0.65,
                "default",
                domain,
            )
            self.assertEqual(len(reps), 1)
            identifiers.append(rows[0]["AAI_Cluster_ID"])
        self.assertEqual(identifiers, ["AAI_B000001", "AAI_A000001"])

    def test_missing_pairs_do_not_invent_singleton_genera(self):
        rows, reps, _, _ = self.run_case([record("A"), record("B")], {})
        self.assertFalse(reps)
        self.assertEqual({r["aai_status"] for r in rows}, {"no_pairwise_evidence"})
        self.assertEqual({r["AAI_Cluster_ID"] for r in rows}, {"NA"})
        rows, reps, _, _ = self.run_case([record("A")], {})
        self.assertEqual(len(reps), 1)
        self.assertEqual(rows[0]["AAI_Is_Representative"], "yes")

    def test_marker_failures_exclusions_and_empty_cohort(self):
        rows, reps, _, _ = self.run_case(
            [
                record("A", status="no_markers", reason="no_usable_markers"),
                record("B", status="excluded", reason="low_quality"),
            ],
            {},
        )
        self.assertFalse(reps)
        self.assertEqual([r["aai_status"] for r in rows], ["no_markers", "excluded"])
        self.assertEqual(self.run_case([], {}), ([], [], [], []))

    def test_missing_metrics_exclude_only_representatives(self):
        a, b = record("A"), record("B")
        a["quality"]["busco"] = "NA"
        rows, reps, _, _ = self.run_case([a, b], {("A", "B"): 1})
        self.assertEqual(reps[0]["Representative_Accession"], "B")
        self.assertEqual(rows[0]["aai_status"], "clustered")
        self.assertEqual(rows[0]["AAI_Score"], "NA")
        self.assertEqual(rows[0]["AAI_to_Representative"], ">90%")
        b["quality"]["n50"] = "NA"
        rows, reps, _, _ = self.run_case([a, b], {("A", "B"): 1})
        self.assertEqual(reps[0]["Representative_Status"], "no_eligible_representative")
        self.assertEqual({r["AAI_Is_Representative"] for r in rows}, {"NA"})

    def test_missing_genome_size_does_not_remove_a_scoreable_candidate(self):
        a = record("A")
        a["quality"]["genome_size"] = "NA"
        rows, reps, _, _ = self.run_case([a], {})
        self.assertEqual(reps[0]["Representative_Accession"], "A")
        self.assertNotEqual(rows[0]["AAI_Score"], "NA")

    def test_invalid_metrics_remove_candidates_without_removing_membership(self):
        for field, value in (
            ("n50", "-1"),
            ("scaffolds", "0"),
            ("checkm2_completeness", "NaN"),
            ("checkm2_contamination", "inf"),
            ("busco", "C:101.0%[S:101.0%,D:0.0%],F:0.0%,M:0.0%,n:200"),
        ):
            with self.subTest(field=field):
                a, b = record("A"), record("B")
                a["quality"][field] = value
                rows, reps, _, _ = self.run_case([a, b], {("A", "B"): 1})
                self.assertEqual(reps[0]["Representative_Accession"], "B")
                self.assertEqual(rows[0]["aai_status"], "clustered")
                self.assertEqual(rows[0]["AAI_Score"], "NA")
                self.assertIn("invalid_", rows[0]["aai_representative_status"])

    def test_quality_ties_and_input_order_are_deterministic(self):
        _, reps, _, trace = self.run_case([record("B"), record("A")], {("B", "A"): 1})
        self.assertEqual(reps[0]["Representative_Accession"], "A")
        self.assertTrue(any("C=0.0" in row["line"] for row in trace))


if __name__ == "__main__":
    unittest.main()
