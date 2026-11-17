import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from assembly_eligibility import classify


class AssemblyEligibilityTests(unittest.TestCase):
    def test_original_order_and_strict_boundary(self):
        samples = [{"accession": name} for name in ("B", "A", "C")]
        rows = classify(
            samples,
            [
                {"accession": name, "genome_size": length}
                for name, length in (("A", "20000"), ("C", "8471"), ("B", "19999"))
            ],
        )
        self.assertEqual([r["accession"] for r in rows], ["B", "A", "C"])
        self.assertEqual(
            [r["Assembly_Exclusion_Reason"] for r in rows],
            ["assembly_below_20000_bp", "NA", "assembly_below_20000_bp"],
        )
