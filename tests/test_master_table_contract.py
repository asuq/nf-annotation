"""Tests for the master-table output contract."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import master_table_contract  # noqa: E402


class MasterTableContractTestCase(unittest.TestCase):
    """Verify the locked derived-column order and asset consistency."""

    def test_default_contract_matches_repository_asset(self) -> None:
        """Keep the asset file aligned with the default code contract."""
        master_table_contract.validate_default_append_columns_asset()
        self.assertEqual(
            master_table_contract.read_append_columns_asset(),
            master_table_contract.build_append_columns(),
        )

    def test_default_sample_status_contract_matches_repository_asset(self) -> None:
        """Keep the sample-status asset aligned with the default code contract."""
        master_table_contract.validate_default_sample_status_columns_asset()
        self.assertEqual(
            master_table_contract.read_sample_status_columns_asset(),
            master_table_contract.build_sample_status_columns(),
        )

    def test_build_master_table_columns_supports_custom_busco_order(self) -> None:
        """Insert BUSCO columns in caller-provided lineage order."""
        columns = master_table_contract.build_master_table_columns(
            metadata_columns=["Accession", "Tax_ID"],
            busco_lineages=["lineage_a", "lineage_b", "lineage_c"],
        )

        self.assertEqual(columns[:2], ["Accession", "Tax_ID"])
        self.assertIn("BUSCO_lineage_a", columns)
        self.assertIn("BUSCO_lineage_b", columns)
        self.assertIn("BUSCO_lineage_c", columns)
        self.assertLess(columns.index("BUSCO_lineage_a"), columns.index("BUSCO_lineage_b"))
        self.assertLess(columns.index("BUSCO_lineage_b"), columns.index("BUSCO_lineage_c"))
        self.assertEqual(columns[-len(master_table_contract.ANNOTATION_COLUMNS):], list(master_table_contract.ANNOTATION_COLUMNS))
        start = columns.index("Cluster_ID")
        self.assertEqual(columns[start:start + 4], list(master_table_contract.ANI_COLUMNS))

    def test_build_sample_status_columns_supports_custom_busco_order(self) -> None:
        """Insert BUSCO status columns in caller-provided lineage order."""
        columns = master_table_contract.build_sample_status_columns(
            ["lineage_a", "lineage_b", "lineage_c"]
        )

        self.assertEqual(
            columns[: len(master_table_contract.SAMPLE_STATUS_PREFIX_COLUMNS)],
            list(master_table_contract.SAMPLE_STATUS_PREFIX_COLUMNS),
        )
        self.assertIn("busco_lineage_a_status", columns)
        self.assertIn("busco_lineage_b_status", columns)
        self.assertIn("busco_lineage_c_status", columns)
        self.assertLess(
            columns.index("busco_lineage_a_status"),
            columns.index("busco_lineage_b_status"),
        )
        self.assertLess(
            columns.index("busco_lineage_b_status"),
            columns.index("busco_lineage_c_status"),
        )
        self.assertEqual(
            columns[-len(master_table_contract.SAMPLE_STATUS_SUFFIX_COLUMNS) :],
            list(master_table_contract.SAMPLE_STATUS_SUFFIX_COLUMNS),
        )

    def test_codetta_columns_follow_gcode_immediately(self) -> None:
        """Keep the Codetta fields directly after the chosen gcode column."""
        columns = master_table_contract.build_append_columns()

        provenance_start = columns.index("Gcode") + 1
        self.assertEqual(
            columns[provenance_start:provenance_start + len(master_table_contract.GCODE_PROVENANCE_COLUMNS)],
            list(master_table_contract.GCODE_PROVENANCE_COLUMNS),
        )
        self.assertEqual(
            columns[columns.index("Codetta_Genetic_Code") + 1],
            "Codetta_NCBI_Table_Candidates",
        )

    def test_gc_content_precedes_gcode(self) -> None:
        """Keep the seqtk-derived GC column before the gcode block."""
        columns = master_table_contract.build_append_columns()

        self.assertEqual(columns[0], "is_new")
        self.assertEqual(
            columns[columns.index("Total_Coding_Sequences_gcode11") + 1],
            "GC_Content",
        )
        self.assertEqual(columns[columns.index("GC_Content") + 1], "Gcode")

    def test_normalise_busco_lineages_rejects_duplicates(self) -> None:
        """Fail when BUSCO lineage names repeat after normalisation."""
        with self.assertRaises(ValueError):
            master_table_contract.normalise_busco_lineages(
                ["bacillota_odb12", "bacillota_odb12"]
            )

    def test_build_master_table_columns_rejects_metadata_overlap(self) -> None:
        """Fail when metadata columns collide with derived output columns."""
        with self.assertRaises(ValueError):
            master_table_contract.build_master_table_columns(["Accession", "Gcode"])

    def test_extract_busco_lineages_from_runtime_contracts(self) -> None:
        """Recover BUSCO lineage order from generated reporting contracts."""
        append_lineages = master_table_contract.extract_busco_lineages_from_append_columns(
            master_table_contract.build_append_columns(["lineage_b", "lineage_a"])
        )
        status_lineages = master_table_contract.extract_busco_lineages_from_sample_status_columns(
            master_table_contract.build_sample_status_columns(["lineage_b", "lineage_a"])
        )

        self.assertEqual(append_lineages, ("lineage_b", "lineage_a"))
        self.assertEqual(status_lineages, ("lineage_b", "lineage_a"))

    def test_resolve_contracts_support_direct_lineage_inputs(self) -> None:
        """Build both output contracts directly from configured lineage names."""
        append_columns, append_lineages = master_table_contract.resolve_append_columns(
            busco_lineages=["lineage_b", "lineage_a"]
        )
        status_columns, status_lineages = master_table_contract.resolve_sample_status_columns(
            busco_lineages=["lineage_b", "lineage_a"]
        )

        self.assertEqual(append_lineages, ("lineage_b", "lineage_a"))
        self.assertEqual(status_lineages, ("lineage_b", "lineage_a"))
        self.assertIn("BUSCO_lineage_b", append_columns)
        self.assertIn("busco_lineage_b_status", status_columns)


if __name__ == "__main__":
    unittest.main()
