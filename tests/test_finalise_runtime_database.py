"""Tests for runtime database finalisation helpers."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import finalise_runtime_database as finalise_module


class FinaliseRuntimeDatabaseTestCase(unittest.TestCase):
    """Lock the runtime database finaliser contract."""

    def test_reuse_of_missing_database_fails_explicitly(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args = finalise_module.parse_args(
                [
                    "--component",
                    "padloc",
                    "--destination",
                    str(root / "missing"),
                    "--report",
                    str(root / "report.tsv"),
                    "--mode",
                    "reuse",
                ]
            )
            with self.assertRaisesRegex(
                finalise_module.PrepareRuntimeDatabasesError, "valid prepared database"
            ):
                finalise_module.finalise_runtime_database(args)
            self.assertFalse((root / "missing").exists())

    def test_prepared_mode_alias_writes_one_new_marker(self) -> None:
        """Accept the native download prep-mode token and treat it as a fresh write."""
        with tempfile.TemporaryDirectory() as temp_dir:
            destination = Path(temp_dir) / "padloc"
            (destination / "hmm").mkdir(parents=True)
            (destination / "hmm" / "padlocdb.hmm").write_text(
                "stub\n", encoding="utf-8"
            )

            args = finalise_module.parse_args(
                [
                    "--component",
                    "padloc",
                    "--destination",
                    str(destination),
                    "--report",
                    str(Path(temp_dir) / "report.tsv"),
                    "--mode",
                    "prepared",
                    "--source",
                    "native_download",
                ]
            )

            record = finalise_module.finalise_runtime_database(args)

            self.assertEqual(record.status, "prepared")
            self.assertEqual(record.component, "padloc")
            self.assertTrue((destination / ".nf_myco_ready.json").exists())


if __name__ == "__main__":
    unittest.main()
