"""Prediction completion, valid zero and missing evidence stay distinct."""

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))

from annotation_common import AnnotationError
from prokka_prediction import read_prediction, status_from_manifest


class ProkkaPredictionTests(unittest.TestCase):
    def test_successful_records_include_zero_and_both_modes(self):
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / "prokka.log"
            for mode, count in (("single", 23), ("meta", 0)):
                log.write_text(
                    f"NF_PROKKA_PREDICTION_V1\t{mode}\t{count}\nexit_code=0\n"
                )
                prediction = read_prediction(log)
                self.assertEqual((prediction.mode, prediction.cds_count), (mode, count))

    def test_failed_missing_partial_and_ambiguous_records_never_become_zero(self):
        cases = (
            "exit_code=0\n",
            "Found 0 CDS\nError: training failed\nexit_code=0\n",
            "NF_PROKKA_PREDICTION_V1\tmeta\t0\nexit_code=7\n",
            "NF_PROKKA_PREDICTION_V1\tmeta\t0\n",
            "NF_PROKKA_PREDICTION_V1\tunknown\t0\nexit_code=0\n",
            "NF_PROKKA_PREDICTION_V1\tmeta\tNA\nexit_code=0\n",
            "NF_PROKKA_PREDICTION_V1\tmeta\t-1\nexit_code=0\n",
            "NF_PROKKA_PREDICTION_V1\tmeta\t0\n" * 2 + "exit_code=0\n",
        )
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / "prokka.log"
            for text in cases:
                with self.subTest(text=text):
                    log.write_text(text)
                    with self.assertRaises(AnnotationError):
                        read_prediction(log)

    def test_only_final_attempt_defines_the_prediction(self):
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / "prokka.log"
            log.write_text(
                "attempt=1/2\nNF_PROKKA_PREDICTION_V1\tsingle\t7\nlater tool failed\nattempt=2/2\nNF_PROKKA_PREDICTION_V1\tmeta\t0\nexit_code=0\n"
            )
            self.assertEqual(read_prediction(log).cds_count, 0)

    def test_status_requires_verified_prediction_and_consistent_outputs(self):
        base = {
            "exit_code": "0",
            "gff_size": "100",
            "faa_size": "50",
            "prodigal_mode": "meta",
            "cds_count": "2",
        }
        for change, expected in (
            ({}, "done"),
            ({"faa_size": "0", "cds_count": "0"}, "no_cds"),
            ({"exit_code": "7", "faa_size": "0", "cds_count": "0"}, "failed"),
            ({"faa_size": "0"}, "failed"),
            ({"gff_size": "0", "faa_size": "0", "cds_count": "0"}, "failed"),
            ({"cds_count": "0"}, "failed"),
            ({"prodigal_mode": "NA"}, "failed"),
        ):
            with self.subTest(change=change):
                status, warnings = status_from_manifest(base | change)
                self.assertEqual(status, expected)
                self.assertEqual(bool(warnings), expected == "failed")
