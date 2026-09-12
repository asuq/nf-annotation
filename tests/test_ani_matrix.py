"""Exact ANI parsing, bounded parser storage and changes between read passes."""

from __future__ import annotations

import io
import sys
import tempfile
import tracemalloc
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

from ani_common import AniInputError, load_matrix
from cluster_ani import (
    assign_cluster_ids,
    build_cluster_rows,
    cluster_complete_linkage,
)


class AniMatrixTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="ani-matrix-", dir="/tmp")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.path = self.root / "matrix"

    def test_exact_float64_values_missingness_names_and_order(self):
        names = ["Z with spaces.fa", "A'quote.fa", "tab\tinside.fa", "B\\slash.fa"]
        tokens = [
            [], ["-0.0"], ["NA", "100"],
            ["9.5e1", "94.99999999999", "95.00000000001"],
        ]
        self.path.write_text(
            "\r\n 4\r\n\r\n"
            + "\r\n\r\n".join(
                "  " + name + "\t" + "\t".join(values) + "  "
                for name, values in zip(names, tokens, strict=True)
            ),
            newline="",
        )
        expected = np.full((4, 4), np.nan, dtype=np.float64)
        np.fill_diagonal(expected, 100.0)
        for i, values in enumerate(tokens):
            for j, value in enumerate(values):
                if value != "NA":
                    expected[i, j] = expected[j, i] = float(value)
        observed_names, matrix, indices = load_matrix(self.path)
        self.assertEqual(observed_names, names)
        self.assertEqual(indices, {name: i for i, name in enumerate(names)})
        self.assertEqual(matrix.dtype, np.float64)
        np.testing.assert_array_equal(
            matrix.view(np.uint64), expected.view(np.uint64)
        )

    def test_singleton(self):
        self.path.write_text("1\nfirst sample.fa\n")
        names, matrix, indices = load_matrix(self.path)
        self.assertEqual(names, ["first sample.fa"])
        self.assertEqual(indices, {"first sample.fa": 0})
        np.testing.assert_array_equal(matrix, [[100.0]])

    def test_invalid_matrix_fails_before_dense_allocation(self):
        for text in (
            "", "\n\t\n", "not-an-integer\n", "1.0\n", "0\n", "-1\n",
            "1000000000\nfirst\n", "3\nA\nB 95\n", "1\nA\nB 95\n",
            "2\nA\nB\n", "3\nA\nB 95\nC 96\n", "2\nA\nA 95\n",
            *[
                f"2\nA\nB {value}\n"
                for value in ("bad", "nan", "inf", "-inf", "101", "-0.1")
            ],
        ):
            with self.subTest(text=text):
                self.path.write_text(text)
                with mock.patch.object(
                    np, "full", side_effect=AssertionError("allocated")
                ):
                    with self.assertRaises(AniInputError):
                        load_matrix(self.path)
        with self.assertRaisesRegex(AniInputError, "not found"):
            load_matrix(self.root / "missing")

    def test_second_pass_rejects_changed_or_truncated_input(self):
        original = "3\nA\nB 96\nC NA 97\n"

        class ChangingInput(io.StringIO):
            def __init__(self, replacement):
                super().__init__(original)
                self.replacement = replacement

            def seek(self, offset, whence=0):
                position = super().seek(offset, whence)
                if offset == 0 and whence == 0:
                    self.truncate()
                    self.write(self.replacement)
                    super().seek(0)
                return position

        self.path.write_text(original)
        for replacement in (
            "2\nA\nB 96\n",                  # changed declared count
            "3\nB\nA 96\nC NA 97\n",       # changed name order
            "3\nA\nB 96\n",                # truncated rows
            "3\nA\nB 96\nC 97\n",          # truncated row width
            "3\nA\nB 96\nC NA inf\n",      # invalid numeric value
            "3\nA\nB 96\nC NA 97\nD\n",   # extra row
            "3\nA\nB 95\nC NA 97\n",       # valid but changed value
            original.replace("\n", "\r\n"),  # changed source bytes
        ):
            with self.subTest(replacement=replacement):
                with mock.patch.object(
                    Path, "open", return_value=ChangingInput(replacement)
                ):
                    with self.assertRaises(AniInputError):
                        load_matrix(self.path)

    def test_threshold_adjacent_values_keep_complete_linkage_membership(self):
        self.path.write_text(
            "4\nB\nA 95.0000000001\nC 94.9999999999 95\nD NA NA NA\n"
        )
        names, matrix, _ = load_matrix(self.path)
        accessions = {name: name for name in names}
        clusters = cluster_complete_linkage(matrix, names, 0.95)
        stable = assign_cluster_ids(clusters, names, accessions)
        self.assertEqual(
            build_cluster_rows(stable, names, accessions),
            [
                {"Accession": "A", "Cluster_ID": "C000001", "Matrix_Name": "A"},
                {"Accession": "B", "Cluster_ID": "C000001", "Matrix_Name": "B"},
                {"Accession": "C", "Cluster_ID": "C000002", "Matrix_Name": "C"},
                {"Accession": "D", "Cluster_ID": "C000003", "Matrix_Name": "D"},
            ],
        )

    def test_parser_overhead_does_not_retain_quadratic_token_lists(self):
        count = 600
        with self.path.open("w") as handle:
            handle.write(f"{count}\n")
            for i in range(count):
                handle.write(f"sample_{i}" + "\t97.2500" * i + "\n")
        tracemalloc.start()
        try:
            _, matrix, _ = load_matrix(self.path)
            _, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        # Allow generous linear parser/index overhead above the required array.
        self.assertLess(peak, matrix.nbytes + 2_000_000)
        self.assertEqual(matrix[599, 0], 97.25)


if __name__ == "__main__":
    unittest.main()
