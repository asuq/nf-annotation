"""Exact ANI parsing, bounded parser storage and changes between read passes."""

from __future__ import annotations

import os
import sys
import tempfile
import tracemalloc
import unittest
from pathlib import Path
from unittest import mock

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import ani_common
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
            [],
            ["-0.0"],
            ["NA", "100"],
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
        np.testing.assert_array_equal(matrix.view(np.uint64), expected.view(np.uint64))

    def test_singleton(self):
        self.path.write_text("1\nfirst sample.fa\n")
        names, matrix, indices = load_matrix(self.path)
        self.assertEqual(names, ["first sample.fa"])
        self.assertEqual(indices, {"first sample.fa": 0})
        np.testing.assert_array_equal(matrix, [[100.0]])

    def test_invalid_matrix_is_rejected(self):
        for text in (
            "",
            "\n\t\n",
            "not-an-integer\n",
            "1.0\n",
            "0\n",
            "-1\n",
            "1000000000\nfirst\n",
            "3\nA\nB 95\n",
            "1\nA\nB 95\n",
            "2\nA\nB\n",
            "3\nA\nB 95\nC 96\n",
            "2\nA\nA 95\n",
            *[
                f"2\nA\nB {value}\n"
                for value in ("bad", "nan", "inf", "-inf", "101", "-0.1")
            ],
        ):
            with self.subTest(text=text):
                self.path.write_text(text)
                with self.assertRaises(AniInputError):
                    load_matrix(self.path)
        with self.assertRaisesRegex(AniInputError, "not found"):
            load_matrix(self.root / "missing")

    def test_impossible_dimensions_fail_before_allocation(self):
        self.path.write_text("1000000000\nfirst\n")
        with mock.patch.object(np, "full", side_effect=AssertionError("allocated")):
            with self.assertRaises(AniInputError):
                load_matrix(self.path)

    def test_mutation_during_single_pass_is_rejected(self):
        original = "3\nA\nB 96\nC NA 97\n"
        rows = ani_common._matrix_rows
        for replace in (False, True):
            self.path.write_text(original)

            def mutate(lines, count):
                for index, item in enumerate(rows(lines, count)):
                    if index == 0:
                        info = self.path.stat()
                        if replace:
                            other = self.path.with_suffix(".new")
                            other.write_text(original)
                            other.replace(self.path)
                        else:
                            self.path.write_text(original.replace("96", "95"))
                            os.utime(self.path, ns=(info.st_atime_ns, info.st_mtime_ns))
                    yield item

            with mock.patch.object(ani_common, "_matrix_rows", side_effect=mutate):
                with self.assertRaisesRegex(AniInputError, "changed while reading"):
                    load_matrix(self.path)

    def test_threshold_adjacent_values_keep_complete_linkage_membership(self):
        self.path.write_text("4\nB\nA 95.0000000001\nC 94.9999999999 95\nD NA NA NA\n")
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
