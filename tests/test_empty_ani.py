"""Distinguish a deliberately empty ANI cohort from missing or corrupt data."""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import ani_common  # noqa: E402
import build_fastani_inputs  # noqa: E402
import cluster_ani  # noqa: E402
import select_ani_representatives  # noqa: E402


class EmptyAniTestCase(unittest.TestCase):
    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory(dir="/tmp", prefix="empty-ani-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.metadata = self.root / "metadata.tsv"
        self.metadata.write_text(
            "\t".join(
                [*build_fastani_inputs.ANI_METADATA_COLUMNS, "BUSCO_bacillota_odb12"]
            )
            + "\n"
        )
        self.matrix = self.root / "matrix"
        self.matrix.touch()
        self.clusters = self.root / "cluster.tsv"
        self.clusters.write_text("Accession\tCluster_ID\tMatrix_Name\n")

    def test_empty_cohort_emits_header_only_cluster_and_representative_tables(
        self,
    ) -> None:
        self.assertEqual(
            cluster_ani.main(
                [
                    "--ani-metadata",
                    str(self.metadata),
                    "--ani-matrix",
                    str(self.matrix),
                    "--outdir",
                    str(self.root),
                ]
            ),
            0,
        )
        summary = self.root / "ani_summary.tsv"
        representatives = self.root / "ani_representatives.tsv"
        self.assertEqual(
            select_ani_representatives.main(
                [
                    "--ani-clusters",
                    str(self.clusters),
                    "--ani-metadata",
                    str(self.metadata),
                    "--ani-matrix",
                    str(self.matrix),
                    "--ani-summary-output",
                    str(summary),
                    "--ani-representatives-output",
                    str(representatives),
                ]
            ),
            0,
        )
        for path in (self.clusters, summary, representatives):
            with path.open() as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                self.assertTrue(reader.fieldnames)
                self.assertEqual(list(reader), [])

    def test_missing_matrix_is_not_an_empty_cohort(self) -> None:
        self.matrix.unlink()
        with self.assertRaises(ani_common.AniInputError):
            ani_common.is_empty_ani_cohort(self.metadata, self.matrix)

    def test_matrix_data_without_metadata_is_rejected(self) -> None:
        self.matrix.write_text("1\nA.fasta\n")
        with self.assertRaisesRegex(ani_common.AniInputError, "non-empty"):
            ani_common.is_empty_ani_cohort(self.metadata, self.matrix)

    def test_metadata_without_matrix_is_not_converted_to_an_empty_cohort(self) -> None:
        self.metadata.write_text(self.metadata.read_text() + "A\tA.fasta\n")
        self.assertFalse(ani_common.is_empty_ani_cohort(self.metadata, self.matrix))
        self.assertNotEqual(
            cluster_ani.main(
                [
                    "--ani-metadata",
                    str(self.metadata),
                    "--ani-matrix",
                    str(self.matrix),
                    "--outdir",
                    str(self.root),
                ]
            ),
            0,
        )

    def test_missing_or_duplicate_metadata_headers_are_rejected(self) -> None:
        for header in ("", "accession\n", "accession\tmatrix_name\taccession\n"):
            self.metadata.write_text(header)
            with self.assertRaises(ani_common.AniInputError):
                ani_common.is_empty_ani_cohort(self.metadata, self.matrix)

    def test_empty_scoring_cohort_still_requires_the_scoring_schema(self) -> None:
        self.metadata.write_text("accession\tmatrix_name\n")
        with self.assertRaises(ani_common.AniInputError):
            ani_common.is_empty_ani_cohort(
                self.metadata, self.matrix, require_scoring=True
            )

    def test_cluster_rows_without_eligible_metadata_are_rejected(self) -> None:
        self.clusters.write_text(self.clusters.read_text() + "A\tC000001\tA.fasta\n")
        with self.assertRaisesRegex(
            select_ani_representatives.RepresentativeSelectionError,
            "clusters contain samples",
        ):
            select_ani_representatives.build_ani_outputs(
                ani_clusters=self.clusters,
                ani_metadata=self.metadata,
                ani_matrix=self.matrix,
                score_profile="default",
            )


if __name__ == "__main__":
    unittest.main()
