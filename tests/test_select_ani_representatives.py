"""Tests for ANI representative selection from raw cluster and matrix inputs."""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import select_ani_representatives  # noqa: E402


def read_tsv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV file into a header and row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames is not None
        return reader.fieldnames, list(reader)


class SelectAniRepresentativesTestCase(unittest.TestCase):
    """Cover ANI representative selection and output contracts."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write a UTF-8 text file and return its path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def test_main_writes_stable_summary_and_representative_tables(self) -> None:
        """Emit stable ANI summary and representative rows from valid raw inputs."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            ani_clusters = self.write_text_file(
                tmpdir / "cluster.tsv",
                "\n".join(
                    [
                        "Accession\tCluster_ID\tMatrix_Name",
                        "ACC1\tC000001\tfastani_inputs/ACC1.fasta",
                        "ACC2\tC000001\tfastani_inputs/ACC2.fasta",
                    ]
                )
                + "\n",
            )
            ani_metadata = self.write_text_file(
                tmpdir / "ani_metadata.tsv",
                "\n".join(
                    [
                        "accession\tmatrix_name\tpath\tassembly_level\tgcode\tcheckm2_completeness\tcheckm2_contamination\tn50\tscaffolds\tgenome_size\torganism_name\tBUSCO_bacillota_odb12",
                        "ACC1\tfastani_inputs/ACC1.fasta\tfastani_inputs/ACC1.fasta\tComplete Genome\t11\t97\t1\t100000\t1\t900000\tOne\tC:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200",
                        "ACC2\tfastani_inputs/ACC2.fasta\tfastani_inputs/ACC2.fasta\tScaffold\t11\t93\t2\t50000\t5\t850000\tTwo\tC:96.0%[S:96.0%,D:0.0%],F:2.0%,M:2.0%,n:200",
                    ]
                )
                + "\n",
            )
            ani_matrix = self.write_text_file(
                tmpdir / "fastani.matrix",
                "\n".join(
                    [
                        "2",
                        "fastani_inputs/ACC1.fasta",
                        "fastani_inputs/ACC2.fasta 97.2500",
                    ]
                )
                + "\n",
            )
            ani_summary = tmpdir / "ani_summary.tsv"
            ani_representatives = tmpdir / "ani_representatives.tsv"

            exit_code = select_ani_representatives.main(
                [
                    "--ani-clusters",
                    str(ani_clusters),
                    "--ani-metadata",
                    str(ani_metadata),
                    "--ani-matrix",
                    str(ani_matrix),
                    "--ani-summary-output",
                    str(ani_summary),
                    "--ani-representatives-output",
                    str(ani_representatives),
                ]
            )

            self.assertEqual(exit_code, 0)

            summary_header, summary_rows = read_tsv_rows(ani_summary)
            representative_header, representative_rows = read_tsv_rows(
                ani_representatives
            )

            self.assertEqual(
                summary_header,
                [
                    "Accession",
                    "Cluster_ID",
                    "Is_Representative",
                    "ANI_to_Representative",
                    "Score",
                ],
            )
            self.assertEqual(
                representative_header,
                [
                    "Cluster_ID",
                    "Representative_Accession",
                    "Organism_Name",
                    "CheckM2_Completeness",
                    "CheckM2_Contamination",
                    "BUSCO",
                    "Assembly_Level",
                    "N50",
                    "Cluster_Size",
                ],
            )

            summary_by_accession = {row["Accession"]: row for row in summary_rows}
            self.assertEqual(summary_by_accession["ACC1"]["Cluster_ID"], "C000001")
            self.assertEqual(summary_by_accession["ACC1"]["Is_Representative"], "yes")
            self.assertEqual(
                summary_by_accession["ACC1"]["ANI_to_Representative"], "100.0000"
            )
            self.assertNotEqual(summary_by_accession["ACC1"]["Score"], "NA")
            self.assertEqual(summary_by_accession["ACC2"]["Is_Representative"], "no")
            self.assertEqual(
                summary_by_accession["ACC2"]["ANI_to_Representative"], "97.2500"
            )

            self.assertEqual(len(representative_rows), 1)
            self.assertEqual(representative_rows[0]["Cluster_ID"], "C000001")
            self.assertEqual(representative_rows[0]["Representative_Accession"], "ACC1")

    def test_main_breaks_full_ties_by_lexicographic_accession(self) -> None:
        """Pick the lexicographically smallest accession after a full tie cascade."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            ani_clusters = self.write_text_file(
                tmpdir / "cluster.tsv",
                "\n".join(
                    [
                        "Accession\tCluster_ID",
                        "ACC1\tC000001",
                        "ACC2\tC000001",
                    ]
                )
                + "\n",
            )
            ani_metadata = self.write_text_file(
                tmpdir / "ani_metadata.tsv",
                "\n".join(
                    [
                        "accession\tmatrix_name\tpath\tassembly_level\tgcode\tcheckm2_completeness\tcheckm2_contamination\tn50\tscaffolds\tgenome_size\torganism_name\tBUSCO_bacillota_odb12",
                        "ACC2\tfastani_inputs/ACC2.fasta\tfastani_inputs/ACC2.fasta\tScaffold\t11\t95\t1\t50000\t4\t900000\tTwo\tC:97.0%[S:97.0%,D:0.0%],F:2.0%,M:1.0%,n:200",
                        "ACC1\tfastani_inputs/ACC1.fasta\tfastani_inputs/ACC1.fasta\tScaffold\t11\t95\t1\t50000\t4\t900000\tOne\tC:97.0%[S:97.0%,D:0.0%],F:2.0%,M:1.0%,n:200",
                    ]
                )
                + "\n",
            )
            ani_matrix = self.write_text_file(
                tmpdir / "fastani.matrix",
                "\n".join(
                    [
                        "2",
                        "fastani_inputs/ACC1.fasta",
                        "fastani_inputs/ACC2.fasta 99.9000",
                    ]
                )
                + "\n",
            )
            ani_summary = tmpdir / "ani_summary.tsv"
            ani_representatives = tmpdir / "ani_representatives.tsv"

            exit_code = select_ani_representatives.main(
                [
                    "--ani-clusters",
                    str(ani_clusters),
                    "--ani-metadata",
                    str(ani_metadata),
                    "--ani-matrix",
                    str(ani_matrix),
                    "--ani-summary-output",
                    str(ani_summary),
                    "--ani-representatives-output",
                    str(ani_representatives),
                ]
            )

            self.assertEqual(exit_code, 0)
            _, representative_rows = read_tsv_rows(ani_representatives)
            self.assertEqual(representative_rows[0]["Representative_Accession"], "ACC1")

    def test_main_fails_when_cluster_accession_is_missing_from_metadata(self) -> None:
        """Fail when cluster.tsv references an accession absent from ANI metadata."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            ani_clusters = self.write_text_file(
                tmpdir / "cluster.tsv",
                "Accession\tCluster_ID\nACC1\tC000001\nACC2\tC000001\n",
            )
            ani_metadata = self.write_text_file(
                tmpdir / "ani_metadata.tsv",
                "\n".join(
                    [
                        "accession\tmatrix_name\tpath\tassembly_level\tgcode\tcheckm2_completeness\tcheckm2_contamination\tn50\tscaffolds\tgenome_size\torganism_name\tBUSCO_bacillota_odb12",
                        "ACC1\tfastani_inputs/ACC1.fasta\tfastani_inputs/ACC1.fasta\tScaffold\t11\t95\t1\t50000\t4\t900000\tOne\tC:97.0%[S:97.0%,D:0.0%],F:2.0%,M:1.0%,n:200",
                    ]
                )
                + "\n",
            )
            ani_matrix = self.write_text_file(
                tmpdir / "fastani.matrix",
                "\n".join(
                    [
                        "1",
                        "fastani_inputs/ACC1.fasta",
                    ]
                )
                + "\n",
            )

            exit_code = select_ani_representatives.main(
                [
                    "--ani-clusters",
                    str(ani_clusters),
                    "--ani-metadata",
                    str(ani_metadata),
                    "--ani-matrix",
                    str(ani_matrix),
                    "--ani-summary-output",
                    str(tmpdir / "ani_summary.tsv"),
                    "--ani-representatives-output",
                    str(tmpdir / "ani_representatives.tsv"),
                ]
            )

            self.assertEqual(exit_code, 1)

    def test_main_fails_when_matrix_name_is_missing_from_metadata(self) -> None:
        """Fail when the FastANI matrix references names missing from ANI metadata."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            ani_clusters = self.write_text_file(
                tmpdir / "cluster.tsv",
                "Accession\tCluster_ID\nACC1\tC000001\n",
            )
            ani_metadata = self.write_text_file(
                tmpdir / "ani_metadata.tsv",
                "\n".join(
                    [
                        "accession\tmatrix_name\tpath\tassembly_level\tgcode\tcheckm2_completeness\tcheckm2_contamination\tn50\tscaffolds\tgenome_size\torganism_name\tBUSCO_bacillota_odb12",
                        "ACC1\tfastani_inputs/ACC1.fasta\tfastani_inputs/ACC1.fasta\tScaffold\t11\t95\t1\t50000\t4\t900000\tOne\tC:97.0%[S:97.0%,D:0.0%],F:2.0%,M:1.0%,n:200",
                    ]
                )
                + "\n",
            )
            ani_matrix = self.write_text_file(
                tmpdir / "fastani.matrix",
                "\n".join(
                    [
                        "2",
                        "fastani_inputs/ACC1.fasta",
                        "fastani_inputs/ACC2.fasta 97.0000",
                    ]
                )
                + "\n",
            )

            exit_code = select_ani_representatives.main(
                [
                    "--ani-clusters",
                    str(ani_clusters),
                    "--ani-metadata",
                    str(ani_metadata),
                    "--ani-matrix",
                    str(ani_matrix),
                    "--ani-summary-output",
                    str(tmpdir / "ani_summary.tsv"),
                    "--ani-representatives-output",
                    str(tmpdir / "ani_representatives.tsv"),
                ]
            )

            self.assertEqual(exit_code, 1)


if __name__ == "__main__":
    unittest.main()
