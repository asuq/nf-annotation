"""Helper-level tests for ANI metadata parsing."""

from __future__ import annotations

import contextlib
import csv
import importlib.util
import io
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import ani_common
import cluster_ani  # noqa: E402

SCIPY_AVAILABLE = importlib.util.find_spec("scipy") is not None


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a TSV file into row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class ClusterAniHelperTestCase(unittest.TestCase):
    """Cover header normalization, BUSCO parsing, and ANI metadata gating."""

    def test_normalize_header_matches_nextflow_style(self) -> None:
        """Normalize mixed punctuation and case into lowercase underscore tokens."""
        self.assertEqual(ani_common.normalize_header(" Assembly Level "), "assembly_level")
        self.assertEqual(ani_common.normalize_header("BUSCO-bacillota_odb12"), "busco_bacillota_odb12")

    def test_try_parse_busco_extracts_complete_and_missing(self) -> None:
        """Read the compact BUSCO string format emitted by the BUSCO summariser."""
        parsed = ani_common.try_parse_busco(
            "C:96.5%[S:95.0%,D:1.5%],F:2.0%,M:1.5%,n:220"
        )
        self.assertEqual(parsed, (96.5, 1.5))
        self.assertIsNone(ani_common.try_parse_busco("NA"))

    def test_build_genome_from_row_accepts_valid_metadata(self) -> None:
        """Produce a Genome record when all ANI-scoring metadata are present."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            genome_path = tmpdir / "genome.fna"
            genome_path.write_text(">a\nACGT\n", encoding="utf-8")
            genome, reasons = ani_common.build_genome_from_row(
                {
                    "path": str(genome_path),
                    "assembly_level": "Scaffolds",
                    "gcode": "11",
                    "checkm2_completeness": "95.0",
                    "checkm2_contamination": "1.5",
                    "n50": "50000",
                    "scaffolds": "4",
                    "genome_size": "1000000",
                    "primary_busco": "C:96.5%[S:95.0%,D:1.5%],F:2.0%,M:1.5%,n:220",
                    "organism_name": "Example",
                },
                acc="ACC1",
                busco_column="primary_busco",
            )

            self.assertEqual(reasons, [])
            assert genome is not None
            self.assertEqual(genome.Accession, "ACC1")
            self.assertEqual(genome.Gcode, 11)
            self.assertEqual(genome.Assembly_Level, "Scaffold")
            self.assertEqual(genome.BUSCO_C, 96.5)
            self.assertEqual(genome.BUSCO_M, 1.5)

    def test_build_genome_from_row_accumulates_exclusion_reasons(self) -> None:
        """Collect all missing or invalid ANI-scoring fields before exclusion."""
        genome, reasons = ani_common.build_genome_from_row(
            {
                "path": "missing.fna",
                "assembly_level": "unknown",
                "gcode": "NA",
                "checkm2_completeness": "NA",
                "checkm2_contamination": "bad",
                "n50": "1.5",
                "scaffolds": "",
                "genome_size": "NA",
                "primary_busco": "NA",
            },
            acc="ACC2",
            busco_column="primary_busco",
        )

        self.assertIsNone(genome)
        self.assertEqual(
            reasons,
            [
                "path_not_found",
                "invalid_assembly_level",
                "invalid_gcode",
                "missing_checkm2_completeness",
                "missing_checkm2_contamination",
                "missing_n50",
                "missing_scaffolds",
                "missing_genome_size",
                "missing_primary_busco",
            ],
        )

    def test_load_ani_metadata_supports_distinct_matrix_name_column(self) -> None:
        """Match PHYLIP matrix rows via a dedicated matrix-name column."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            genome_path = tmpdir / "genome.fna"
            genome_path.write_text(">a\nACGT\n", encoding="utf-8")
            metadata_tsv = tmpdir / "ani_metadata.tsv"
            metadata_tsv.write_text(
                "\n".join(
                    [
                        "accession\tmatrix_name\tpath\tassembly_level\tgcode\tcheckm2_completeness\tcheckm2_contamination\tn50\tscaffolds\tgenome_size\tBUSCO_bacillota_odb12",
                        f"ACC1\tfastani_inputs/ACC1.fasta\t{genome_path}\tScaffold\t11\t95\t1.5\t50000\t4\t1000000\tC:96.5%[S:95.0%,D:1.5%],F:2.0%,M:1.5%,n:220",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            metadata, eligible_names = ani_common.load_ani_metadata(
                ani_metadata=metadata_tsv,
                matrix_names=["fastani_inputs/ACC1.fasta"],
                busco_column="BUSCO_bacillota_odb12",
                matrix_name_column="matrix_name",
            )

            self.assertEqual(eligible_names, ["fastani_inputs/ACC1.fasta"])
            self.assertIn("fastani_inputs/ACC1.fasta", metadata)
            self.assertEqual(metadata["fastani_inputs/ACC1.fasta"].Accession, "ACC1")

    def test_main_fails_clearly_when_scipy_is_unavailable(self) -> None:
        """Report a clear runtime contract failure when SciPy cannot be imported."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            matrix = tmpdir / "fastani.matrix"
            matrix.write_text(
                "\n".join(
                    [
                        "2",
                        "fastani_inputs/ACC1.fasta",
                        "fastani_inputs/ACC2.fasta 97.5000",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            metadata_tsv = tmpdir / "ani_metadata.tsv"
            metadata_tsv.write_text(
                "\n".join(
                    [
                        "accession\tmatrix_name\tpath\tassembly_level\tgcode\tcheckm2_completeness\tcheckm2_contamination\tn50\tscaffolds\tgenome_size\torganism_name\tBUSCO_bacillota_odb12",
                        "ACC1\tfastani_inputs/ACC1.fasta\tfastani_inputs/ACC1.fasta\tScaffold\t11\t94\t1.5\t49000\t5\t1000000\tOne\tC:97.0%[S:97.0%,D:0.0%],F:2.0%,M:1.0%,n:200",
                        "ACC2\tfastani_inputs/ACC2.fasta\tfastani_inputs/ACC2.fasta\tScaffold\t11\t95\t1.0\t50000\t4\t1000000\tTwo\tC:96.0%[S:96.0%,D:0.0%],F:2.0%,M:2.0%,n:200",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            outdir = tmpdir / "out"
            stderr = io.StringIO()

            with (
                contextlib.redirect_stderr(stderr),
                mock.patch.object(
                    cluster_ani,
                    "load_scipy_complete_linkage",
                    side_effect=cluster_ani.AniInputError(
                        "SciPy is required for ANI clustering. Build "
                        "params.python_container with scipy and numpy before "
                        "running the pipeline."
                    ),
                ),
            ):
                exit_code = cluster_ani.main(
                    [
                        "--ani-matrix",
                        str(matrix),
                        "--ani-metadata",
                        str(metadata_tsv),
                        "--threshold",
                        "0.95",
                        "--outdir",
                        str(outdir),
                    ]
                )

            self.assertEqual(exit_code, 1)
            self.assertIn("SciPy is required for ANI clustering.", stderr.getvalue())

    @unittest.skipUnless(SCIPY_AVAILABLE, "SciPy is required for ANI clustering tests.")
    def test_main_writes_cluster_memberships_without_representatives(self) -> None:
        """Emit stable cluster assignments only from the ANI matrix and metadata."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            matrix = tmpdir / "fastani.matrix"
            matrix.write_text(
                "\n".join(
                    [
                        "3",
                        "fastani_inputs/ACC1.fasta",
                        "fastani_inputs/ACC2.fasta 97.5000",
                        "fastani_inputs/ACC3.fasta 80.0000 79.5000",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            metadata_tsv = tmpdir / "ani_metadata.tsv"
            metadata_tsv.write_text(
                "\n".join(
                    [
                        "accession\tmatrix_name\tpath\tassembly_level\tgcode\tcheckm2_completeness\tcheckm2_contamination\tn50\tscaffolds\tgenome_size\torganism_name\tBUSCO_bacillota_odb12",
                        "ACC2\tfastani_inputs/ACC2.fasta\tfastani_inputs/ACC2.fasta\tScaffold\t11\t95\t1.0\t50000\t4\t1000000\tTwo\tC:96.0%[S:96.0%,D:0.0%],F:2.0%,M:2.0%,n:200",
                        "ACC1\tfastani_inputs/ACC1.fasta\tfastani_inputs/ACC1.fasta\tScaffold\t11\t94\t1.5\t49000\t5\t1000000\tOne\tC:97.0%[S:97.0%,D:0.0%],F:2.0%,M:1.0%,n:200",
                        "ACC3\tfastani_inputs/ACC3.fasta\tfastani_inputs/ACC3.fasta\tScaffold\t4\t93\t2.0\t48000\t6\t1000000\tThree\tC:95.0%[S:95.0%,D:0.0%],F:3.0%,M:2.0%,n:200",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            outdir = tmpdir / "out"

            exit_code = cluster_ani.main(
                [
                    "--ani-matrix",
                    str(matrix),
                    "--ani-metadata",
                    str(metadata_tsv),
                    "--threshold",
                    "0.95",
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            rows = read_tsv(outdir / "cluster.tsv")
            self.assertEqual(
                rows,
                [
                    {
                        "Accession": "ACC1",
                        "Cluster_ID": "C000001",
                        "Matrix_Name": "fastani_inputs/ACC1.fasta",
                    },
                    {
                        "Accession": "ACC2",
                        "Cluster_ID": "C000001",
                        "Matrix_Name": "fastani_inputs/ACC2.fasta",
                    },
                    {
                        "Accession": "ACC3",
                        "Cluster_ID": "C000002",
                        "Matrix_Name": "fastani_inputs/ACC3.fasta",
                    },
                ],
            )
            self.assertFalse((outdir / "representatives.tsv").exists())


if __name__ == "__main__":
    unittest.main()
