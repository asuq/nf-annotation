"""Tests for Barrnap-derived 16S summarisation and cohort concatenation."""

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

import concat_best_16s
import summarise_16s

CohortInputRow = tuple[str, str, Path, Path] | tuple[str, str, Path, Path, Path]


def read_status_row(path: Path) -> dict[str, str]:
    """Read a single-row status TSV."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return rows[0]


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    """Read a TSV file into row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class Summarise16STestCase(unittest.TestCase):
    """Test Barrnap-only 16S status calling and best-hit selection."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write a UTF-8 text file and return its path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def barrnap_header(
        self,
        rrna_name: str,
        contig: str,
        start0: int,
        end0: int,
        strand: str,
    ) -> str:
        """Build one Barrnap-style FASTA header."""
        return f"{rrna_name}::{contig}:{start0}-{end0}({strand})"

    def write_barrnap_fasta(
        self,
        path: Path,
        records: list[tuple[str, str]],
    ) -> Path:
        """Write a Barrnap FASTA file from header and sequence pairs."""
        lines: list[str] = []
        for header, sequence in records:
            lines.extend((f">{header}", sequence))
        return self.write_text_file(path, "\n".join(lines) + "\n")

    def test_main_prefers_lowest_score_then_longest_sequence(self) -> None:
        """Choose the best intact 16S hit with the locked tie-break order."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            best_header = self.barrnap_header("16S_rRNA", "contig1", 119, 260, "+")
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                (
                    "contig1\tbarrnap\trRNA\t1\t100\t5.0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA\n"
                    "contig1\tbarrnap\trRNA\t120\t260\t5.0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA\n"
                    "contig1\tbarrnap\trRNA\t300\t360\t7.0\t+\t.\tName=23S_rRNA;product=23S ribosomal RNA"
                )
                + "\n",
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [
                    (
                        self.barrnap_header("16S_rRNA", "contig1", 0, 100, "+"),
                        "A" * 100,
                    ),
                    (best_header, "C" * 141),
                    (
                        self.barrnap_header("23S_rRNA", "contig1", 299, 360, "+"),
                        "G" * 61,
                    ),
                ],
            )
            outdir = tmpdir / "out"

            exit_code = summarise_16s.main(
                [
                    "--accession",
                    "ACC1",
                    "--rrna-gff",
                    str(gff),
                    "--rrna-fasta",
                    str(fasta),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            status_row = read_status_row(outdir / "16S_status.tsv")
            self.assertEqual(
                tuple(status_row),
                ("accession", "16S", "best_16S_header", "best_16S_length", "warnings"),
            )
            self.assertEqual(status_row["16S"], "Yes")
            self.assertEqual(status_row["best_16S_header"], best_header)
            self.assertEqual(status_row["best_16S_length"], "141")
            self.assertEqual(status_row["warnings"], "")
            best_fasta = (outdir / "best_16S.fna").read_text(encoding="utf-8")
            self.assertTrue(best_fasta.startswith(f">{best_header}\n"))

    def test_main_reports_partial_when_only_partial_hits_exist(self) -> None:
        """Emit `partial` when all 16S hits are partial."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            partial_header = self.barrnap_header(
                "16S_rRNA", "KB890753.1", 102, 1222, "+"
            )
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                (
                    "KB890753.1\tbarrnap\trRNA\t103\t1222\t5.7e-275\t+\t.\t"
                    "Name=16S_rRNA;product=16S ribosomal RNA (partial);"
                    "note=aligned only 70 percent of the 16S ribosomal RNA\n"
                ),
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [(partial_header, "A" * 1120)],
            )
            outdir = tmpdir / "out"

            exit_code = summarise_16s.main(
                [
                    "--accession",
                    "ACC2",
                    "--rrna-gff",
                    str(gff),
                    "--rrna-fasta",
                    str(fasta),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            status_row = read_status_row(outdir / "16S_status.tsv")
            self.assertEqual(status_row["16S"], "partial")
            self.assertEqual(status_row["best_16S_header"], partial_header)

    def test_main_reports_no_when_no_16s_hits_exist(self) -> None:
        """Emit `No` when Barrnap found only non-16S rRNA."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "contig1\tbarrnap\trRNA\t1\t100\t5.0\t+\t.\tName=23S_rRNA;product=23S ribosomal RNA\n",
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [(self.barrnap_header("23S_rRNA", "contig1", 0, 100, "+"), "A" * 100)],
            )
            outdir = tmpdir / "out"

            summarise_16s.main(
                [
                    "--accession",
                    "ACC5",
                    "--rrna-gff",
                    str(gff),
                    "--rrna-fasta",
                    str(fasta),
                    "--outdir",
                    str(outdir),
                ]
            )

            status_row = read_status_row(outdir / "16S_status.tsv")
            self.assertEqual(status_row["16S"], "No")
            self.assertEqual(status_row["best_16S_header"], "NA")

    def test_main_returns_na_when_gff_and_fasta_do_not_match(self) -> None:
        """Gracefully degrade to NA when Barrnap 16S GFF and FASTA differ."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                "contig1\tbarrnap\trRNA\t1\t100\t5.0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA\n",
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [(self.barrnap_header("23S_rRNA", "contig1", 0, 100, "+"), "A" * 100)],
            )
            outdir = tmpdir / "out"

            exit_code = summarise_16s.main(
                [
                    "--accession",
                    "ACC6",
                    "--rrna-gff",
                    str(gff),
                    "--rrna-fasta",
                    str(fasta),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            status_row = read_status_row(outdir / "16S_status.tsv")
            self.assertEqual(status_row["16S"], "NA")
            self.assertEqual(status_row["warnings"], "invalid_barrnap_output")
            self.assertEqual((outdir / "best_16S.fna").read_text(encoding="utf-8"), "")

    def test_main_matches_mixed_rrna_order_by_key(self) -> None:
        """Match 16S records by Barrnap key even when mixed rRNA order differs."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            best_header = self.barrnap_header("16S_rRNA", "contig1", 319, 470, "+")
            gff = self.write_text_file(
                tmpdir / "rrna.gff",
                (
                    "contig1\tbarrnap\trRNA\t1\t100\t8.0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA\n"
                    "contig1\tbarrnap\trRNA\t120\t220\t1.0\t+\t.\tName=23S_rRNA;product=23S ribosomal RNA\n"
                    "contig1\tbarrnap\trRNA\t240\t300\t1.0\t+\t.\tName=5S_rRNA;product=5S ribosomal RNA\n"
                    "contig1\tbarrnap\trRNA\t320\t470\t4.0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA"
                )
                + "\n",
            )
            fasta = self.write_barrnap_fasta(
                tmpdir / "rrna.fa",
                [
                    (
                        self.barrnap_header("23S_rRNA", "contig1", 119, 220, "+"),
                        "G" * 101,
                    ),
                    (best_header, "C" * 151),
                    (
                        self.barrnap_header("5S_rRNA", "contig1", 239, 300, "+"),
                        "T" * 61,
                    ),
                    (
                        self.barrnap_header("16S_rRNA", "contig1", 0, 100, "+"),
                        "A" * 100,
                    ),
                ],
            )
            outdir = tmpdir / "out"

            exit_code = summarise_16s.main(
                [
                    "--accession",
                    "ACC7",
                    "--rrna-gff",
                    str(gff),
                    "--rrna-fasta",
                    str(fasta),
                    "--outdir",
                    str(outdir),
                ]
            )

            self.assertEqual(exit_code, 0)
            status_row = read_status_row(outdir / "16S_status.tsv")
            self.assertEqual(status_row["16S"], "Yes")
            self.assertEqual(status_row["best_16S_header"], best_header)


class ConcatBest16STestCase(unittest.TestCase):
    """Test cohort-level 16S concatenation and atypical inclusion rules."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write a UTF-8 text file and return the path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def write_status_and_fasta(
        self,
        base_dir: Path,
        *,
        accession: str,
        status_value: str,
        header: str,
        length: int,
        warnings: str = "",
    ) -> tuple[Path, Path]:
        """Write one per-sample 16S status file and best FASTA."""
        status_path = self.write_text_file(
            base_dir / "16S_status.tsv",
            "\n".join(
                [
                    "accession\t16S\tbest_16S_header\tbest_16S_length\twarnings",
                    f"{accession}\t{status_value}\t{header}\t{length}\t{warnings}",
                ]
            )
            + "\n",
        )
        best_fasta = self.write_text_file(
            base_dir / "best_16S.fna",
            f">{header}\n" + "A" * length + "\n",
        )
        return status_path, best_fasta

    def write_qc_summary(
        self,
        path: Path,
        *,
        accession: str,
        low_quality: str = "false",
    ) -> Path:
        """Write one per-sample CheckM2 QC summary."""
        return self.write_text_file(
            path,
            "\n".join(
                [
                    "accession\tLow_quality",
                    f"{accession}\t{low_quality}",
                ]
            )
            + "\n",
        )

    def write_cohort_inputs(
        self,
        path: Path,
        rows: list[CohortInputRow],
    ) -> Path:
        """Write one cohort input manifest."""
        lines = ["accession\tinternal_id\tstatus_tsv\tbest_16s_fasta\tqc_tsv"]
        for row in rows:
            accession, internal_id, status_path, best_fasta = row[:4]
            if len(row) == 5:
                qc_path = row[4]
            else:
                qc_path = self.write_qc_summary(
                    path.parent / f"{internal_id}_checkm2_summary.tsv",
                    accession=accession,
                )
            lines.append(
                f"{accession}\t{internal_id}\t{status_path}\t{best_fasta}\t{qc_path}"
            )
        return self.write_text_file(path, "\n".join(lines) + "\n")

    def test_concat_best_16s_includes_non_atypical_intact_rows(self) -> None:
        """Build the intact cohort from ordinary non-atypical samples."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            include_dir = tmpdir / "include"
            include_dir.mkdir()
            status_path, best_fasta = self.write_status_and_fasta(
                include_dir,
                accession="ACC_OK",
                status_value="Yes",
                header="hit_ok",
                length=100,
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tAtypical_Warnings\nACC_OK\tNA\n",
            )
            cohort_inputs = self.write_cohort_inputs(
                tmpdir / "cohort_inputs.tsv",
                [("ACC_OK", "acc_ok_id", status_path, best_fasta)],
            )
            output_fasta = tmpdir / "all_best_16S.fna"
            output_manifest = tmpdir / "all_best_16S_manifest.tsv"

            exit_code = concat_best_16s.main(
                [
                    "--inputs",
                    str(cohort_inputs),
                    "--metadata",
                    str(metadata),
                    "--output-fasta",
                    str(output_fasta),
                    "--output-manifest",
                    str(output_manifest),
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertEqual(
                output_fasta.read_text(encoding="utf-8"),
                ">acc_ok_id\n" + "A" * 100 + "\n",
            )
            manifest_rows = read_tsv_rows(output_manifest)
            self.assertEqual(manifest_rows[0]["accession"], "ACC_OK")
            self.assertEqual(manifest_rows[0]["best_16S_header"], "hit_ok")
            self.assertEqual(manifest_rows[0]["include_in_all_best_16S"], "true")

    def test_concat_best_16s_allows_unverified_source_exception(self) -> None:
        """Keep the locked unverified-source exception in the intact cohort."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            sample_dir = tmpdir / "sample"
            sample_dir.mkdir()
            status_path, best_fasta = self.write_status_and_fasta(
                sample_dir,
                accession="ACC_EXCEPTION",
                status_value="Yes",
                header="hit_exception",
                length=90,
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tAtypical_Warnings\nACC_EXCEPTION\tUnverified source organism\n",
            )
            cohort_inputs = self.write_cohort_inputs(
                tmpdir / "cohort_inputs.tsv",
                [("ACC_EXCEPTION", "acc_exception_id", status_path, best_fasta)],
            )

            exit_code = concat_best_16s.main(
                [
                    "--inputs",
                    str(cohort_inputs),
                    "--metadata",
                    str(metadata),
                    "--output-fasta",
                    str(tmpdir / "all_best_16S.fna"),
                    "--output-manifest",
                    str(tmpdir / "all_best_16S_manifest.tsv"),
                ]
            )

            self.assertEqual(exit_code, 0)
            manifest_rows = read_tsv_rows(tmpdir / "all_best_16S_manifest.tsv")
            self.assertEqual(
                [row["accession"] for row in manifest_rows], ["ACC_EXCEPTION"]
            )

    def test_concat_best_16s_excludes_other_atypical_intact_rows(self) -> None:
        """Exclude intact 16S rows when the atypical reason is not the exception."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            sample_dir = tmpdir / "sample"
            sample_dir.mkdir()
            status_path, best_fasta = self.write_status_and_fasta(
                sample_dir,
                accession="ACC_ATYP",
                status_value="Yes",
                header="hit_atyp",
                length=80,
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tAtypical_Warnings\nACC_ATYP\tcell culture adapted\n",
            )
            cohort_inputs = self.write_cohort_inputs(
                tmpdir / "cohort_inputs.tsv",
                [("ACC_ATYP", "acc_atyp_id", status_path, best_fasta)],
            )

            exit_code = concat_best_16s.main(
                [
                    "--inputs",
                    str(cohort_inputs),
                    "--metadata",
                    str(metadata),
                    "--output-fasta",
                    str(tmpdir / "all_best_16S.fna"),
                    "--output-manifest",
                    str(tmpdir / "all_best_16S_manifest.tsv"),
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertEqual(
                (tmpdir / "all_best_16S.fna").read_text(encoding="utf-8"), ""
            )
            self.assertEqual(read_tsv_rows(tmpdir / "all_best_16S_manifest.tsv"), [])

    def test_concat_best_16s_excludes_mixed_atypical_reasons(self) -> None:
        """Exclude intact 16S rows when unverified source organism is not the sole reason."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            sample_dir = tmpdir / "sample"
            sample_dir.mkdir()
            status_path, best_fasta = self.write_status_and_fasta(
                sample_dir,
                accession="ACC_MIXED",
                status_value="Yes",
                header="hit_mixed",
                length=85,
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tAtypical_Warnings\nACC_MIXED\tgenome length too small, unverified source organism\n",
            )
            cohort_inputs = self.write_cohort_inputs(
                tmpdir / "cohort_inputs.tsv",
                [("ACC_MIXED", "acc_mixed_id", status_path, best_fasta)],
            )

            exit_code = concat_best_16s.main(
                [
                    "--inputs",
                    str(cohort_inputs),
                    "--metadata",
                    str(metadata),
                    "--output-fasta",
                    str(tmpdir / "all_best_16S.fna"),
                    "--output-manifest",
                    str(tmpdir / "all_best_16S_manifest.tsv"),
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertEqual(
                (tmpdir / "all_best_16S.fna").read_text(encoding="utf-8"), ""
            )
            self.assertEqual(read_tsv_rows(tmpdir / "all_best_16S_manifest.tsv"), [])

    def test_concat_best_16s_builds_partial_cohort_including_atypical_samples(
        self,
    ) -> None:
        """Build the partial cohort from all partial rows, including atypical ones."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            partial_dir = tmpdir / "partial"
            atypical_partial_dir = tmpdir / "atypical_partial"
            partial_dir.mkdir()
            atypical_partial_dir.mkdir()

            partial_status, partial_best = self.write_status_and_fasta(
                partial_dir,
                accession="ACC_PART",
                status_value="partial",
                header="hit_partial",
                length=90,
            )
            atypical_partial_status, atypical_partial_best = (
                self.write_status_and_fasta(
                    atypical_partial_dir,
                    accession="ACC_ATYP_PART",
                    status_value="partial",
                    header="hit_atyp_partial",
                    length=80,
                    warnings="atypical",
                )
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tAtypical_Warnings\nACC_PART\tNA\nACC_ATYP_PART\tcell culture adapted"
                + "\n",
            )
            cohort_inputs = self.write_cohort_inputs(
                tmpdir / "cohort_inputs.tsv",
                [
                    ("ACC_PART", "acc_part_id", partial_status, partial_best),
                    (
                        "ACC_ATYP_PART",
                        "acc_atyp_part_id",
                        atypical_partial_status,
                        atypical_partial_best,
                    ),
                ],
            )

            exit_code = concat_best_16s.main(
                [
                    "--inputs",
                    str(cohort_inputs),
                    "--metadata",
                    str(metadata),
                    "--cohort-kind",
                    "partial",
                    "--output-fasta",
                    str(tmpdir / "all_partial_16S.fna"),
                    "--output-manifest",
                    str(tmpdir / "all_partial_16S_manifest.tsv"),
                ]
            )

            self.assertEqual(exit_code, 0)
            manifest_rows = read_tsv_rows(tmpdir / "all_partial_16S_manifest.tsv")
            self.assertEqual(
                [row["accession"] for row in manifest_rows],
                ["ACC_PART", "ACC_ATYP_PART"],
            )
            fasta_text = (tmpdir / "all_partial_16S.fna").read_text(encoding="utf-8")
            self.assertTrue(fasta_text.startswith(">acc_part_id\n"))
            self.assertIn(">acc_atyp_part_id\n", fasta_text)
            self.assertEqual(
                [row["best_16S_header"] for row in manifest_rows],
                ["hit_partial", "hit_atyp_partial"],
            )

    def test_concat_best_16s_moves_low_quality_intact_rows_out(self) -> None:
        """Exclude low-quality intact rows from the standard intact cohort."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            sample_dir = tmpdir / "sample"
            sample_dir.mkdir()
            status_path, best_fasta = self.write_status_and_fasta(
                sample_dir,
                accession="ACC_LOW",
                status_value="Yes",
                header="hit_low",
                length=70,
            )
            qc_path = self.write_qc_summary(
                sample_dir / "ACC_LOW_checkm2_summary.tsv",
                accession="ACC_LOW",
                low_quality="true",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tAtypical_Warnings\nACC_LOW\tNA\n",
            )
            cohort_inputs = self.write_cohort_inputs(
                tmpdir / "cohort_inputs.tsv",
                [("ACC_LOW", "acc_low_id", status_path, best_fasta, qc_path)],
            )

            intact_exit_code = concat_best_16s.main(
                [
                    "--inputs",
                    str(cohort_inputs),
                    "--metadata",
                    str(metadata),
                    "--output-fasta",
                    str(tmpdir / "all_best_16S.fna"),
                    "--output-manifest",
                    str(tmpdir / "all_best_16S_manifest.tsv"),
                ]
            )
            low_quality_exit_code = concat_best_16s.main(
                [
                    "--inputs",
                    str(cohort_inputs),
                    "--metadata",
                    str(metadata),
                    "--cohort-kind",
                    "low_quality_intact",
                    "--output-fasta",
                    str(tmpdir / "low_quality_best_16S.fna"),
                    "--output-manifest",
                    str(tmpdir / "low_quality_best_16S_manifest.tsv"),
                ]
            )

            self.assertEqual(intact_exit_code, 0)
            self.assertEqual(low_quality_exit_code, 0)
            self.assertEqual(read_tsv_rows(tmpdir / "all_best_16S_manifest.tsv"), [])
            manifest_rows = read_tsv_rows(tmpdir / "low_quality_best_16S_manifest.tsv")
            self.assertEqual([row["accession"] for row in manifest_rows], ["ACC_LOW"])
            self.assertEqual(manifest_rows[0]["low_quality"], "true")
            self.assertEqual(manifest_rows[0]["include_in_all_best_16S"], "false")

    def test_concat_best_16s_moves_low_quality_partial_rows_out(self) -> None:
        """Exclude low-quality partial rows from the standard partial cohort."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            sample_dir = tmpdir / "sample"
            sample_dir.mkdir()
            status_path, best_fasta = self.write_status_and_fasta(
                sample_dir,
                accession="ACC_LOW_PART",
                status_value="partial",
                header="hit_low_part",
                length=60,
            )
            qc_path = self.write_qc_summary(
                sample_dir / "ACC_LOW_PART_checkm2_summary.tsv",
                accession="ACC_LOW_PART",
                low_quality="true",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tAtypical_Warnings\nACC_LOW_PART\tNA\n",
            )
            cohort_inputs = self.write_cohort_inputs(
                tmpdir / "cohort_inputs.tsv",
                [("ACC_LOW_PART", "acc_low_part_id", status_path, best_fasta, qc_path)],
            )

            partial_exit_code = concat_best_16s.main(
                [
                    "--inputs",
                    str(cohort_inputs),
                    "--metadata",
                    str(metadata),
                    "--cohort-kind",
                    "partial",
                    "--output-fasta",
                    str(tmpdir / "all_partial_16S.fna"),
                    "--output-manifest",
                    str(tmpdir / "all_partial_16S_manifest.tsv"),
                ]
            )
            low_quality_exit_code = concat_best_16s.main(
                [
                    "--inputs",
                    str(cohort_inputs),
                    "--metadata",
                    str(metadata),
                    "--cohort-kind",
                    "low_quality_partial",
                    "--output-fasta",
                    str(tmpdir / "low_quality_partial_16S.fna"),
                    "--output-manifest",
                    str(tmpdir / "low_quality_partial_16S_manifest.tsv"),
                ]
            )

            self.assertEqual(partial_exit_code, 0)
            self.assertEqual(low_quality_exit_code, 0)
            self.assertEqual(read_tsv_rows(tmpdir / "all_partial_16S_manifest.tsv"), [])
            manifest_rows = read_tsv_rows(
                tmpdir / "low_quality_partial_16S_manifest.tsv"
            )
            self.assertEqual(
                [row["accession"] for row in manifest_rows], ["ACC_LOW_PART"]
            )
            self.assertEqual(manifest_rows[0]["low_quality"], "true")

    def test_concat_best_16s_includes_low_quality_atypical_rows(self) -> None:
        """Keep low-quality audit outputs independent of atypical filtering."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            intact_dir = tmpdir / "intact"
            partial_dir = tmpdir / "partial"
            intact_dir.mkdir()
            partial_dir.mkdir()
            intact_status, intact_fasta = self.write_status_and_fasta(
                intact_dir,
                accession="ACC_ATYP_LOW",
                status_value="Yes",
                header="hit_atyp_low",
                length=75,
            )
            partial_status, partial_fasta = self.write_status_and_fasta(
                partial_dir,
                accession="ACC_ATYP_LOW_PART",
                status_value="partial",
                header="hit_atyp_low_part",
                length=65,
            )
            intact_qc = self.write_qc_summary(
                intact_dir / "ACC_ATYP_LOW_checkm2_summary.tsv",
                accession="ACC_ATYP_LOW",
                low_quality="true",
            )
            partial_qc = self.write_qc_summary(
                partial_dir / "ACC_ATYP_LOW_PART_checkm2_summary.tsv",
                accession="ACC_ATYP_LOW_PART",
                low_quality="true",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tAtypical_Warnings\nACC_ATYP_LOW\tcell culture adapted\nACC_ATYP_LOW_PART\tcell culture adapted"
                + "\n",
            )
            cohort_inputs = self.write_cohort_inputs(
                tmpdir / "cohort_inputs.tsv",
                [
                    (
                        "ACC_ATYP_LOW",
                        "acc_atyp_low_id",
                        intact_status,
                        intact_fasta,
                        intact_qc,
                    ),
                    (
                        "ACC_ATYP_LOW_PART",
                        "acc_atyp_low_part_id",
                        partial_status,
                        partial_fasta,
                        partial_qc,
                    ),
                ],
            )

            intact_exit_code = concat_best_16s.main(
                [
                    "--inputs",
                    str(cohort_inputs),
                    "--metadata",
                    str(metadata),
                    "--cohort-kind",
                    "low_quality_intact",
                    "--output-fasta",
                    str(tmpdir / "low_quality_best_16S.fna"),
                    "--output-manifest",
                    str(tmpdir / "low_quality_best_16S_manifest.tsv"),
                ]
            )
            partial_exit_code = concat_best_16s.main(
                [
                    "--inputs",
                    str(cohort_inputs),
                    "--metadata",
                    str(metadata),
                    "--cohort-kind",
                    "low_quality_partial",
                    "--output-fasta",
                    str(tmpdir / "low_quality_partial_16S.fna"),
                    "--output-manifest",
                    str(tmpdir / "low_quality_partial_16S_manifest.tsv"),
                ]
            )

            self.assertEqual(intact_exit_code, 0)
            self.assertEqual(partial_exit_code, 0)
            low_quality_best_rows = read_tsv_rows(
                tmpdir / "low_quality_best_16S_manifest.tsv"
            )
            low_quality_partial_rows = read_tsv_rows(
                tmpdir / "low_quality_partial_16S_manifest.tsv"
            )
            self.assertEqual(
                [row["accession"] for row in low_quality_best_rows],
                ["ACC_ATYP_LOW"],
            )
            self.assertEqual(
                [row["accession"] for row in low_quality_partial_rows],
                ["ACC_ATYP_LOW_PART"],
            )

    def test_concat_best_16s_keeps_na_low_quality_in_standard_outputs(self) -> None:
        """Treat Low_quality=NA as a standard eligible QC value."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            sample_dir = tmpdir / "sample"
            sample_dir.mkdir()
            status_path, best_fasta = self.write_status_and_fasta(
                sample_dir,
                accession="ACC_NA_QC",
                status_value="Yes",
                header="hit_na_qc",
                length=70,
            )
            qc_path = self.write_qc_summary(
                sample_dir / "ACC_NA_QC_checkm2_summary.tsv",
                accession="ACC_NA_QC",
                low_quality="NA",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tAtypical_Warnings\nACC_NA_QC\tNA\n",
            )
            cohort_inputs = self.write_cohort_inputs(
                tmpdir / "cohort_inputs.tsv",
                [("ACC_NA_QC", "acc_na_qc_id", status_path, best_fasta, qc_path)],
            )

            exit_code = concat_best_16s.main(
                [
                    "--inputs",
                    str(cohort_inputs),
                    "--metadata",
                    str(metadata),
                    "--output-fasta",
                    str(tmpdir / "all_best_16S.fna"),
                    "--output-manifest",
                    str(tmpdir / "all_best_16S_manifest.tsv"),
                ]
            )

            self.assertEqual(exit_code, 0)
            manifest_rows = read_tsv_rows(tmpdir / "all_best_16S_manifest.tsv")
            self.assertEqual([row["accession"] for row in manifest_rows], ["ACC_NA_QC"])
            self.assertEqual(manifest_rows[0]["low_quality"], "NA")

    def test_concat_best_16s_fails_when_qc_accession_mismatches(self) -> None:
        """Fail when the cohort input and QC summary accessions disagree."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            status_path, best_fasta = self.write_status_and_fasta(
                tmpdir,
                accession="ACC_BAD_QC",
                status_value="Yes",
                header="hit_bad_qc",
                length=70,
            )
            qc_path = self.write_qc_summary(
                tmpdir / "ACC_OTHER_checkm2_summary.tsv",
                accession="ACC_OTHER",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tAtypical_Warnings\nACC_BAD_QC\tNA\n",
            )
            cohort_inputs = self.write_cohort_inputs(
                tmpdir / "cohort_inputs.tsv",
                [("ACC_BAD_QC", "acc_bad_qc_id", status_path, best_fasta, qc_path)],
            )

            exit_code = concat_best_16s.main(
                [
                    "--inputs",
                    str(cohort_inputs),
                    "--metadata",
                    str(metadata),
                    "--output-fasta",
                    str(tmpdir / "all_best_16S.fna"),
                    "--output-manifest",
                    str(tmpdir / "all_best_16S_manifest.tsv"),
                ]
            )

            self.assertEqual(exit_code, 1)

    def test_concat_best_16s_fails_when_qc_file_is_missing(self) -> None:
        """Fail when a cohort input row points to a missing QC summary."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            status_path, best_fasta = self.write_status_and_fasta(
                tmpdir,
                accession="ACC_MISSING_QC",
                status_value="Yes",
                header="hit_missing_qc",
                length=70,
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tAtypical_Warnings\nACC_MISSING_QC\tNA\n",
            )
            cohort_inputs = self.write_cohort_inputs(
                tmpdir / "cohort_inputs.tsv",
                [
                    (
                        "ACC_MISSING_QC",
                        "acc_missing_qc_id",
                        status_path,
                        best_fasta,
                        tmpdir / "missing_qc.tsv",
                    )
                ],
            )

            exit_code = concat_best_16s.main(
                [
                    "--inputs",
                    str(cohort_inputs),
                    "--metadata",
                    str(metadata),
                    "--output-fasta",
                    str(tmpdir / "all_best_16S.fna"),
                    "--output-manifest",
                    str(tmpdir / "all_best_16S_manifest.tsv"),
                ]
            )

            self.assertEqual(exit_code, 1)

    def test_concat_best_16s_fails_when_included_fasta_is_empty(self) -> None:
        """Fail when a cohort-included sample has no FASTA content to append."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            status_path = self.write_text_file(
                tmpdir / "16S_status.tsv",
                "accession\t16S\tbest_16S_header\tbest_16S_length\twarnings\nACC_BAD\tYes\thit_bad\t100\t"
                + "\n",
            )
            best_fasta = self.write_text_file(tmpdir / "best_16S.fna", "")
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tAtypical_Warnings\nACC_BAD\tNA\n",
            )
            cohort_inputs = self.write_cohort_inputs(
                tmpdir / "cohort_inputs.tsv",
                [("ACC_BAD", "acc_bad_id", status_path, best_fasta)],
            )

            exit_code = concat_best_16s.main(
                [
                    "--inputs",
                    str(cohort_inputs),
                    "--metadata",
                    str(metadata),
                    "--output-fasta",
                    str(tmpdir / "all_best_16S.fna"),
                    "--output-manifest",
                    str(tmpdir / "all_best_16S_manifest.tsv"),
                ]
            )

            self.assertEqual(exit_code, 1)

    def test_concat_best_16s_fails_when_internal_id_column_is_missing(self) -> None:
        """Fail when the cohort input manifest lacks the internal_id column."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            status_path, best_fasta = self.write_status_and_fasta(
                tmpdir,
                accession="ACC_BAD",
                status_value="Yes",
                header="hit_bad",
                length=100,
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tAtypical_Warnings\nACC_BAD\tNA\n",
            )
            cohort_inputs = self.write_text_file(
                tmpdir / "cohort_inputs.tsv",
                "\n".join(
                    [
                        "accession\tstatus_tsv\tbest_16s_fasta",
                        f"ACC_BAD\t{status_path}\t{best_fasta}",
                    ]
                )
                + "\n",
            )

            exit_code = concat_best_16s.main(
                [
                    "--inputs",
                    str(cohort_inputs),
                    "--metadata",
                    str(metadata),
                    "--output-fasta",
                    str(tmpdir / "all_best_16S.fna"),
                    "--output-manifest",
                    str(tmpdir / "all_best_16S_manifest.tsv"),
                ]
            )

            self.assertEqual(exit_code, 1)


if __name__ == "__main__":
    unittest.main()
