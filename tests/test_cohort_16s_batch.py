"""Byte equivalence to the original four 16S passes, including input order."""

import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from annotation_common import write_tsv
from build_all_cohort_16s import OUTPUTS, build
from concat_best_16s import REQUIRED_INPUT_COLUMNS, Cohort16SError, build_cohort_fasta


class Cohort16SBatchTests(unittest.TestCase):
    def test_all_four_outputs_match_the_original_rules_and_order(self):
        with tempfile.TemporaryDirectory(dir="/tmp") as directory:
            root = Path(directory)
            qc_dir = root / "qc"
            qc_dir.mkdir()
            statuses = []
            inputs = []
            metadata = []
            # Deliberately non-sorted input order; include all statuses and qualities.
            for index, (status, quality, warning) in enumerate(
                (
                    ("Yes", "false", "NA"),
                    ("partial", "false", "NA"),
                    ("Yes", "true", "NA"),
                    ("partial", "true", "NA"),
                    ("No", "false", "NA"),
                    ("NA", "NA", "NA"),
                    ("Yes", "false", "contaminated"),
                )
            ):
                internal = f"S{9 - index}"
                accession = f'Accession "{index}"'
                status_file = root / (internal + "_16S_status.tsv")
                fasta = root / (internal + "_best_16S.fna")
                qc_file = qc_dir / (internal + "_checkm2_summary.tsv")
                row = {
                    "accession": accession,
                    "16S": status,
                    "best_16S_header": "native",
                    "best_16S_length": "4",
                    "warnings": "",
                    "include_in_all_best_16S": "true",
                }
                write_tsv(status_file, list(row), [row])
                write_tsv(
                    qc_file,
                    ["accession", "Low_quality"],
                    [{"accession": accession, "Low_quality": quality}],
                )
                fasta.write_text(
                    ">native\nACGT\n" if status in ("Yes", "partial") else ""
                )
                statuses.append(status_file)
                inputs.append(
                    {
                        "accession": accession,
                        "internal_id": internal,
                        "status_tsv": str(status_file),
                        "best_16s_fasta": str(fasta),
                        "qc_tsv": str(qc_file),
                    }
                )
                metadata.append({"Accession": accession, "Atypical_Warnings": warning})
            manifest = root / "inputs.tsv"
            meta = root / "metadata.tsv"
            write_tsv(manifest, REQUIRED_INPUT_COLUMNS, inputs)
            write_tsv(meta, ["Accession", "Atypical_Warnings"], metadata)
            actual = root / "actual"
            build(statuses, qc_dir, meta, actual)
            expected = root / "expected"
            expected.mkdir()
            for kind, name in OUTPUTS.items():
                fasta, table = (
                    expected / (name + ".fna"),
                    expected / (name + "_manifest.tsv"),
                )
                build_cohort_fasta(manifest, meta, fasta, table, kind)
                for path in (fasta, table):
                    self.assertEqual(
                        path.read_bytes(), (actual / path.name).read_bytes()
                    )
            with self.assertRaisesRegex(Cohort16SError, "Duplicate"):
                build(statuses * 2, qc_dir, meta, root / "duplicate")

    def test_empty_cohort_retains_empty_fastas_and_header_only_tables(self):
        with tempfile.TemporaryDirectory(dir="/tmp") as directory:
            root = Path(directory)
            build([], root, None, root)
            for name in OUTPUTS.values():
                self.assertEqual((root / (name + ".fna")).read_bytes(), b"")
                self.assertEqual(
                    len((root / (name + "_manifest.tsv")).read_text().splitlines()), 1
                )
