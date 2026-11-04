"""Removed work stays removed while genome reports and integrity stay exact."""

import sys
import threading
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import aggregate_annotations
import annotation_common
import report_genomes
from annotation_common import AnnotationError, read_json, write_json
from fingerprinting.source_fingerprint import FingerprintScanner, checked_digest
from reuse_source_validation import validate, validate_or_reuse
import test_annotation_aggregation_scaling as aggregation_fixture
import test_parallel_fingerprinting as fingerprint_fixture
import test_annotation_source as source_fixture


class StreamliningTests(unittest.TestCase):
    def test_ani_report_preserves_genome_tables_without_reading_protein_payloads(self):
        case = aggregation_fixture.AnnotationAggregationScalingTests()
        case.setUp()
        self.addCleanup(case.doCleanups)
        args = case.cohort(["B", "A"], tools=())
        full, ani = case.root / "full", case.root / "ani"
        aggregate_annotations.aggregate(*args, full)
        with patch.object(
            annotation_common,
            "bundle_proteins",
            side_effect=AssertionError("protein read"),
        ):
            manifest = report_genomes.report(
                *args[:3], [p / "bundle.json" for p in args[3]], ani
            )
        for name in (
            "master_table.tsv",
            "sample_status.tsv",
            "annotation_status.tsv",
            "annotation_provenance.tsv",
        ):
            self.assertEqual(
                (full / "tables" / name).read_bytes(),
                (ani / "tables" / name).read_bytes(),
                name,
            )
        self.assertEqual(len(manifest["tables"]), 4)
        self.assertFalse((ani / "tables/protein_manifest.tsv").exists())
        self.assertEqual(
            read_json(ani / "annotation_results.json")["report_scope"], "genome"
        )

    def test_scan_digests_supply_semantic_validation_without_second_payload_reads(self):
        case = source_fixture.AnnotationSourceTests()
        case.setUp()
        self.addCleanup(case.doCleanups)
        with FingerprintScanner(2) as scanner:
            record = scanner.fingerprint(case.source)
            with patch(
                "annotation_common.hashlib.file_digest",
                side_effect=AssertionError("second payload read"),
            ):
                validate([case.source], [record], checksum=scanner.checksum)
            target = case.source / "tables/master_table.tsv"
            target.write_bytes(target.read_bytes() + b"\n")
            with self.assertRaisesRegex(
                AnnotationError, "changed after fingerprinting"
            ):
                scanner.checksum(target)

    def test_semantic_receipt_reuse_requires_fresh_content_and_code(self):
        case = source_fixture.AnnotationSourceTests()
        case.setUp()
        self.addCleanup(case.doCleanups)
        previous = case.root / "receipt.json"
        with FingerprintScanner(2) as scanner:
            records = [scanner.fingerprint(case.source)]
            receipts = validate([case.source], records, checksum=scanner.checksum)
        write_json(previous, {"receipts": receipts})
        with FingerprintScanner(2) as scanner:
            records = [scanner.fingerprint(case.source)]
            with patch("reuse_source_validation.validate", side_effect=AssertionError("revalidated")):
                _, reused = validate_or_reuse([case.source], records, previous, checksum=scanner.checksum)
                self.assertTrue(reused)
        target = case.source / "tables/master_table.tsv"
        target.write_bytes(target.read_bytes() + b"\n")
        with FingerprintScanner(2) as scanner:
            records = [scanner.fingerprint(case.source)]
            with self.assertRaises(AnnotationError):
                validate_or_reuse([case.source], records, previous, checksum=scanner.checksum)

    def test_slow_first_file_does_not_stop_reader_replenishment(self):
        case = fingerprint_fixture.ParallelFingerprintTests()
        case.setUp()
        self.addCleanup(case.doCleanups)
        paths = sorted(case.sources[0].glob("samples/A/*"))
        completed = threading.Event()
        entered = []
        lock = threading.Lock()

        def reader(path, before):
            if path == paths[0]:
                if not completed.wait(5):
                    raise AssertionError("Readers stopped behind first file")
            else:
                with lock:
                    entered.append(path)
                    if len(entered) >= 12:
                        completed.set()
            return checked_digest(path, before)

        with patch(
            "fingerprinting.source_fingerprint.checked_digest", side_effect=reader
        ):
            with FingerprintScanner(2) as scanner:
                scanner.metrics.update(current_files_hashed=0, current_bytes_hashed=0)
                observed = list(scanner.hashes([(p, p.stat()) for p in paths]))
                self.assertEqual([p for p, _, _ in observed], paths)
                self.assertLessEqual(scanner.metrics["peak_pending_reads"], 4)


if __name__ == "__main__":
    unittest.main()
