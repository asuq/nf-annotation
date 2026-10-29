"""Parallel reads retain legacy identities, bounds and rejection of changing inputs."""

import hashlib
import json
import os
import sys
import tempfile
import threading
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from annotation_common import AnnotationError, digest, identity
from fingerprinting.source_fingerprint import (
    checked_digest,
    fingerprint,
    fingerprint_sources,
)


class ParallelFingerprintTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(dir="/tmp")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        self.sources = []
        for index in range(3):
            root = self.root / f"source{index}"
            (root / "tables").mkdir(parents=True)
            (root / "samples/A").mkdir(parents=True)
            (root / "annotation_results.json").write_text("{}\n")
            for number in range(24):
                (root / f"samples/A/{number:02}.fasta").write_bytes(
                    bytes([65 + index]) * (number + 1) * 1024
                )
            self.sources.append(root)

    def legacy(self, root):
        # v0.4.5 byte contract: sorted relative paths, sizes and complete SHA-256s.
        paths = sorted([root / "annotation_results.json", *root.glob("samples/A/*")])
        tree = hashlib.sha256()
        for path in paths:
            tree.update(
                identity(
                    [
                        path.relative_to(root).as_posix(),
                        path.stat().st_size,
                        digest(path),
                    ]
                ).encode()
            )
        return {
            "root": str(root),
            "sha256": tree.hexdigest(),
            "files": len(paths),
            "bytes": sum(path.stat().st_size for path in paths),
        }

    def test_workers_preserve_legacy_records_bytes_and_source_order(self):
        order = [self.sources[2], self.sources[0], self.sources[1]]
        expected = [self.legacy(root) for root in order]
        for workers in (1, 2, 4, 8):
            with self.subTest(workers=workers):
                records, metrics = fingerprint_sources(order, workers=workers)
                self.assertEqual(records, expected)
                self.assertEqual(json.dumps(records), json.dumps(expected))
                self.assertEqual(metrics["workers"], workers)
                self.assertEqual(metrics["files_hashed"], 75)
                self.assertLessEqual(metrics["peak_pending_reads"], 2 * workers)

    def test_readers_are_parallel_and_bounded_without_nested_pools(self):
        lock = threading.Lock()
        active = peak = 0
        barrier = threading.Barrier(4)
        entered = 0

        def reader(path, before):
            nonlocal active, peak, entered
            with lock:
                active += 1
                entered += 1
                initial = entered <= 4
                peak = max(peak, active)
            try:
                if initial:
                    barrier.wait(timeout=5)
                return checked_digest(path, before)
            finally:
                with lock:
                    active -= 1

        with patch(
            "fingerprinting.source_fingerprint.checked_digest", side_effect=reader
        ):
            records, metrics = fingerprint_sources(self.sources, workers=4)
        self.assertEqual(peak, 4)
        self.assertEqual(len(records), 3)
        self.assertEqual(metrics["peak_pending_reads"], 8)

    def test_same_size_restored_mtime_edit_during_read_fails(self):
        victim = self.sources[0] / "samples/A/00.fasta"
        original = hashlib.file_digest

        def mutate(handle, algorithm):
            result = original(handle, algorithm)
            info = os.fstat(handle.fileno())
            if info.st_ino == victim.stat().st_ino:
                victim.write_bytes(b"Z" * info.st_size)
                os.utime(victim, ns=(info.st_atime_ns, info.st_mtime_ns))
            return result

        with (
            patch(
                "fingerprinting.source_fingerprint.hashlib.file_digest",
                side_effect=mutate,
            ),
            self.assertRaisesRegex(AnnotationError, "changed during fingerprinting"),
        ):
            fingerprint(self.sources[0], workers=4)

    def test_changed_directory_membership_and_leaf_replacement_fail(self):
        root = self.sources[0]

        def add_file(path, before):
            result = checked_digest(path, before)
            (root / "samples/A/new").write_text("new")
            return result

        with (
            patch(
                "fingerprinting.source_fingerprint.checked_digest", side_effect=add_file
            ),
            self.assertRaisesRegex(AnnotationError, "changed during fingerprinting"),
        ):
            fingerprint(root, workers=4)
        victim = root / "samples/A/00.fasta"
        before = victim.stat()
        victim.unlink()
        victim.symlink_to(root / "samples/A/01.fasta")
        with self.assertRaises(OSError):
            checked_digest(victim, before)

    def test_progress_reports_while_a_read_is_blocked_and_finishes_atomically(self):
        progress = self.root / "progress.json"
        entered = threading.Event()
        release = threading.Event()
        original = checked_digest

        def blocked(path, before):
            entered.set()
            release.wait(timeout=5)
            return original(path, before)

        results = []
        with patch(
            "fingerprinting.source_fingerprint.checked_digest", side_effect=blocked
        ):
            worker = threading.Thread(
                target=lambda: results.append(
                    fingerprint_sources(self.sources, workers=2, progress=progress)
                )
            )
            worker.start()
            self.assertTrue(entered.wait(timeout=5))
            report = json.loads(progress.read_text())
            self.assertEqual(report["phase"], "hashing")
            self.assertEqual(report["workers"], 2)
            self.assertGreater(report["current_bytes_total"], 0)
            release.set()
            worker.join(timeout=10)
        self.assertFalse(worker.is_alive())
        report = json.loads(progress.read_text())
        self.assertEqual(report["phase"], "complete")
        self.assertEqual(report["completed_roots"], 3)
        self.assertEqual(report["files_hashed"], 75)
        self.assertEqual(
            report["bytes_hashed"], sum(item["bytes"] for item in results[0][0])
        )
        self.assertFalse(progress.with_name(progress.name + ".tmp").exists())

    def test_no_cross_run_cache_and_worker_limits(self):
        before, _ = fingerprint_sources(self.sources, workers=4)
        victim = self.sources[0] / "samples/A/00.fasta"
        info = victim.stat()
        victim.write_bytes(b"Z" * info.st_size)
        os.utime(victim, ns=(info.st_atime_ns, info.st_mtime_ns))
        after, _ = fingerprint_sources(self.sources, workers=4)
        self.assertNotEqual(before[0]["sha256"], after[0]["sha256"])
        for workers in (0, 9, -1, True, 1.5, "4"):
            with self.subTest(workers=workers), self.assertRaises(AnnotationError):
                fingerprint_sources(self.sources, workers=workers)
