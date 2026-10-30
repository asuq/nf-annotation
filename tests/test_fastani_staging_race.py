"""Native staging must not change the QC inode while a reference scan reads it."""

import os
import threading
import unittest
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from unittest.mock import patch

import test_reference_publication as fixture
from annotation_common import AnnotationError, write_json
from build_fastani_inputs import ensure_fastani_input
from fingerprinting.source_fingerprint import checked_digest, fingerprint, signature
from verify_reference_acceptance import verify


class FastaniStagingRaceTests(unittest.TestCase):
    def setUp(self):
        self.case = fixture.ReferencePublicationTests()
        self.case.setUp()
        self.addCleanup(self.case.doCleanups)
        self.view = self.case.publish(self.case.source, "view")
        accession = self.case.accessions[0]
        self.source = self.case.source / f"samples/{accession}/staged/{accession}.fasta"
        self.acceptance = self.case.root / "acceptance.json"
        write_json(self.acceptance, {"complete": True})
        self.expected = fingerprint(self.view)

    def during_reference_scan(self, action):
        entered, release = threading.Event(), threading.Event()

        def paused(path, before):
            if path == self.source:
                entered.set()
                if not release.wait(10):
                    raise AssertionError("staging did not finish")
            return checked_digest(path, before)

        with (
            patch(
                "fingerprinting.source_fingerprint.checked_digest", side_effect=paused
            ),
            ThreadPoolExecutor(max_workers=1) as executor,
        ):
            future = executor.submit(verify, self.view, self.acceptance, workers=4)
            try:
                self.assertTrue(
                    entered.wait(10), "reference scan did not reach staged FASTA"
                )
                action()
            finally:
                release.set()
            future.result(timeout=10)

    def exercise_independent_staging(self):
        before = self.source.stat()
        incoming = self.case.root / "incoming.fasta"
        incoming.symlink_to(self.source)

        def stage():
            for name in ("cohort-a", "cohort-b"):
                directory = self.case.root / name / "fastani_inputs"
                with patch.object(
                    Path,
                    "hardlink_to",
                    side_effect=AssertionError("QC hard link created"),
                ):
                    relative, absolute = ensure_fastani_input(
                        "incoming.fasta", self.case.root, "GENOME", directory
                    )
                target = Path(absolute)
                self.assertEqual(relative, "fastani_inputs/GENOME.fasta")
                self.assertEqual(target.read_bytes(), self.source.read_bytes())
                self.assertFalse(target.is_symlink())
                self.assertFalse(target.samefile(self.source))
                # Link/unlink activity on the independent staging inode cannot touch QC.
                alias = directory / "second-link.fasta"
                alias.hardlink_to(target)
                alias.unlink()
                ensure_fastani_input(
                    "incoming.fasta", self.case.root, "GENOME", directory
                )
                target.unlink()

        self.during_reference_scan(stage)
        self.assertEqual(signature(before), signature(self.source.stat()))
        self.assertEqual(before.st_nlink, self.source.stat().st_nlink)
        self.assertEqual(self.expected, fingerprint(self.view, workers=4))

    def test_concurrent_clone_staging_and_cleanup_preserve_source(self):
        self.exercise_independent_staging()

    def test_concurrent_copy_fallback_preserves_source(self):
        with patch("independent_copy.sys.platform", "no_clone_primitive"):
            self.exercise_independent_staging()

    def test_legacy_link_creation_and_removal_reproduce_the_race(self):
        before = self.source.stat()

        def legacy_stage():
            target = self.case.root / "legacy.fasta"
            target.hardlink_to(self.source)
            target.unlink()

        with self.assertRaisesRegex(AnnotationError, "changed during fingerprinting"):
            self.during_reference_scan(legacy_stage)
        self.assertEqual(before.st_mtime_ns, self.source.stat().st_mtime_ns)
        self.assertEqual(before.st_nlink, self.source.stat().st_nlink)
        self.assertEqual(self.expected, fingerprint(self.view))

    def test_real_edit_with_restored_size_and_mtime_still_fails(self):
        before = self.source.stat()

        def mutate():
            content = self.source.read_bytes()
            changed = content.replace(b"ATGTGA", b"ATGAGA", 1)
            self.assertNotEqual(content, changed)
            self.source.write_bytes(changed)
            os.utime(self.source, ns=(before.st_atime_ns, before.st_mtime_ns))

        with self.assertRaisesRegex(AnnotationError, "changed during fingerprinting"):
            self.during_reference_scan(mutate)

    def test_identical_file_replacement_still_fails(self):
        before = self.source.stat()

        def replace():
            replacement = self.case.root / "replacement.fasta"
            replacement.write_bytes(self.source.read_bytes())
            os.utime(replacement, ns=(before.st_atime_ns, before.st_mtime_ns))
            replacement.replace(self.source)

        with self.assertRaisesRegex(AnnotationError, "changed during fingerprinting"):
            self.during_reference_scan(replace)

    def test_symlink_retarget_still_fails(self):
        def retarget():
            elsewhere = self.case.root / "elsewhere.fasta"
            elsewhere.write_bytes(self.source.read_bytes())
            self.source.unlink()
            self.source.symlink_to(elsewhere)

        with self.assertRaises((OSError, AnnotationError)):
            self.during_reference_scan(retarget)
