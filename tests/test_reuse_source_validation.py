"""Validate source receipt identity and detect nested edits independently of stat caches."""

import copy
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from annotation_common import AnnotationError, identity
from reuse_source_validation import checked_roots, fingerprint, validate


class ReuseSourceValidationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(dir="/tmp")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        (self.root / "tables").mkdir()
        (self.root / "samples/A").mkdir(parents=True)
        (self.root / "annotation_results.json").write_text("{}\n")
        self.data = self.root / "samples/A/sequence.fasta"
        self.data.write_text(">A\nACGT\n")

    def test_nested_content_changes_invalidate_even_with_preserved_size_and_mtime(self):
        first = fingerprint(self.root)
        before = self.data.stat()
        self.data.write_text(">A\nTGCA\n")
        os.utime(self.data, ns=(before.st_atime_ns, before.st_mtime_ns))
        changed = fingerprint(self.root)
        self.assertEqual(first["bytes"], changed["bytes"])
        self.assertNotEqual(first["sha256"], changed["sha256"])
        (self.root / "pipeline_info").mkdir()
        (self.root / "pipeline_info/trace.tsv").write_text("execution-only\n")
        self.assertEqual(changed, fingerprint(self.root))

    def test_receipts_are_bound_to_content_code_and_complete_root_membership(self):
        records = [fingerprint(self.root)]
        with patch("reuse_source_validation.validate_source") as checker:
            receipts = validate([self.root], records)
        checker.assert_called_once_with(self.root)
        self.assertEqual(checked_roots(receipts, records), {self.root})
        bad = copy.deepcopy(receipts)
        bad[0]["validator"] = "old-validator"
        bad[0]["receipt_id"] = identity(
            {key: value for key, value in bad[0].items() if key != "receipt_id"}
        )
        with self.assertRaises(AnnotationError):
            checked_roots(bad, records)
        self.data.write_text("changed bytes\n")
        with self.assertRaises(AnnotationError):
            checked_roots(receipts, [fingerprint(self.root)])
        with self.assertRaises(AnnotationError):
            checked_roots(receipts * 2, records)

    def test_fingerprint_rejects_links_in_scientific_payload(self):
        (self.root / "samples/A/link").symlink_to(self.data)
        with self.assertRaisesRegex(AnnotationError, "symlink"):
            fingerprint(self.root)
