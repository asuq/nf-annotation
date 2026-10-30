"""Shared providers get one semantic validation per fresh source-group pass."""

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import annotation_source
import reuse_source_validation
import test_reference_publication as fixture
from annotation_common import AnnotationError
from fingerprinting.source_fingerprint import fingerprint_sources

ROOT = Path(__file__).resolve().parents[1]


class ValidationSharingTests(unittest.TestCase):
    def test_provider_validation_is_shared_but_views_are_still_checked(self):
        case = fixture.ReferencePublicationTests()
        case.setUp()
        self.addCleanup(case.doCleanups)
        first = case.publish(case.source, "first")
        second = case.publish(case.source, "second")
        roots = [first, second]
        records, _ = fingerprint_sources(roots, workers=2)
        original = annotation_source.source_batches
        with patch("annotation_source.source_batches", wraps=original) as batches:
            reuse_source_validation.validate(roots, records)
            checked = [call.args[0].resolve() for call in batches.call_args_list]
            self.assertEqual(checked.count(case.source), 1)
            self.assertEqual(checked.count(first), 1)
            self.assertEqual(checked.count(second), 1)
            reuse_source_validation.validate(roots, records)
            checked = [call.args[0].resolve() for call in batches.call_args_list]
            self.assertEqual(checked.count(case.source), 2)  # No cross-run trust cache.
        table = second / "tables/master_table.tsv"
        table.write_text(table.read_text() + "\n")
        records, _ = fingerprint_sources(roots, workers=2)
        with self.assertRaises(AnnotationError):
            reuse_source_validation.validate(roots, records)

    @unittest.skipUnless(
        shutil.which("nextflow"), "Nextflow required for code-stamp agreement"
    )
    def test_code_identity_covers_nested_helpers_and_matches_nextflow(self):
        with tempfile.TemporaryDirectory(dir="/tmp") as directory:
            root = Path(directory)
            (root / "bin/fingerprinting").mkdir(parents=True)
            (root / "bin/fingerprinting.py").write_text("first\n")
            helper = root / "bin/fingerprinting/reader.py"
            helper.write_text("nested\n")
            with patch.object(
                reuse_source_validation,
                "__file__",
                str(root / "bin/reuse_source_validation.py"),
            ):
                first = reuse_source_validation.code_id()
                helper.write_text("changed nested\n")
                second = reuse_source_validation.code_id()
            self.assertNotEqual(first, second)
            (root / "lib").mkdir()
            shutil.copy2(
                ROOT / "lib/SamplePolicy.groovy", root / "lib/SamplePolicy.groovy"
            )
            (root / "main.nf").write_text(
                "workflow { assert SamplePolicy.reuseCodeHash(params.project) == params.expected }\n"
            )
            result = subprocess.run(
                [
                    shutil.which("nextflow"),
                    "run",
                    str(root / "main.nf"),
                    "--project",
                    str(root),
                    "--expected",
                    second,
                ],
                cwd=root,
                capture_output=True,
                text=True,
                check=False,
                timeout=120,
                env=dict(
                    os.environ, NXF_DISABLE_CHECK_LATEST="true", NXF_ANSI_LOG="false"
                ),
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
