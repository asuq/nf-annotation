"""SIF receipt reuse does not read image payloads during runtime preflight."""

import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from annotation_tasks import runtime_identity
from prepare_annotation_runtime import prepare, read_receipt


class RuntimeReceiptTests(unittest.TestCase):
    def test_preparation_reuse_and_changed_image(self):
        with tempfile.TemporaryDirectory() as directory:
            image = Path(directory) / "tool.sif"
            image.write_bytes(b"frozen image")
            with self.assertRaisesRegex(ValueError, "Prepare the SIF receipt"):
                runtime_identity(str(image))
            receipt = prepare(image)
            with patch(
                "hashlib.file_digest", side_effect=AssertionError("image reread")
            ):
                self.assertEqual(
                    runtime_identity(str(image)), "sif-md5:" + receipt["md5"]
                )
            os.utime(
                image,
                ns=(image.stat().st_atime_ns, image.stat().st_mtime_ns + 1000000000),
            )
            with self.assertRaisesRegex(ValueError, "changed image"):
                read_receipt(image)

    def test_preflight_never_hashes_database_or_sif_payloads(self):
        import test_annotation_resources as fixture_module
        from annotation_tasks import preflight
        from prepare_annotation_runtime import prepare

        fixture = fixture_module.ResourceTests()
        fixture.setUp()
        self.addCleanup(fixture.tearDown)
        self.root = fixture.root
        source = fixture.prepared()
        image = self.root / "helper.sif"
        image.write_bytes(b"runtime fixture")
        prepare(image)
        original_open = Path.open

        def guarded_open(path, *args, **kwargs):
            if path in (source / "resource.dat", image):
                raise AssertionError("Runtime preflight reread a payload")
            return original_open(path, *args, **kwargs)

        config = {
            "annotation_tools": "pfam",
            "helper_container": str(image),
            "tools": {
                "pfam": {
                    "container": str(image),
                    "resource": str(source),
                    "cpus": 1,
                    "memory_gib": 32,
                }
            },
        }
        original_stat = Path.stat

        def guarded_stat(path, *args, **kwargs):
            if path == source / "resource.dat":
                raise AssertionError("Runtime preflight traversed the payload tree")
            return original_stat(path, *args, **kwargs)

        with (
            patch.object(Path, "open", guarded_open),
            patch.object(Path, "stat", guarded_stat),
        ):
            result = preflight(config)
            self.assertEqual(result["enabled_tools"], ["pfam"])
