"""Independent publication, exact bundle reuse and per-sample batch failures."""

import shutil
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import test_annotation_bundle as fixture
from annotation_common import (
    AnnotationError,
    bundle_proteins,
    digest,
    read_json,
    write_json,
)
from import_ani_batch import copy_independent, import_sample, regular_files, run


class AniBatchImportTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixture.AnnotationBundleTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root
        self.publication = self.root / "source"
        self.source = self.publication / "samples/A"
        self.destination = self.root / "destination"
        self.destination.mkdir()
        self.channels = self.root / "channels"
        self.meta = {"accession": "A", "internal_id": "A"}
        for original, name in (
            (self.fixture.genome, "staged/A.fasta"),
            (self.fixture.faa, "prokka/prokka.faa"),
            (self.fixture.gff, "prokka/prokka.gff"),
            (self.fixture.gbk, "prokka/prokka.gbk"),
        ):
            target = self.source / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(original, target)
        for name in (
            "16s/best_16S.fna",
            "16s/16S_status.tsv",
            "checkm2/checkm2_summary.tsv",
        ):
            target = self.source / name
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text("fixture\n")
        (self.source / "prokka/prokka.log").write_text("exit_code=0\n")
        self.bundle = self.source / "annotation/bundle"
        shutil.copytree(self.fixture.build("A"), self.bundle)
        record = read_json(self.bundle / "bundle.json")
        record["prokka_log_sha256"] = digest(self.source / "prokka/prokka.log")
        write_json(self.bundle / "bundle.json", record)

    def test_reuse_keeps_native_bundle_bytes_and_publication_is_independent(self):
        expected = {
            name: digest(path) for name, path in regular_files(self.source).items()
        }
        with patch(
            "import_ani_batch.build_bundle",
            side_effect=AssertionError("unnecessary bundle rebuild"),
        ):
            result = import_sample(
                self.source, self.meta, "4", self.destination, self.channels
            )
        target = self.destination / "samples/A"
        self.assertEqual(result["bundle"], "reused")
        self.assertEqual(result["files_copied"] + result["files_cloned"], len(expected))
        self.assertEqual(
            {name: digest(path) for name, path in regular_files(target).items()},
            expected,
        )
        self.assertNotEqual(
            (target / "prokka/prokka.faa").stat().st_ino, self.fixture.faa.stat().st_ino
        )
        self.assertEqual(
            bundle_proteins(target / "annotation/bundle"), bundle_proteins(self.bundle)
        )
        again = import_sample(
            self.source, self.meta, "4", self.destination, self.channels
        )
        self.assertEqual(again["publication"], "verified_existing")
        self.assertEqual(again["files_copied"] + again["files_cloned"], 0)
        (self.source / "prokka/prokka.faa").write_text("changed\n")
        self.assertEqual(
            digest(target / "prokka/prokka.faa"), expected["prokka/prokka.faa"]
        )
        shutil.rmtree(self.source)
        bundle_proteins(target / "annotation/bundle")

    def test_old_producer_is_rebuilt_with_the_unchanged_adapter(self):
        record = read_json(self.bundle / "bundle.json")
        record["producer_sha256"] = "0" * 64
        write_json(self.bundle / "bundle.json", record)
        result = import_sample(
            self.source, self.meta, "4", self.destination, self.channels
        )
        self.assertEqual(result["bundle"], "rebuilt")
        bundle_proteins(self.destination / "samples/A/annotation/bundle")

    def test_copy_fallback_preserves_bytes_without_hardlinks(self):
        target = self.root / "copied"
        with patch("import_ani_batch.sys.platform", "unsupported"):
            self.assertEqual(copy_independent(self.fixture.faa, target), "copied")
        self.assertEqual(target.read_bytes(), self.fixture.faa.read_bytes())
        self.assertNotEqual(target.stat().st_ino, self.fixture.faa.stat().st_ino)

    def test_changed_destination_and_source_links_fail_with_sample_identity(self):
        import_sample(self.source, self.meta, "4", self.destination, self.channels)
        (self.destination / "samples/A/prokka/prokka.faa").write_text("changed\n")
        with self.assertRaisesRegex(
            AnnotationError, "Existing publication differs for A"
        ):
            import_sample(self.source, self.meta, "4", self.destination, self.channels)
        (self.source / "escape").symlink_to(self.fixture.genome)
        audit = self.root / "audit.json"
        with self.assertRaisesRegex(AnnotationError, "failed samples"):
            run(
                [{"meta": self.meta, "gcode": "4", "source_index": 0}],
                [self.publication],
                self.destination,
                self.channels,
                audit,
            )
        report = read_json(audit)
        self.assertFalse(report["complete"])
        self.assertEqual(report["samples"][0]["accession"], "A")
        self.assertEqual(report["samples"][0]["status"], "failed")

    def test_no_code_omits_bundle_and_foreign_samples_directory_is_rejected(self):
        result = import_sample(
            self.source, self.meta, "NA", self.destination, self.channels
        )
        self.assertEqual(result["bundle"], "not_requested")
        self.assertFalse((self.destination / "samples/A/annotation").exists())
        other = self.root / "other-destination"
        other.mkdir()
        (other / "samples").symlink_to(self.source, target_is_directory=True)
        with self.assertRaisesRegex(AnnotationError, "ordinary directory"):
            import_sample(self.source, self.meta, "4", other, self.channels)
