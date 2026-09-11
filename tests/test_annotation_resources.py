"""Checks for immutable annotation resource identities and profile selection."""

from __future__ import annotations

import gzip
import io
import json
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import annotation_resources as resources
import prepare_runtime_databases as preparation


class ResourceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def prepared(self, name="source"):
        root = self.root / name
        root.mkdir()
        (root / "resource.dat").write_text("fixed resource content\n")
        contract = {
            "component": "pfam",
            "version": "38.2",
            "settings": {"synthetic_fixture": True},
            "files": resources.file_records(root),
        }
        record = {
            "schema_version": 1,
            "component": "pfam",
            "resource_id": resources.identity(contract),
            "contract": contract,
            "acquisition": {"files": [{"file": "resource.dat.gz", "sha256": "0" * 64}]},
        }
        (root / resources.RESOURCE_FILE).write_text(json.dumps(record))
        return root

    def test_prepared_source_copy_retains_identity_and_reuse(self):
        source = self.prepared()
        destination = self.root / "destination"
        first = resources.prepare_resource("pfam", "38.2", source, destination)
        self.assertEqual(first, resources.validate_resource(destination, "pfam"))
        self.assertEqual(
            first, resources.prepare_resource("pfam", "38.2", source, destination)
        )
        self.assertEqual(
            first["resource_id"],
            resources.validate_resource(source, "pfam")["resource_id"],
        )

    def test_tampered_content_or_contract_fails(self):
        source = self.prepared()
        (source / "resource.dat").write_text("modified resource data\n")
        with self.assertRaises(resources.AnnotationResourceError):
            resources.validate_resource(source, "pfam")
        record = json.loads((source / resources.RESOURCE_FILE).read_text())
        record["contract"]["version"] = "other"
        (source / resources.RESOURCE_FILE).write_text(json.dumps(record))
        with self.assertRaisesRegex(
            resources.AnnotationResourceError, "identity mismatch"
        ):
            resources.validate_resource(source, "pfam")

    def test_wrong_component_and_unprepared_sources_fail(self):
        source = self.prepared()
        with self.assertRaises(resources.AnnotationResourceError):
            resources.validate_resource(source, "eggnog")
        (source / resources.RESOURCE_FILE).unlink()
        with self.assertRaises(resources.AnnotationResourceError):
            resources.prepare_resource("pfam", "38.2", source, self.root / "out")

    def test_same_version_with_different_source_is_not_reused(self):
        source = self.prepared()
        destination = self.root / "destination"
        resources.prepare_resource("pfam", "38.2", source, destination)
        record = json.loads((source / resources.RESOURCE_FILE).read_text())
        record["acquisition"]["files"][0]["sha256"] = "1" * 64
        (source / resources.RESOURCE_FILE).write_text(json.dumps(record))
        with self.assertRaisesRegex(
            resources.AnnotationResourceError, "different source content"
        ):
            resources.prepare_resource("pfam", "38.2", source, destination)
        with self.assertRaises(resources.AnnotationResourceError):
            resources.prepare_resource("pfam", "future", destination, destination)

    def test_path_traversal_inventory_fails_even_with_recomputed_identity(self):
        source = self.prepared()
        record = json.loads((source / resources.RESOURCE_FILE).read_text())
        record["contract"]["files"][0]["path"] = "../outside"
        record["resource_id"] = resources.identity(record["contract"])
        (source / resources.RESOURCE_FILE).write_text(json.dumps(record))
        with self.assertRaisesRegex(
            resources.AnnotationResourceError, "Invalid or duplicate"
        ):
            resources.validate_resource(source, "pfam")

    def test_unsafe_archive_is_rejected_before_extraction(self):
        source = self.root / "unsafe.tar.gz"
        with tarfile.open(source, "w:gz") as archive:
            info = tarfile.TarInfo("../outside")
            info.size = 3
            archive.addfile(info, io.BytesIO(b"bad"))
        with self.assertRaisesRegex(resources.AnnotationResourceError, "Unsafe"):
            resources.extract(source, self.root / "out")
        self.assertFalse((self.root / "outside").exists())

    def test_canonical_cli_accepts_prepared_annotation_source(self):
        source = self.prepared()
        destination = self.root / "pfam"
        args = preparation.parse_args(
            ["--pfam-source", str(source), "--pfam-dest", str(destination)]
        )
        records = preparation.build_component_records(args)
        self.assertEqual(records[0].component, "pfam")
        self.assertEqual(records[0].status, "prepared")
        self.assertIn(
            ("--pfam_db", str(destination.resolve())),
            preparation.build_nextflow_arguments(records),
        )
        self.assertTrue((destination / preparation.MARKER_FILE_NAME).is_file())

    def kofam_source(self, threshold="20.0", missing=False):
        source = self.root / "acquired"
        source.mkdir()
        files = {
            "profiles/prokaryote.hal": b"K00001.hmm\n",
            "profiles/K00001.hmm": b"HMMER3/f\nNAME K00001\n//\n",
            "profiles/K00002.hmm": b"HMMER3/f\nNAME K00002\n//\n",
            "profiles/eukaryote.hal": b"K00002.hmm\n",
        }
        with tarfile.open(source / "profiles.tar.gz", "w:gz") as archive:
            for name, data in files.items():
                if missing and name.endswith("K00001.hmm"):
                    continue
                info = tarfile.TarInfo(name)
                info.size = len(data)
                archive.addfile(info, io.BytesIO(data))
        with gzip.open(source / "ko_list.gz", "wt") as handle:
            handle.write(
                f"knum\tthreshold\tscore_type\tdefinition\nK00001\t{threshold}\tfull\tsynthetic\nK00002\t-\t-\tsynthetic\n"
            )
        return source

    def test_kofam_uses_exact_native_prokaryotic_subset(self):
        source = self.kofam_source()
        destination = self.root / "prepared"
        destination.mkdir()
        with patch.object(resources, "validate_hmm_profiles") as validation:
            settings = resources.prepare_kofam(source, destination)
        validation.assert_called_once_with(
            destination, [destination / "profiles/K00001.hmm"]
        )
        self.assertEqual(settings["profile_count"], 1)
        self.assertEqual(
            sorted(path.name for path in (destination / "profiles").iterdir()),
            ["K00001.hmm", "prokaryote.hal"],
        )
        self.assertIn(
            "\t20.0\tfull\n", (destination / "prokaryote_profiles.tsv").read_text()
        )

    def test_kofam_missing_profile_and_nonfinite_threshold_fail(self):
        for index, (threshold, missing) in enumerate((("nan", False), ("20", True))):
            with self.subTest(threshold=threshold, missing=missing):
                self.root = Path(self.temporary.name) / str(index)
                self.root.mkdir()
                source = self.kofam_source(threshold, missing)
                destination = self.root / "prepared"
                destination.mkdir()
                with self.assertRaises(resources.AnnotationResourceError):
                    resources.prepare_kofam(source, destination)


if __name__ == "__main__":
    unittest.main()
