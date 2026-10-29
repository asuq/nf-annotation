"""Reference publications retain identities and reject changed dependencies."""

import copy
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from annotation_common import (
    AnnotationError,
    digest,
    identity,
    read_json,
    read_tsv,
    write_json,
)
from annotation_source import import_source, validate_source
from benchmark_ani_reuse import make_source
from complete_ani_references import complete
from import_ani_batch import run
from publication_references import REFERENCE_CONTRACT, REFERENCE_FILE, Publication
from fingerprinting.source_fingerprint import fingerprint, fingerprint_sources
from verify_reference_acceptance import verify


class ReferencePublicationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(dir="/tmp")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        self.source, self.accessions, _ = make_source(self.root, 2)

    def publish(self, source: Path, name: str, mode: str = "reference") -> Path:
        output = self.root / name
        output.mkdir()
        records = [
            {"meta": row, "gcode": "4", "source_index": 0}
            for row in read_tsv(source / "tables/validated_samples.tsv")
        ]
        audit = self.root / (name + ".json")
        run(
            records,
            [source],
            output,
            self.root / (name + "-channels"),
            audit,
            mode=mode,
            fingerprints=[fingerprint(source)],
        )
        shutil.copytree(source / "tables", output / "tables")
        manifest = copy.deepcopy(read_json(source / "annotation_results.json"))
        if mode == "reference":
            complete([audit], output / "tables/validated_samples.tsv", output)
            manifest["pipeline_contract"] = REFERENCE_CONTRACT
            manifest["reference_manifest_sha256"] = digest(output / REFERENCE_FILE)
        else:
            manifest["pipeline_contract"] = "nf-annotation-v0.4"
            manifest.pop("reference_manifest_sha256", None)
        for accession, bundle in manifest.get("bundles", {}).items():
            bundle["manifest_sha256"] = digest(output / f"samples/{accession}/annotation/bundle/bundle.json")
        write_json(output / "annotation_results.json", manifest)
        validate_source(output)
        return output

    def test_parallel_shared_dependencies_are_read_once_with_stable_records(self):
        first = self.publish(self.source, "first")
        second = self.publish(self.source, "second")
        roots = [second, self.source, first, second]
        expected = [fingerprint(root) for root in roots]
        records, metrics = fingerprint_sources(roots, workers=4)
        self.assertEqual(records, expected)
        self.assertEqual(metrics["unique_roots_completed"], 3)
        self.assertEqual(metrics["shared_root_hits"], 3)
        self.assertEqual(metrics["files_hashed"], sum(item["files"] for item in expected[:3]))
        self.assertLessEqual(metrics["peak_pending_reads"], 8)

    def test_reference_chain_is_flattened_and_can_be_materialised(self):
        first = self.publish(self.source, "first")
        audit = read_json(self.root / "first.json")
        self.assertTrue(
            all(
                row["files_copied"] == row["files_cloned"] == 0
                for row in audit["samples"]
            )
        )
        self.assertTrue(all(row["links_created"] > 0 for row in audit["samples"]))
        second = self.publish(first, "second")
        self.assertEqual(
            [
                source["root"]
                for source in read_json(second / REFERENCE_FILE)["sources"]
            ],
            [str(self.source)],
        )
        shutil.rmtree(first)
        validate_source(second)
        portable = self.publish(second, "portable", "self_contained")
        expected = {
            path.relative_to(portable): digest(path)
            for path in portable.rglob("*")
            if path.is_file()
        }
        shutil.rmtree(second)
        shutil.rmtree(self.source)
        validate_source(portable)
        self.assertFalse(any(path.is_symlink() for path in portable.rglob("*")))
        self.assertEqual(
            expected,
            {
                path.relative_to(portable): digest(path)
                for path in portable.rglob("*")
                if path.is_file()
            },
        )

    def test_changed_source_is_rejected_even_with_size_and_mtime_preserved(self):
        view = self.publish(self.source, "view")
        path = (
            self.source
            / f"samples/{self.accessions[0]}/staged/{self.accessions[0]}.fasta"
        )
        before = path.stat()
        content = path.read_text()
        path.write_text(content.replace("ATGTGA", "ATGAGA", 1))
        os.utime(path, ns=(before.st_atime_ns, before.st_mtime_ns))
        with self.assertRaisesRegex(AnnotationError, "Referenced source changed"):
            fingerprint(view, workers=4)
        with self.assertRaisesRegex(AnnotationError, "Referenced source changed"):
            validate_source(view)

    def test_missing_and_retargeted_dependencies_fail(self):
        view = self.publish(self.source, "view")
        link = view / f"samples/{self.accessions[0]}/staged"
        link.unlink()
        link.symlink_to(self.root, target_is_directory=True)
        with self.assertRaisesRegex(AnnotationError, "changed publication reference"):
            validate_source(view)
        link.unlink()
        link.symlink_to(
            self.source / f"samples/{self.accessions[0]}/staged",
            target_is_directory=True,
        )
        self.source.rename(self.root / "moved-source")
        with self.assertRaises((OSError, AnnotationError)):
            validate_source(view)

    def test_circular_reference_is_not_followed(self):
        view = self.publish(self.source, "view")
        record = read_json(view / REFERENCE_FILE)
        record["sources"][0]["root"] = str(view)
        record["reference_id"] = identity(
            {key: value for key, value in record.items() if key != "reference_id"}
        )
        write_json(view / REFERENCE_FILE, record)
        manifest = read_json(view / "annotation_results.json")
        manifest["reference_manifest_sha256"] = digest(view / REFERENCE_FILE)
        write_json(view / "annotation_results.json", manifest)
        with self.assertRaisesRegex(AnnotationError, "circular"):
            Publication(view)
        with self.assertRaisesRegex(AnnotationError, "circular"):
            fingerprint_sources([view], workers=4)

    def test_final_acceptance_detects_input_metadata_and_source_changes(self):
        view = self.publish(self.source, "view")
        acceptance = self.root / "acceptance.json"
        write_json(acceptance, {"complete": True})
        snapshot = self.root / "snapshot.json"
        write_json(snapshot, fingerprint(view))
        verify(view, acceptance, snapshot)
        (view / "tables/master_table.tsv").write_text(
            (view / "tables/master_table.tsv").read_text() + "\n"
        )
        with self.assertRaisesRegex(AnnotationError, "changed during annotation"):
            verify(view, acceptance, snapshot)
        (self.source / f"samples/{self.accessions[0]}/prokka/prokka.log").write_text(
            "exit_code=1\n"
        )
        with self.assertRaisesRegex(AnnotationError, "changed during analysis"):
            verify(view, acceptance)

    def test_locally_rebuilt_bundle_remains_an_explicit_dependency(self):
        accession = self.accessions[0]
        path = self.source / f"samples/{accession}/annotation/bundle/bundle.json"
        data = read_json(path)
        data["producer_sha256"] = "0" * 64
        write_json(path, data)
        manifest = read_json(self.source / "annotation_results.json")
        manifest["bundles"][accession]["manifest_sha256"] = digest(path)
        write_json(self.source / "annotation_results.json", manifest)
        first = self.publish(self.source, "rebuilt-view")
        self.assertEqual(
            read_json(first / REFERENCE_FILE)["samples"][accession]["bundle"],
            {"local": True},
        )
        second = self.publish(first, "linked-rebuilt-view")
        self.assertIn(
            str(first),
            [
                record["root"]
                for record in read_json(second / REFERENCE_FILE)["sources"]
            ],
        )
        shutil.rmtree(first)
        with self.assertRaises((OSError, AnnotationError)):
            validate_source(second)

    def test_reannotation_import_materialises_declared_links(self):
        view = self.publish(self.source, "view")
        imported = self.root / "reannotation-import"
        import_source(view, imported)
        self.assertTrue((imported / "source_snapshot.json").is_file())
        self.assertFalse(any(path.is_symlink() for path in imported.rglob("*")))
        for accession in self.accessions:
            source_bundle = self.source / f"samples/{accession}/annotation/bundle"
            target_bundle = (
                imported / f"inherited/samples/{accession}/annotation/bundle"
            )
            self.assertEqual(
                digest(source_bundle / "bundle.json"),
                digest(target_bundle / "bundle.json"),
            )
