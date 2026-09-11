"""Portable-source integrity checks for native v0.4 reannotation."""

import shutil
import unittest

import test_annotation_bundle as bundle_fixture
from annotation_common import AnnotationError, digest, read_json, write_json, write_tsv
from annotation_source import import_source, validate_source
from annotation_workflow_fixture import publish_disabled


class AnnotationSourceTests(unittest.TestCase):
    def setUp(self):
        fixture = bundle_fixture.AnnotationBundleTests()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        self.root = fixture.root
        self.source = self.root / "source"
        (self.source / "tables").mkdir(parents=True)
        self.bundle = self.source / "samples/A/annotation/bundle"
        shutil.copytree(fixture.build("A"), self.bundle)
        write_tsv(
            self.source / "tables/master_table.tsv",
            ("accession", "Gcode"),
            [dict(accession="A", Gcode=4)],
        )
        write_tsv(
            self.source / "tables/sample_status.tsv",
            ("accession",),
            [dict(accession="A")],
        )
        write_tsv(
            self.source / "tables/validated_samples.tsv",
            ("accession",),
            [dict(accession="A")],
        )
        publish_disabled(self.source, ["A"], [self.bundle])

    def test_native_source_import_preserves_bundle_and_rejects_changed_tables(self):
        self.assertTrue(validate_source(self.source)["complete"])
        output = self.root / "imported"
        import_source(self.source, output)
        self.assertEqual(
            (self.bundle / "bundle.json").read_bytes(),
            (output / "inherited/samples/A/annotation/bundle/bundle.json").read_bytes(),
        )
        (self.source / "tables/protein_function_summary.tsv").write_text("changed\n")
        with self.assertRaisesRegex(AnnotationError, "table changed"):
            validate_source(self.source)

    def test_linked_bundle_manifest_and_upstream_evidence_are_rejected(self):
        original = self.root / "bundle.json"
        shutil.move(self.bundle / "bundle.json", original)
        (self.bundle / "bundle.json").symlink_to(original)
        with self.assertRaisesRegex(AnnotationError, "symlink"):
            validate_source(self.source)
        (self.bundle / "bundle.json").unlink()
        shutil.move(original, self.bundle / "bundle.json")
        (self.source / "samples/A/work_reference").symlink_to(self.root / "absent_work")
        with self.assertRaisesRegex(AnnotationError, "symlink"):
            import_source(self.source, self.root / "invalid_import")
        self.assertFalse((self.root / "invalid_import").exists())

    def test_detailed_evidence_tables_must_remain_in_the_manifest(self):
        path = self.source / "annotation_results.json"
        manifest = read_json(path)
        for name in (
            "eggnog_seed_hits",
            "eggnog_annotations",
            "cog_assignments",
            "pfam_domains",
            "kofam_hits",
            "defence_systems",
            "defence_genes",
        ):
            with self.subTest(name=name):
                table = f"tables/{name}.tsv"
                checksum = manifest["tables"].pop(table)
                write_json(path, manifest)
                with self.assertRaisesRegex(AnnotationError, "omits required tables"):
                    validate_source(self.source)
                manifest["tables"][table] = checksum

    def test_invalid_failed_bundle_and_undeclared_accessions_are_rejected(self):
        bundle = read_json(self.bundle / "bundle.json")
        bundle.update(status="upstream_failed", accession="wrong")
        write_json(self.bundle / "bundle.json", bundle)
        manifest = read_json(self.source / "annotation_results.json")
        manifest["bundles"]["A"]["manifest_sha256"] = digest(
            self.bundle / "bundle.json"
        )
        write_json(self.source / "annotation_results.json", manifest)
        with self.assertRaisesRegex(AnnotationError, "Invalid published bundle"):
            validate_source(self.source)
        manifest["results"]["undeclared"] = {}
        write_json(self.source / "annotation_results.json", manifest)
        with self.assertRaisesRegex(AnnotationError, "undeclared accession"):
            validate_source(self.source)
