"""Portable-source integrity checks for native v0.4 reannotation."""

import csv
import json
import shutil
import unittest

import test_annotation_bundle as bundle_fixture
from annotation_common import (
    AnnotationError,
    digest,
    read_json,
    read_tsv,
    write_json,
    write_tsv,
)
from annotation_source import import_source, validate_source
from annotation_workflow_fixture import publish_disabled, refresh_tables


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
            [{"accession": "A", "Gcode": 4}],
        )
        write_tsv(
            self.source / "tables/sample_status.tsv",
            ("accession",),
            [{"accession": "A"}],
        )
        write_tsv(
            self.source / "tables/validated_samples.tsv",
            ("accession",),
            [{"accession": "A"}],
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

    def test_large_validation_cells_roundtrip_and_restore_csv_limit(self):
        previous_limit = csv.field_size_limit()
        errors = [
            {"query": f"protein_{index}", "field": "categories", "value": "COG0123"}
            for index in range(3000)
        ]
        expected = [{
            "accession": "A",
            "field_errors": json.dumps(errors),
            "note": 'Quoted "native" evidence\tand a\nsecond line',
        }]
        self.assertGreater(len(expected[0]["field_errors"]), previous_limit)
        path = self.root / "large_status.tsv"
        write_tsv(path, expected[0], expected)
        self.assertEqual(read_tsv(path), expected)
        self.assertEqual(csv.field_size_limit(), previous_limit)
        for content, message in (
            (
                "accession\taccession\nA\t" + expected[0]["field_errors"] + "\n",
                "header",
            ),
            (
                "accession\tfield_errors\n" + "A" * (previous_limit + 1) + "\n",
                "Malformed table row",
            ),
        ):
            with self.subTest(message=message):
                path.write_text(content)
                with self.assertRaisesRegex(AnnotationError, message):
                    read_tsv(path)
                self.assertEqual(csv.field_size_limit(), previous_limit)

    def test_portable_source_retains_large_metadata_cells(self):
        path = self.source / "tables/master_table.tsv"
        rows = read_tsv(path)
        rows[0]["analysis_note"] = "native evidence " * 10000
        self.assertGreater(len(rows[0]["analysis_note"]), csv.field_size_limit())
        write_tsv(path, rows[0], rows)
        refresh_tables(self.source)
        self.assertTrue(validate_source(self.source)["complete"])
        output = self.root / "large_import"
        import_source(self.source, output)
        self.assertEqual(
            path.read_bytes(), (output / "upstream_master.tsv").read_bytes()
        )

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
