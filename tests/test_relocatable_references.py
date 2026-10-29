"""Relative references survive a collection move without repairing publications."""

import copy
import os
import shutil
import tempfile
import unittest
from pathlib import Path

import test_reference_publication as fixture
from annotation_common import AnnotationError, digest, identity, read_json, write_json
from annotation_source import import_source, validate_source
from fingerprinting.source_fingerprint import fingerprint
from publication_references import REFERENCE_CONTRACT, REFERENCE_FILE, Publication


class RelocatableReferenceTests(unittest.TestCase):
    def setUp(self):
        self.case = fixture.ReferencePublicationTests()
        self.case.setUp()
        self.addCleanup(self.case.doCleanups)

    def reseal(self, view, record):
        record["reference_id"] = identity(
            {key: value for key, value in record.items() if key != "reference_id"}
        )
        write_json(view / REFERENCE_FILE, record)
        manifest = read_json(view / "annotation_results.json")
        manifest["reference_manifest_sha256"] = digest(view / REFERENCE_FILE)
        write_json(view / "annotation_results.json", manifest)

    def test_move_collection_then_reuse_materialise_and_import(self):
        view = self.case.publish(self.case.source, "view")
        before = fingerprint(view, workers=1)
        links = {
            p.relative_to(view).as_posix(): os.readlink(p)
            for p in (view / "samples").rglob("*")
            if p.is_symlink()
        }
        self.assertTrue(links)
        self.assertTrue(all(not Path(value).is_absolute() for value in links.values()))
        document = (view / REFERENCE_FILE).read_bytes()
        with tempfile.TemporaryDirectory(dir="/tmp") as relocated:
            old = self.case.root
            moved = Path(relocated).resolve() / "different parent" / "collection"
            moved.parent.mkdir()
            old.rename(moved)
            self.assertFalse(old.exists())
            view = moved / "view"
            self.case.root = moved
            self.case.source = moved / "source"
            validate_source(view)
            self.assertEqual(before["sha256"], fingerprint(view, workers=4)["sha256"])
            self.assertEqual(document, (view / REFERENCE_FILE).read_bytes())
            self.assertEqual(
                links,
                {
                    p.relative_to(view).as_posix(): os.readlink(p)
                    for p in (view / "samples").rglob("*")
                    if p.is_symlink()
                },
            )
            successor = self.case.publish(view, "successor")
            shutil.rmtree(view)
            validate_source(successor)
            portable = self.case.publish(successor, "materialised", "self_contained")
            imported = moved / "annotation-import"
            import_source(successor, imported)
            shutil.rmtree(successor)
            shutil.rmtree(self.case.source)
            validate_source(portable)
            self.assertFalse(any(p.is_symlink() for p in imported.rglob("*")))

    def test_escaped_absolute_and_aliased_provider_locations_are_rejected(self):
        view = self.case.publish(self.case.source, "view")
        original = read_json(view / REFERENCE_FILE)
        for location in ("../source", str(self.case.source), "nested/../source"):
            with self.subTest(location=location):
                record = copy.deepcopy(original)
                record["sources"][0]["path"] = location
                self.reseal(view, record)
                with self.assertRaises((AnnotationError, OSError)):
                    Publication(view)
        record = copy.deepcopy(original)
        record["qc_publication_root"] = str(self.case.source)
        self.reseal(view, record)
        with self.assertRaises(AnnotationError):
            Publication(view)
        record = copy.deepcopy(original)
        (self.case.source / "alias").symlink_to(
            self.case.source, target_is_directory=True
        )
        record["sources"][0]["path"] = "alias"
        self.reseal(view, record)
        with self.assertRaisesRegex(AnnotationError, "symlink"):
            Publication(view)

    def test_final_gate_rechecks_transient_view_inputs_without_retaining_them(self):
        from verify_reference_acceptance import verify

        first = self.case.publish(self.case.source, "first")
        second = self.case.publish(first, "second")
        checks = self.case.root / "source_validation.json"
        write_json(checks, {"fingerprints": [fingerprint(first)]})
        acceptance = self.case.root / "acceptance.json"
        write_json(acceptance, {"complete": True})
        for workers in (1, 4, 8):
            verify(second, acceptance, workers=workers, source_validation=checks)
        table = first / "tables/master_table.tsv"
        table.write_text(table.read_text() + "\n")
        with self.assertRaisesRegex(AnnotationError, "Reference input changed"):
            verify(second, acceptance, workers=4, source_validation=checks)
        shutil.rmtree(first)
        validate_source(second)

    def test_reference_graph_cycle_fails_before_following_payloads(self):
        from publication_references import REFERENCE_V2_CONTRACT

        view = self.case.publish(self.case.source, "view")
        record = copy.deepcopy(read_json(view / REFERENCE_FILE))
        record["qc_publication_root"] = "../view"
        record["sources"] = [{"path": ".", "sha256": "0" * 64}]
        self.reseal(self.case.source, record)
        manifest = read_json(self.case.source / "annotation_results.json")
        manifest["pipeline_contract"] = REFERENCE_V2_CONTRACT
        write_json(self.case.source / "annotation_results.json", manifest)
        with self.assertRaisesRegex(AnnotationError, "circular"):
            fingerprint(view, workers=4)

    def test_dangling_reference_descriptor_cannot_masquerade_as_qc(self):
        (self.case.source / REFERENCE_FILE).symlink_to(self.case.root / "missing")
        with self.assertRaises(AnnotationError):
            Publication(self.case.source)
        with self.assertRaises(AnnotationError):
            fingerprint(self.case.source, workers=4)

    def test_unknown_v2_fields_and_nonboolean_local_bundle_fail(self):
        view = self.case.publish(self.case.source, "view")
        original = read_json(view / REFERENCE_FILE)
        altered = copy.deepcopy(original)
        altered["legacy_root"] = str(self.case.source)
        self.reseal(view, altered)
        with self.assertRaisesRegex(AnnotationError, "Unknown v2"):
            Publication(view)
        altered = copy.deepcopy(original)
        altered["samples"][self.case.accessions[0]]["bundle"] = {"local": 1}
        self.reseal(view, altered)
        with self.assertRaisesRegex(AnnotationError, "Invalid bundle"):
            Publication(view)

    def test_wrong_boundary_and_derived_provider_are_rejected(self):
        from annotation_common import read_tsv
        from import_ani_batch import run

        records = [
            {"meta": row, "gcode": "4", "source_index": 0}
            for row in read_tsv(self.case.source / "tables/validated_samples.tsv")
        ]
        target = self.case.root / "wrong-boundary-output"
        target.mkdir()
        boundary = self.case.root / "other-qc"
        boundary.mkdir()
        with self.assertRaises(AnnotationError):
            run(
                records,
                [self.case.source],
                target,
                self.case.root / "channels",
                self.case.root / "audit.json",
                mode="reference",
                qc_root=boundary,
                fingerprints=[fingerprint(self.case.source)],
            )
        (self.case.source / "cohort/fastani").mkdir(parents=True)
        with self.assertRaises(AnnotationError):
            self.case.publish(self.case.source, "derived")

    def test_legacy_absolute_view_remains_readable_but_is_not_a_new_provider(self):
        view = self.case.publish(self.case.source, "legacy")
        document = read_json(view / REFERENCE_FILE)
        document.pop("qc_publication_root")
        document["schema_version"] = 1
        document["sources"] = [
            {
                "root": str(self.case.source),
                "sha256": fingerprint(self.case.source)["sha256"],
            }
        ]
        for link in list((view / "samples").rglob("*")):
            if link.is_symlink():
                target = link.resolve()
                link.unlink()
                link.symlink_to(target, target_is_directory=target.is_dir())
        self.reseal(view, document)
        manifest = read_json(view / "annotation_results.json")
        manifest["pipeline_contract"] = REFERENCE_CONTRACT
        write_json(view / "annotation_results.json", manifest)
        validate_source(view)
        with self.assertRaises(AnnotationError):
            self.case.publish(view, "not-migrated")
