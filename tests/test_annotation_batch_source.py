"""Portable source manifests must retain exact and complete shared batch evidence."""

from __future__ import annotations

import copy
import shutil
import unittest
from unittest.mock import patch

from Bio import SeqIO

import test_annotation_batch_results as fixture
from annotation_common import (
    AnnotationError,
    digest,
    identity,
    read_json,
    read_tsv,
    write_json,
    write_tsv,
)
from annotation_result import inventory, validate_native_batch, validate_raw_evidence
from annotation_source import import_source, validate_source
from annotation_workflow_fixture import publish_disabled, refresh_tables
from eggnog_batches import prepare_batches


class AnnotationBatchSourceTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixture.AnnotationBatchResultTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root
        self.source = self.root / "published"
        self.native = self.fixture.native
        self.packet = self.fixture.native_record
        self.batch_id = self.fixture.batch["batch_id"]
        self.accessions = [
            member["accession"] for member in self.fixture.batch["members"]
        ]
        self.results = {
            record["accession"]: (path, copy.deepcopy(record))
            for path, record in self.fixture.results
        }
        bundles = []
        available = {
            read_json(path / "bundle.json")["accession"]: path
            for path in self.fixture.fixture.bundles
        }
        for accession in self.accessions:
            bundle = self.source / "samples" / accession / "annotation/bundle"
            shutil.copytree(available[accession], bundle)
            bundles.append(bundle)
        (self.source / "tables").mkdir()
        write_tsv(
            self.source / "tables/master_table.tsv",
            ("accession", "Gcode"),
            [dict(accession=accession, Gcode=4) for accession in self.accessions],
        )
        for name in ("sample_status", "validated_samples"):
            write_tsv(
                self.source / "tables" / (name + ".tsv"),
                ("accession",),
                [dict(accession=accession) for accession in self.accessions],
            )
        publish_disabled(self.source, self.accessions, bundles)
        self.manifest = read_json(self.source / "annotation_results.json")
        self.manifest.update(
            enabled_tools=["eggnog"],
            native_batches={
                self.batch_id: dict(
                    path="annotation_batches/" + self.batch_id,
                    batch_result_id=self.packet["batch_result_id"],
                )
            },
        )
        self.manifest["results"] = {
            accession: {
                "eggnog": dict(
                    path=f"samples/{accession}/annotation/eggnog",
                    result_id=record["result_id"],
                )
            }
            for accession, (_, record) in self.results.items()
        }
        self.save_manifest()
        self.set_statuses({accession: "success" for accession in self.accessions})

    def save_manifest(self):
        write_json(self.source / "annotation_results.json", self.manifest)

    def save_result(self, accession, record):
        path = self.results[accession][0]
        record["result_id"] = identity(
            {key: value for key, value in record.items() if key != "result_id"}
        )
        write_json(path / "result.json", record)
        self.results[accession] = path, record
        self.manifest["results"][accession]["eggnog"]["result_id"] = record["result_id"]
        self.save_manifest()

    def set_statuses(self, statuses):
        for name in ("master_table", "sample_status"):
            path = self.source / "tables" / (name + ".tsv")
            rows = read_tsv(path)
            for row in rows:
                row["eggnog_status"] = statuses[row["accession"]]
            write_tsv(path, rows[0], rows)
        path = self.source / "tables/annotation_status.tsv"
        rows = read_tsv(path)
        for row in rows:
            if row["tool"] == "eggnog":
                row["status"] = statuses[row["accession"]]
        write_tsv(path, rows[0], rows)
        self.manifest["complete"] = all(
            status == "success" for status in statuses.values()
        )
        self.save_manifest()
        refresh_tables(self.source)
        self.manifest = read_json(self.source / "annotation_results.json")

    def test_archive_checked_once_for_all_members_and_portable_copy(self):
        copied = self.root / "copied"
        shutil.copytree(self.source, copied)
        self.source.rename(self.root / "unavailable")
        with (
            patch(
                "annotation_result.validate_native_batch", wraps=validate_native_batch
            ) as validate,
            patch(
                "annotation_source.validate_raw_evidence", wraps=validate_raw_evidence
            ) as raw,
        ):
            manifest, batches = validate_source(copied)
        self.assertTrue(manifest["complete"])
        self.assertEqual(validate.call_count, 1)
        self.assertEqual(raw.call_count, 2)
        self.assertTrue(all(call.args[2] is batches for call in raw.call_args_list))
        self.assertEqual(
            batches[self.batch_id].root, copied / "annotation_batches" / self.batch_id
        )
        for accession in self.accessions:
            self.assertFalse(
                (copied / "samples" / accession / "annotation/eggnog/raw").exists()
            )

    def test_import_preserves_upstream_only_and_does_not_republish_batches(self):
        upstream = self.source / "samples/A/upstream.txt"
        upstream.write_text("retained upstream evidence\n")
        output = self.root / "imported"
        import_source(self.source, output)
        self.assertEqual(
            (output / "inherited/samples/A/upstream.txt").read_bytes(),
            upstream.read_bytes(),
        )
        self.assertTrue(
            (output / "inherited/samples/A/annotation/bundle/bundle.json").is_file()
        )
        self.assertFalse((output / "inherited/samples/A/annotation/eggnog").exists())
        self.assertFalse((output / "annotation_batches").exists())
        self.assertFalse((output / "inherited/annotation_batches").exists())

    def test_shared_raw_input_and_packet_identity_tampering_fail(self):
        raw = self.native / "raw/native.txt"
        original = raw.read_bytes()
        raw.write_text("tampered\n")
        with self.assertRaisesRegex(AnnotationError, "evidence changed"):
            validate_source(self.source)
        raw.write_bytes(original)
        original_id = self.manifest["native_batches"][self.batch_id]["batch_result_id"]
        self.manifest["native_batches"][self.batch_id]["batch_result_id"] = "0" * 64
        self.save_manifest()
        with self.assertRaisesRegex(AnnotationError, "result identity changed"):
            validate_source(self.source)
        self.manifest["native_batches"][self.batch_id]["batch_result_id"] = original_id
        self.save_manifest()
        (self.native / "inputs/input.faa").write_text(">altered\nMWA\n")
        with self.assertRaisesRegex(AnnotationError, "input files changed"):
            validate_source(self.source)

    def test_missing_foreign_and_undeclared_batch_archives_fail(self):
        moved = self.root / "missing_archive"
        self.native.rename(moved)
        with self.assertRaisesRegex(AnnotationError, "archives differ"):
            validate_source(self.source)
        moved.rename(self.native)
        entry = self.manifest.pop("native_batches")
        self.save_manifest()
        with self.assertRaisesRegex(AnnotationError, "archives differ"):
            validate_source(self.source)
        self.manifest["native_batches"] = entry
        self.save_manifest()
        self.manifest["native_batches"][self.batch_id]["path"] = "../outside"
        self.save_manifest()
        with self.assertRaisesRegex(AnnotationError, "reference or path"):
            validate_source(self.source)

    def test_packet_id_must_match_its_declared_archive_name(self):
        false_id = "0" * 64
        self.native.rename(self.native.with_name(false_id))
        entry = self.manifest["native_batches"].pop(self.batch_id)
        entry["path"] = "annotation_batches/" + false_id
        self.manifest["native_batches"][false_id] = entry
        self.save_manifest()
        with self.assertRaisesRegex(AnnotationError, "packet IDs differ"):
            validate_source(self.source)

    def test_missing_batch_mapping_cannot_be_treated_as_individual_raw(self):
        self.native.parent.rename(self.root / "unavailable_batches")
        del self.manifest["native_batches"]
        self.save_manifest()
        with self.assertRaisesRegex(
            AnnotationError, "Shared native eggNOG batch is missing"
        ):
            validate_source(self.source)

    def test_linked_packet_manifest_is_rejected_before_indexing(self):
        packet = self.native / "batch_result.json"
        moved = self.root / "linked_packet.json"
        packet.rename(moved)
        packet.symlink_to(moved)
        with patch(
            "annotation_result.validate_native_batch", wraps=validate_native_batch
        ) as validate:
            with self.assertRaisesRegex(AnnotationError, "symlink"):
                validate_source(self.source)
        validate.assert_not_called()

    def test_unused_and_missing_member_references_fail(self):
        original = copy.deepcopy(self.manifest["results"])
        self.manifest["results"] = {}
        self.set_statuses({accession: "failed" for accession in self.accessions})
        with self.assertRaisesRegex(AnnotationError, "unused"):
            validate_source(self.source)
        self.manifest["results"] = {"A": original["A"]}
        self.set_statuses({"A": "success", "B": "failed"})
        with self.assertRaisesRegex(AnnotationError, "omits a member"):
            validate_source(self.source)

    def test_member_reference_must_match_the_batch_packet(self):
        original = copy.deepcopy(self.results["A"][1])
        for change in (
            lambda record: record.update(input_id="0" * 64),
            lambda record: record["native_batch"].update(batch_result_id="0" * 64),
            lambda record: record["search"]["batch"].update(input_id="0" * 64),
            lambda record: record.update(raw_files={}),
        ):
            with self.subTest(change=change):
                record = copy.deepcopy(original)
                change(record)
                self.save_result("A", record)
                with self.assertRaisesRegex(AnnotationError, "Member result differs"):
                    validate_source(self.source)

    def test_foreign_accession_in_an_actual_packed_archive_fails(self):
        inputs, batch = prepare_batches(
            self.fixture.fixture.bundles,
            self.fixture.fixture.method,
            self.root / "foreign_inputs",
            target_fasta_bytes=1000,
        )[0]
        native = self.source / "annotation_batches" / batch["batch_id"]
        native.mkdir()
        shutil.copytree(inputs, native / "inputs")
        shutil.copytree(self.native / "raw", native / "raw")
        packet = dict(
            schema_version=1,
            kind="eggnog_native_batch",
            batch=batch,
            exit_code=0,
            raw_files=inventory(native / "raw"),
        )
        packet["batch_result_id"] = identity(packet)
        write_json(native / "batch_result.json", packet)
        self.manifest["native_batches"][batch["batch_id"]] = dict(
            path="annotation_batches/" + batch["batch_id"],
            batch_result_id=packet["batch_result_id"],
        )
        self.save_manifest()
        with self.assertRaisesRegex(AnnotationError, "foreign accession"):
            validate_source(self.source)

    def test_published_bundle_genome_must_match_the_native_batch_member(self):
        builder = self.fixture.fixture.fixture
        for path in (builder.genome, builder.gff):
            path.write_text(path.read_text().replace("ATGTGAGCTTAA", "ATGTGAGCCTAA"))
        genbank = SeqIO.read(builder.gbk, "genbank")
        genbank.seq = genbank.seq.replace("ATGTGAGCTTAA", "ATGTGAGCCTAA")
        SeqIO.write(genbank, builder.gbk, "genbank")
        updated = builder.build("A", "updated_A")
        bundle = self.source / "samples/A/annotation/bundle"
        bundle.rename(self.root / "previous_A")
        shutil.copytree(updated, bundle)
        self.manifest["bundles"]["A"]["manifest_sha256"] = digest(
            bundle / "bundle.json"
        )
        self.save_manifest()
        with self.assertRaisesRegex(
            AnnotationError, "bundle differs from its native batch member"
        ):
            validate_source(self.source)

    def test_failed_native_archive_is_importable_but_cannot_claim_success(self):
        (self.native / "raw/exit_code.txt").write_text("17\n")
        packet = dict(
            self.packet, exit_code=17, raw_files=inventory(self.native / "raw")
        )
        packet["batch_result_id"] = identity(
            {key: value for key, value in packet.items() if key != "batch_result_id"}
        )
        write_json(self.native / "batch_result.json", packet)
        self.manifest["native_batches"][self.batch_id]["batch_result_id"] = packet[
            "batch_result_id"
        ]
        for accession, (_, original) in list(self.results.items()):
            record = dict(
                original,
                status="failed",
                exit_code=17,
                raw_files=packet["raw_files"],
                native_batch=dict(
                    batch_id=self.batch_id, batch_result_id=packet["batch_result_id"]
                ),
            )
            self.save_result(accession, record)
        self.set_statuses({accession: "failed" for accession in self.accessions})
        manifest, batches = validate_source(self.source)
        self.assertFalse(manifest["complete"])
        self.assertEqual(batches[self.batch_id].record["exit_code"], 17)
        import_source(self.source, self.root / "failed_import")
        self.assertTrue(
            (self.root / "failed_import/inherited/samples/A/annotation/bundle").is_dir()
        )
        record = dict(self.results["A"][1], status="success")
        self.save_result("A", record)
        self.set_statuses({"A": "success", "B": "failed"})
        with self.assertRaisesRegex(AnnotationError, "nonzero or missing native exit"):
            validate_source(self.source)


if __name__ == "__main__":
    unittest.main()
