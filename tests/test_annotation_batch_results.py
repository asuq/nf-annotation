"""Portable shared raw evidence must remain exact, complete and member-specific."""

from __future__ import annotations

import copy
import shutil
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import test_eggnog_batches as fixture
from annotation_common import AnnotationError, identity, write_json
from annotation_result import (
    batch_search_identity,
    inventory,
    native_batch_index,
    validate_native_batch,
    validate_raw_evidence,
    validate_result,
)
from eggnog_native_fixture import phased_raw


class AnnotationBatchResultTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixture.EggnogBatchTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root
        inputs, self.batch = self.fixture.batches()[0]
        self.native = (
            self.root / "published" / "annotation_batches" / self.batch["batch_id"]
        )
        self.native.mkdir(parents=True)
        shutil.copytree(inputs, self.native / "inputs")
        (self.native / "raw").mkdir()
        (self.native / "raw/exit_code.txt").write_text("0\n")
        (self.native / "raw/native.txt").write_text("unchanged native batch evidence\n")
        phased_raw(self.native / "inputs", self.native / "raw")
        self.native_record = {
            "schema_version": 1,
            "kind": "eggnog_native_batch",
            "batch": self.batch,
            "exit_code": 0,
            "raw_files": inventory(self.native / "raw"),
        }
        self.native_record["batch_result_id"] = identity(self.native_record)
        write_json(self.native / "batch_result.json", self.native_record)
        self.results = []
        for member in self.batch["members"]:
            result = (
                self.root
                / "published/samples"
                / member["accession"]
                / "annotation/eggnog"
            )
            (result / "normalized").mkdir(parents=True)
            (result / "normalized/evidence.tsv").write_text("gene_id\n")
            record = {
                "schema_version": 1,
                "tool": "eggnog",
                "status": "success",
                "accession": member["accession"],
                "input_id": member["input_id"],
                "input_proteins": member["protein_count"],
                "search": {
                    **{
                        key: member[key]
                        for key in ("input_id", "genetic_code", "source_genome_sha256")
                    },
                    "method": self.batch["search_method"],
                    "batch": batch_search_identity(self.batch),
                },
                "exit_code": 0,
                "raw_files": self.native_record["raw_files"],
                "native_batch": {
                    "batch_id": self.batch["batch_id"],
                    "batch_result_id": self.native_record["batch_result_id"],
                },
                "normalized_files": inventory(result / "normalized"),
            }
            self.save_result(result, record)
            self.results.append((result, record))

    def save_result(self, path, record):
        record["result_id"] = identity(
            {key: value for key, value in record.items() if key != "result_id"}
        )
        write_json(path / "result.json", record)

    def test_shared_native_files_are_checked_once_and_copy_is_portable(self):
        with patch(
            "annotation_result.validate_native_batch", wraps=validate_native_batch
        ) as validate:
            batches = native_batch_index([self.native])
            for result, record in self.results:
                self.assertEqual(validate_result(result, batches=batches), record)
            self.assertEqual(validate.call_count, 1)
        copied = self.root / "copied"
        shutil.copytree(self.root / "published", copied)
        (self.root / "published").rename(self.root / "unavailable")
        copied_batches = native_batch_index(
            [copied / "annotation_batches" / self.batch["batch_id"]]
        )
        for _, record in self.results:
            result = copied / "samples" / record["accession"] / "annotation/eggnog"
            self.assertEqual(validate_result(result, batches=copied_batches), record)
            self.assertFalse((result / "raw").exists())

    def test_missing_duplicate_and_changed_shared_evidence_fail(self):
        with self.assertRaisesRegex(AnnotationError, "missing"):
            validate_result(self.results[0][0])
        with self.assertRaisesRegex(AnnotationError, "Duplicate"):
            native_batch_index([self.native, self.native])
        (self.native / "raw/native.txt").write_text("tampered\n")
        with self.assertRaisesRegex(AnnotationError, "evidence changed"):
            native_batch_index([self.native])

    def test_foreign_member_method_and_raw_reference_fail(self):
        batches = native_batch_index([self.native])
        path, original = self.results[0]
        for change in (
            lambda record: record.update(accession="foreign"),
            lambda record: record.update(input_id="0" * 64),
            lambda record: record.update(input_proteins=True),
            lambda record: record.update(exit_code=False),
            lambda record: record.update(raw_files={}),
            lambda record: record.update(native_batch=None),
            lambda record: record["native_batch"].update(batch_result_id="0" * 64),
            lambda record: record["search"].update(genetic_code=11),
            lambda record: record["search"]["method"].update(resource_id="0" * 64),
            lambda record: record["search"]["batch"].update(
                search_fingerprint="0" * 64
            ),
        ):
            with self.subTest(change=change):
                record = copy.deepcopy(original)
                change(record)
                self.save_result(path, record)
                with self.assertRaises(AnnotationError):
                    validate_result(path, batches=batches)

    def test_shared_search_cannot_be_relabelled_as_individual_raw(self):
        path, record = self.results[0]
        del record["native_batch"]
        shutil.copytree(self.native / "raw", path / "raw")
        with self.assertRaisesRegex(AnnotationError, "reference is missing"):
            validate_raw_evidence(path, record)

    def test_linked_raw_and_inputs_fail(self):
        for child in ("raw", "inputs"):
            with self.subTest(child=child):
                original = self.native / child
                moved = self.root / ("linked_" + child)
                original.rename(moved)
                original.symlink_to(moved)
                with self.assertRaisesRegex(AnnotationError, "Linked"):
                    validate_native_batch(self.native)
                original.unlink()
                moved.rename(original)

    def test_failed_native_execution_cannot_be_declared_successful(self):
        (self.native / "raw/exit_code.txt").write_text("4\n")
        self.native_record["exit_code"] = 4
        self.native_record["raw_files"] = inventory(self.native / "raw")
        self.native_record["batch_result_id"] = identity(
            {
                key: value
                for key, value in self.native_record.items()
                if key != "batch_result_id"
            }
        )
        write_json(self.native / "batch_result.json", self.native_record)
        batches = native_batch_index([self.native])
        path, record = self.results[0]
        record["exit_code"] = 4
        record["raw_files"] = self.native_record["raw_files"]
        record["native_batch"]["batch_result_id"] = self.native_record[
            "batch_result_id"
        ]
        self.save_result(path, record)
        with self.assertRaisesRegex(AnnotationError, "nonzero"):
            validate_result(path, batches=batches)
        record["status"] = "failed"
        validate_raw_evidence(path, record, batches)


if __name__ == "__main__":
    unittest.main()
