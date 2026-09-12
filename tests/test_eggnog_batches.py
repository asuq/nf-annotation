"""Identity, conservation and memory-bounded packing for eggNOG query batches."""

from __future__ import annotations

import copy
import shutil
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import test_annotation_bundle as bundle_fixture
from annotation_commands import commands
from annotation_common import (
    AnnotationError,
    bundle_proteins,
    digest,
    identity,
    read_json,
    write_json,
)
from eggnog_batches import prepare_batches, validate_batch
from eggnog_native import NATIVE_CODE_FILES


class EggnogBatchTests(unittest.TestCase):
    def setUp(self):
        self.fixture = bundle_fixture.AnnotationBundleTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root
        self.bundles = [
            self.fixture.build(name, f"bundle_{name}") for name in ("A", "B", "C")
        ]
        self.method = {
            "tool": "eggnog",
            "runtime_id": "sha256:" + "1" * 64,
            "resource_id": "2" * 64,
            "command": commands("eggnog", 16, 64),
            "native_code": {
                name: digest(Path(__file__).resolve().parents[1] / "bin" / name)
                for name in NATIVE_CODE_FILES
            },
        }
        self.size = (self.bundles[0] / "proteins.faa").stat().st_size

    def batches(self, name="batches", bundles=None, target=None, method=None):
        return prepare_batches(
            self.bundles if bundles is None else bundles,
            self.method if method is None else method,
            self.root / name,
            target_fasta_bytes=2 * self.size if target is None else target,
        )

    def reseal(self, path, change):
        """Create internally hashed malformed receipts to exercise semantic checks."""
        record = read_json(path / "batch.json")
        change(record)
        inputs = {
            key: value
            for key, value in record.items()
            if key
            not in ("input_id", "search_method", "search_fingerprint", "batch_id")
        }
        record["input_id"] = identity(inputs)
        record["search_fingerprint"] = identity(
            {"input_id": record["input_id"], "method": record["search_method"]}
        )
        record["batch_id"] = identity(
            {key: value for key, value in record.items() if key != "batch_id"}
        )
        write_json(path / "batch.json", record)

    def test_order_independent_whole_proteome_packing_and_portability(self):
        original = self.batches()
        reordered = self.batches("reordered", list(reversed(self.bundles)))
        self.assertEqual(
            [record["batch_id"] for _, record in original],
            [record["batch_id"] for _, record in reordered],
        )
        self.assertEqual(
            [
                [member["accession"] for member in record["members"]]
                for _, record in original
            ],
            [["A", "B"], ["C"]],
        )
        self.assertEqual(
            (original[0][0] / "input.faa").read_bytes(),
            b"".join(
                (bundle / "proteins.faa").read_bytes() for bundle in self.bundles[:2]
            ),
        )
        copied = self.root / "relocated"
        shutil.copytree(original[0][0], copied)
        record, proteins = validate_batch(
            copied, expected_search_fingerprint=original[0][1]["search_fingerprint"]
        )
        self.assertEqual(record, original[0][1])
        # Original protein IDs are sample-local and must not be rewritten.
        self.assertEqual([row["protein_id"] for row in proteins], ["gene_1", "gene_1"])
        self.assertEqual(len({row["tool_id"] for row in proteins}), 2)
        self.assertNotIn(str(self.root), (copied / "batch.json").read_text())

    def test_large_proteome_is_explicit_unmodified_singleton(self):
        batches = self.batches(target=self.size - 1)
        self.assertEqual(len(batches), 3)
        for (path, record), bundle in zip(batches, self.bundles, strict=True):
            self.assertTrue(record["oversized"])
            self.assertEqual(len(record["members"]), 1)
            self.assertEqual(
                (path / "input.faa").read_bytes(),
                (bundle / "proteins.faa").read_bytes(),
            )

    def test_membership_packing_policy_and_method_invalidate_search(self):
        original = self.batches("original", self.bundles[:2])[0][1]
        standalone = self.batches("standalone", self.bundles[:1])[0][1]
        target = self.batches("target", self.bundles[:2], target=3 * self.size)[0][1]
        method = copy.deepcopy(self.method)
        method["resource_id"] = "3" * 64
        new_method = self.batches("method", self.bundles[:2], method=method)[0][1]
        self.assertEqual(original["input_id"], new_method["input_id"])
        self.assertEqual(
            len(
                {
                    record["search_fingerprint"]
                    for record in (original, standalone, target, new_method)
                }
            ),
            4,
        )
        with self.assertRaisesRegex(AnnotationError, "fingerprint"):
            validate_batch(
                self.root / "original/batch00000000",
                expected_search_fingerprint=standalone["search_fingerprint"],
            )

    def test_changed_sequence_fails_even_after_bundle_checksum_refresh(self):
        bundle = self.bundles[0]
        path = bundle / "proteins.faa"
        path.write_text(path.read_text().replace("MWA", "MAA"))
        record = read_json(bundle / "bundle.json")
        record["files"]["proteins.faa"] = digest(path)
        write_json(bundle / "bundle.json", record)
        with self.assertRaisesRegex(AnnotationError, "checksum"):
            self.batches()

    def test_multiline_native_fasta_bytes_are_preserved(self):
        bundle = self.bundles[0]
        path = bundle / "proteins.faa"
        path.write_text(path.read_text().replace("MWA\n", "M\nWA\n"))
        record = read_json(bundle / "bundle.json")
        record["files"]["proteins.faa"] = digest(path)
        write_json(bundle / "bundle.json", record)
        batches = self.batches(bundles=[bundle], target=self.size + 10)
        self.assertEqual((batches[0][0] / "input.faa").read_bytes(), path.read_bytes())

    def test_duplicate_cohort_and_invalid_target_fail(self):
        with self.assertRaisesRegex(AnnotationError, "duplicate accessions"):
            self.batches(bundles=[self.bundles[0], self.bundles[0]])
        for target in (0, -1, True, 1.5):
            with self.assertRaisesRegex(AnnotationError, "positive integer"):
                self.batches(target=target)
        with self.assertRaisesRegex(AnnotationError, "empty"):
            self.batches(bundles=[])

    def test_missing_or_mutable_search_provenance_fails(self):
        for key, value in (
            ("runtime_id", "eggnog:latest"),
            ("resource_id", "resource"),
            ("command", {}),
        ):
            method = copy.deepcopy(self.method)
            method[key] = value
            with self.assertRaises(AnnotationError):
                self.batches(method=method)

    def test_global_id_collision_rejected_across_singleton_batches(self):
        a, b = (bundle_proteins(bundle) for bundle in self.bundles[:2])
        b[1][0]["tool_id"] = a[1][0]["tool_id"]
        # Inject a synthetic hash collision after canonical validation. This
        # exercises cohort-wide detection across already emitted batch files.
        with (
            patch("eggnog_batches.bundle_proteins", side_effect=[a, b]),
            patch("eggnog_batches.validate_proteins"),
            self.assertRaisesRegex(AnnotationError, "cross-sample"),
        ):
            self.batches(bundles=self.bundles[:2], target=self.size - 1)

    def test_changed_mapping_and_linked_input_fail(self):
        batches = self.batches()
        path = batches[0][0]
        mapping = path / "mapping.tsv"
        mapping.write_text(mapping.read_text().replace("A::gene_1", "B::gene_1"))
        self.reseal(
            path,
            lambda record: record["files"].update({"mapping.tsv": digest(mapping)}),
        )
        with self.assertRaisesRegex(AnnotationError, "mapping"):
            validate_batch(path)
        path = batches[1][0]
        fasta = path / "input.faa"
        fasta.unlink()
        fasta.symlink_to(self.bundles[2] / "proteins.faa")
        with self.assertRaisesRegex(AnnotationError, "linked"):
            validate_batch(path)

    def test_archival_code_integrity_and_current_plan_invalidation(self):
        path, original = self.batches()[0]
        updated_code = {
            "eggnog_batches.py": "8" * 64,
            "annotation_common.py": "9" * 64,
        }
        with patch("eggnog_batches.packing_code_identity", return_value=updated_code):
            archived, _ = validate_batch(path)
            self.assertEqual(archived, original)
            new_record = self.batches("new_code")[0][1]
        self.assertNotEqual(original["input_id"], new_record["input_id"])
        self.assertNotEqual(
            original["search_fingerprint"], new_record["search_fingerprint"]
        )
        with self.assertRaisesRegex(AnnotationError, "planned fingerprint"):
            validate_batch(
                path, expected_search_fingerprint=new_record["search_fingerprint"]
            )
        self.reseal(
            path,
            lambda record: record["packing_code"].update({"eggnog_batches.py": "bad"}),
        )
        with self.assertRaisesRegex(AnnotationError, "packing code identity"):
            validate_batch(path)

    def test_malformed_member_boundaries_and_oversize_policy_fail(self):
        for name, change, message in (
            (
                "offset",
                lambda record: record["members"][1].update(byte_offset=0),
                "boundaries",
            ),
            (
                "member",
                lambda record: record["members"][0].update(protein_count=0),
                "boundaries",
            ),
            (
                "boolean_offset",
                lambda record: record["members"][0].update(byte_offset=False),
                "boundaries",
            ),
            ("oversized", lambda record: record.update(oversized=True), "oversized"),
            (
                "schema",
                lambda record: record["members"][0].pop("input_id"),
                "membership",
            ),
        ):
            path = self.batches(name)[0][0]
            self.reseal(path, change)
            with self.assertRaisesRegex(AnnotationError, message):
                validate_batch(path)


if __name__ == "__main__":
    unittest.main()
