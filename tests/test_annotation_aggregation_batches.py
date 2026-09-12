"""Shared native evidence must be complete and portable at cohort publication."""

from __future__ import annotations

import copy
import shutil
import subprocess
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import aggregate_annotations
import annotation_result
import test_annotation_batch_results as batch_fixture
from annotation_common import (
    TOOLS,
    AnnotationError,
    bundle_proteins,
    identity,
    read_json,
    read_tsv,
    write_json,
    write_tsv,
)
from annotation_result import batch_search_identity, inventory


class AnnotationAggregationBatchTests(unittest.TestCase):
    def setUp(self):
        self.fixture = batch_fixture.AnnotationBatchResultTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root
        self.bundles = self.fixture.fixture.bundles
        self.native = self.fixture.native
        self.batch = self.fixture.batch
        self.results = [path for path, _ in self.fixture.results]
        self.accessions = ["B", "C", "A"]
        self.master, self.status, self.plan_path = [
            self.root / name for name in ("master.tsv", "status.tsv", "plan.json")
        ]
        interpretation = {
            "schema_version": 1,
            "policy": "synthetic aggregation fixture",
        }
        for path, original in self.fixture.results:
            record = copy.deepcopy(original)
            record["interpretation"] = interpretation
            gene = record["accession"] + "::gene_1"
            hit = record["accession"] == "A"
            record.update(
                reason="validated_native_output",
                action="run",
                reported_hits=int(hit),
                evidence=dict(
                    mapped_genes=[gene] if hit else [],
                    accepted_genes=[gene] if hit else [],
                    features={"ko": {gene: ["K00001"]}} if hit else {},
                    definitions={},
                    field_errors=[],
                ),
            )
            self.save_result(path, record)
        # A third sample retains individual native evidence alongside the batch.
        bundle, proteins = bundle_proteins(self.bundles[2])
        path = self.root / "published/samples/C/annotation/eggnog"
        (path / "raw").mkdir(parents=True)
        (path / "raw/native.txt").write_text("individual native evidence\n")
        (path / "normalized").mkdir()
        (path / "normalized/evidence.tsv").write_text("gene_id\n")
        gene = proteins[0]["gene_id"]
        record = dict(
            schema_version=1,
            tool="eggnog",
            status="success",
            accession="C",
            input_id=bundle["input_id"],
            input_proteins=1,
            search={
                **{
                    key: bundle[key]
                    for key in ("input_id", "genetic_code", "source_genome_sha256")
                },
                "method": self.batch["search_method"],
            },
            exit_code=0,
            raw_files=inventory(path / "raw"),
            normalized_files=inventory(path / "normalized"),
            interpretation=interpretation,
            reason="validated_native_output",
            action="run",
            reported_hits=1,
            evidence=dict(
                mapped_genes=[gene],
                accepted_genes=[gene],
                features={"ko": {gene: ["K00002"]}},
                definitions={},
                field_errors=[],
            ),
        )
        self.save_result(path, record)
        self.results.append(path)
        self.save_plan()

    def save_result(self, path, record):
        record["search_fingerprint"] = identity(record["search"])
        record["normalization_fingerprint"] = identity(
            dict(
                search_fingerprint=record["search_fingerprint"],
                interpretation=record["interpretation"],
            )
        )
        record["method_id"] = identity(
            dict(
                search=record["search"]["method"],
                interpretation=record["interpretation"],
            )
        )
        self.fixture.save_result(path, record)

    def save_plan(self):
        records = {
            read_json(path / "result.json")["accession"]: read_json(
                path / "result.json"
            )
            for path in self.results
        }
        tasks = []
        for accession in self.accessions:
            for tool in TOOLS:
                task = dict(
                    accession=accession,
                    tool=tool,
                    input_proteins=1,
                    action="run" if tool == "eggnog" else "skip",
                    status="planned" if tool == "eggnog" else "skipped_disabled",
                    reason="fixture",
                )
                if tool == "eggnog":
                    task.update(
                        {
                            name: records[accession][name]
                            for name in (
                                "search_fingerprint",
                                "normalization_fingerprint",
                                "method_id",
                            )
                        }
                    )
                tasks.append(task)
        plan = dict(
            schema_version=1,
            accessions=self.accessions,
            enabled_tools=["eggnog"],
            tasks=tasks,
        )
        plan["plan_id"] = identity(plan)
        write_json(self.plan_path, plan)
        write_tsv(
            self.master,
            ("Accession",),
            [dict(Accession=acc) for acc in self.accessions],
        )
        write_tsv(
            self.status,
            ("accession",),
            [dict(accession=acc) for acc in reversed(self.accessions)],
        )

    def aggregate(self, name, *, results=None, batches=None):
        return aggregate_annotations.aggregate(
            self.plan_path,
            self.master,
            self.status,
            self.bundles,
            self.results if results is None else results,
            self.root / name,
            batch_dirs=[self.native] if batches is None else batches,
        )

    def make_individual(self, path):
        record = read_json(path / "result.json")
        del record["native_batch"]
        del record["search"]["batch"]
        shutil.copytree(self.native / "raw", path / "raw")
        self.save_result(path, record)

    def test_shared_batch_is_verified_once_and_manifest_is_portable(self):
        with patch.object(
            annotation_result,
            "validate_native_batch",
            wraps=annotation_result.validate_native_batch,
        ) as validate:
            manifest = self.aggregate("report")
        self.assertEqual(validate.call_count, 1)
        batch_id = self.batch["batch_id"]
        self.assertEqual(
            manifest["native_batches"],
            {
                batch_id: dict(
                    path=f"annotation_batches/{batch_id}",
                    batch_result_id=self.fixture.native_record["batch_result_id"],
                ),
            },
        )
        self.assertNotIn(
            str(self.root), (self.root / "report/annotation_results.json").read_text()
        )
        self.assertEqual(
            (
                self.root / "report/tables/functional_matrices/eggnog_ko_counts.tsv"
            ).read_text(),
            "accession\tK00001\tK00002\nB\t0\t0\nC\t0\t1\nA\t1\t0\n",
        )
        self.assertTrue(manifest["complete"])
        self.assertFalse((self.results[0] / "raw").exists())
        self.assertFalse((self.results[1] / "raw").exists())

    def test_duplicate_missing_unused_and_unreferenced_members_fail_before_publication(
        self,
    ):
        for name, batches, results, error in (
            ("duplicate", [self.native, self.native], self.results, "Duplicate"),
            ("missing", [], self.results, "batch is missing"),
            ("member", [self.native], self.results[1:], "Missing result references"),
        ):
            with (
                self.subTest(name=name),
                self.assertRaisesRegex(AnnotationError, error),
            ):
                self.aggregate(name, batches=batches, results=results)
            self.assertFalse((self.root / name).exists())
        for path in self.results[:2]:
            self.make_individual(path)
        self.save_plan()
        with self.assertRaisesRegex(AnnotationError, "Unused native"):
            self.aggregate("unused")
        self.assertFalse((self.root / "unused").exists())
        self.assertNotIn("native_batches", self.aggregate("individual", batches=[]))
        command = self.cli("individual-cli")
        write_json(self.root / "batch-list with spaces.json", [])
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        for path in (self.root / "individual").rglob("*"):
            if path.is_file():
                self.assertEqual(
                    path.read_bytes(),
                    (
                        self.root
                        / "individual-cli"
                        / path.relative_to(self.root / "individual")
                    ).read_bytes(),
                )

    def test_foreign_declared_member_and_foreign_reference_fail(self):
        original = self.accessions
        self.accessions = ["C", "A"]
        self.save_plan()
        with self.assertRaisesRegex(AnnotationError, "undeclared cohort members"):
            self.aggregate("foreign-member")
        self.accessions = original
        self.save_plan()
        path = self.results[2]
        record = read_json(path / "result.json")
        record["native_batch"] = copy.deepcopy(
            read_json(self.results[0] / "result.json")["native_batch"]
        )
        record["search"]["batch"] = batch_search_identity(self.batch)
        record["raw_files"] = self.fixture.native_record["raw_files"]
        shutil.rmtree(path / "raw")
        self.save_result(path, record)
        self.save_plan()
        with self.assertRaisesRegex(
            AnnotationError, "differs from its shared native batch"
        ):
            self.aggregate("foreign-reference")
        self.assertFalse((self.root / "foreign-member").exists())
        self.assertFalse((self.root / "foreign-reference").exists())

    def test_mismatched_reference_and_tampered_native_or_normalized_evidence_fail(self):
        path = self.results[0]
        original = read_json(path / "result.json")
        changed = copy.deepcopy(original)
        changed["native_batch"]["batch_result_id"] = "0" * 64
        self.save_result(path, changed)
        with self.assertRaisesRegex(
            AnnotationError, "differs from its shared native batch"
        ):
            self.aggregate("reference")
        self.save_result(path, original)
        normalized = path / "normalized/evidence.tsv"
        normalized.write_text("changed\n")
        with self.assertRaisesRegex(AnnotationError, "Normalized evidence changed"):
            self.aggregate("normalized")
        (self.native / "raw/native.txt").write_text("changed raw\n")
        with self.assertRaisesRegex(AnnotationError, "Native batch evidence changed"):
            self.aggregate("native")
        self.assertFalse(list(self.root.glob(".annotation-aggregation-*")))

    def test_failed_member_remains_referenced_and_counts_unavailable(self):
        path = self.results[1]
        record = read_json(path / "result.json")
        record.update(
            status="failed",
            reason="normalization_failed",
            evidence=None,
            normalized_files={},
        )
        self.save_result(path, record)
        manifest = self.aggregate("failed")
        self.assertFalse(manifest["complete"])
        self.assertIn("native_batches", manifest)
        rows = read_tsv(
            self.root / "failed/tables/functional_matrices/eggnog_ko_counts.tsv"
        )
        self.assertEqual(rows[0], dict(accession="B", K00001="NA", K00002="NA"))

    def cli(self, output="cli"):
        command = [
            sys.executable,
            str(Path(aggregate_annotations.__file__)),
            "--plan",
            str(self.plan_path),
            "--master",
            str(self.master),
            "--sample-status",
            str(self.status),
            "--output",
            str(self.root / output),
        ]
        for name, paths in (
            ("bundle-list", self.bundles),
            ("result-list", list(reversed(self.results))),
            ("batch-list", [self.native]),
        ):
            path = self.root / (name + " with spaces.json")
            write_json(path, [str(item) for item in paths])
            command.extend(["--" + name, str(path)])
        return command

    def test_cli_explicit_lists_match_api_and_reject_legacy_or_invalid_lists(self):
        self.aggregate("api")
        command = self.cli()
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        for path in (self.root / "api").rglob("*"):
            if path.is_file():
                self.assertEqual(
                    path.read_bytes(),
                    (
                        self.root / "cli" / path.relative_to(self.root / "api")
                    ).read_bytes(),
                )
        for flag in ("--bundle", "--result", "--batch"):
            with self.subTest(flag=flag):
                result = subprocess.run(
                    command + [flag, "unused"], capture_output=True, text=True
                )
                self.assertEqual(result.returncode, 2)
                self.assertIn("unrecognized arguments", result.stderr)
        command = self.cli("invalid")
        write_json(
            self.root / "batch-list with spaces.json", {"paths": [str(self.native)]}
        )
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 1)
        self.assertIn("Expected a JSON array", result.stderr)
        self.assertFalse((self.root / "invalid").exists())
        command = self.cli("required")
        index = command.index("--batch-list")
        del command[index : index + 2]
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("--batch-list", result.stderr)
        self.assertFalse((self.root / "required").exists())


if __name__ == "__main__":
    unittest.main()
