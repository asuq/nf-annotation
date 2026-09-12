"""Native batch task transitions, immutable copies and all-member completion."""

from __future__ import annotations

import copy
import shutil
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import annotation_batch_tasks as engine
import normalize_eggnog as egg
import test_annotation_normalizers as normalizer_fixture
import test_eggnog_batches as batch_fixture
from annotation_commands import commands, shell_script
from annotation_common import (
    AnnotationError,
    bundle_proteins,
    digest,
    identity,
    read_json,
    write_json,
)
from annotation_resources import RESOURCE_FILE, file_records
from annotation_result import inventory, validate_native_batch, validate_result
from annotation_tasks import POLICY, code_identity, task_identity
from eggnog_batches import prepare_batches


class AnnotationBatchTaskTests(unittest.TestCase):
    def setUp(self):
        self.fixture = batch_fixture.EggnogBatchTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root
        self.normalizer = normalizer_fixture.NormalizerTests()
        self.normalizer.setUp()
        self.addCleanup(self.normalizer.tearDown)
        self.raw, self.resource = self.root / "raw", self.root / "resource"
        self.raw.mkdir()
        self.resource.mkdir()
        self.normalizer.raw, self.normalizer.resource = self.raw, self.resource
        self.manifests = [bundle_proteins(path)[0] for path in self.fixture.bundles]
        proteins = [bundle_proteins(path)[1][0] for path in self.fixture.bundles]
        a, b = proteins[0]["tool_id"], proteins[1]["tool_id"]
        rows = [
            self.annotation(a, KEGG_ko="K00001"),
            self.annotation(b, GOs="GO:0000003", EC="ec:bad", KEGG_ko="K00002"),
        ]
        namespaces = [
            self.namespaces(a),
            self.namespaces(b, gos_mf="GO:0000003"),
        ]
        self.normalizer.write_eggnog(rows, namespaces)
        seeds = self.raw / "eggnog.emapper.seed_orthologs"
        seeds.write_text(
            seeds.read_text().replace(
                "\t1\t100\t1\t100\t90\t50\t50\n",
                "\t1\t3\t1\t100\t90\t100\t50\n",
            )
        )
        (self.raw / "exit_code.txt").write_text("0\n")
        contract = dict(
            component="eggnog",
            version="synthetic-v7",
            settings=dict(fixture=True),
            files=file_records(self.resource),
        )
        self.resource_record = dict(
            schema_version=1,
            component="eggnog",
            resource_id=identity(contract),
            contract=contract,
        )
        write_json(self.resource / RESOURCE_FILE, self.resource_record)
        method = copy.deepcopy(self.fixture.method)
        method["resource_id"] = self.resource_record["resource_id"]
        interpretation = dict(
            schema_version=1, code=code_identity("eggnog"), policy=POLICY["eggnog"]
        )
        self.entry = dict(
            search_method=method,
            interpretation=interpretation,
            method_id=identity(dict(search=method, interpretation=interpretation)),
            container="sha256:" + "1" * 64,
            resource=self.resource_record,
            resource_manifest_sha256=digest(self.resource / RESOURCE_FILE),
        )
        self.batchdir, self.batch = prepare_batches(
            self.fixture.bundles,
            method,
            self.root / "inputs",
            target_fasta_bytes=100_000,
        )[0]
        self.members = [
            engine.member_task_identity(manifest, self.entry, self.batch)
            for manifest in self.manifests
        ]

    def annotation(self, query, **fields):
        row = dict.fromkeys(egg.HEADER, "-")
        row.update(
            query=query,
            seed_ortholog="1234",
            evalue="0",
            score="100",
            COG_category="J",
            **fields,
        )
        row["annotation_confidence"] = "".join(
            "h" if name in fields else "-" for name in egg.FIELDS
        )
        return row

    def namespaces(self, query, **fields):
        row = dict.fromkeys(egg.GO_HEADER, "-")
        row.update(query=query, **fields)
        for name in fields:
            row[name + "_confidence"] = "high"
        return row

    def plan(self, name, previous=None, previous_results=None, entry=None):
        entry = self.entry if entry is None else entry
        members = [
            engine.member_task_identity(manifest, entry, self.batch)
            for manifest in self.manifests
        ]
        return engine.plan_batch(
            self.batchdir,
            entry,
            members,
            self.root / name,
            previous,
            {} if previous_results is None else previous_results,
        )

    def complete(self, name, task="run-task", raw=None):
        return engine.complete_batch(
            self.root / task,
            self.raw if raw is None else raw,
            self.batchdir,
            self.resource,
            self.root / name,
        )

    def published(self, name):
        root = self.root / name
        native = validate_native_batch(
            root / "annotation_batches" / self.batch["batch_id"]
        )
        results = {
            member["accession"]: root
            / "samples"
            / member["accession"]
            / "annotation/eggnog"
            for member in self.members
        }
        return native, results

    def first_run(self):
        self.plan("run-task")
        self.complete("first")
        return self.published("first")

    def test_run_parses_once_and_preserves_empty_go_no_hits_and_raw_footer(self):
        planned = self.plan("run-task")
        self.assertEqual(planned["action"], "run")
        self.assertEqual(
            (self.root / "run-task/run.sh").read_text(),
            shell_script(self.entry["search_method"]["command"]),
        )
        self.assertEqual(
            (self.root / "run-task/input.faa").read_bytes(),
            (self.batchdir / "input.faa").read_bytes(),
        )
        before = inventory(self.raw)
        with patch.object(egg, "normalize_batch", wraps=egg.normalize_batch) as parser:
            result = self.complete("first")
        self.assertEqual(parser.call_count, 1)
        self.assertEqual(list(result["results"]), ["A", "B", "C"])
        native, paths = self.published("first")
        self.assertEqual(native.record["raw_files"], before)
        self.assertNotIn("action", native.record)
        for accession, path in paths.items():
            record = validate_result(path, batches={self.batch["batch_id"]: native})
            self.assertEqual(record, result["results"][accession])
            self.assertFalse((path / "raw").exists())
            self.assertEqual(
                record["input_id"],
                self.manifests[["A", "B", "C"].index(accession)]["input_id"],
            )
            self.assertEqual(
                record["coordinate_id"],
                self.manifests[["A", "B", "C"].index(accession)]["coordinate_id"],
            )
        self.assertEqual(result["results"]["A"]["evidence"]["definitions"], {"go": {}})
        self.assertEqual(result["results"]["C"]["reported_hits"], 0)
        self.assertEqual(result["results"]["C"]["evidence"]["features"], {})
        self.assertEqual(
            {
                error["field"]
                for error in result["results"]["B"]["evidence"]["field_errors"]
            },
            {"go", "ec"},
        )
        for name in egg.TABLE_COLUMNS:
            self.assertEqual(
                (paths["C"] / "normalized" / name).read_text(),
                "\t".join(egg.TABLE_COLUMNS[name]) + "\n",
            )
        for name in ("eggnog.emapper.annotations", "eggnog.emapper.seed_orthologs"):
            self.assertEqual(
                (native.root / "raw" / name).read_bytes(),
                (self.raw / name).read_bytes(),
            )

    def test_reuse_is_uniform_portable_and_retains_native_packet_identity(self):
        previous, paths = self.first_run()
        planned = self.plan("reuse-task", previous, paths)
        self.assertEqual(planned["action"], "reuse")
        self.assertEqual({member["action"] for member in planned["members"]}, {"reuse"})
        self.assertEqual(set(planned["previous_result_ids"]), {"A", "B", "C"})
        self.assertNotIn(
            str(self.root), (self.root / "reuse-task/task.json").read_text()
        )
        portable = self.root / "portable"
        portable.mkdir()
        for source, destination in (
            (self.root / "reuse-task", portable / "task"),
            (self.batchdir, portable / "inputs"),
            (self.resource, portable / "resource"),
        ):
            shutil.copytree(source, destination)
        (self.root / "first").rename(self.root / "unavailable-first")
        self.batchdir.rename(self.root / "unavailable-inputs")
        with patch.object(egg, "normalize_batch") as parser:
            result = engine.complete_batch(
                portable / "task",
                portable / "task/previous_batch/raw",
                portable / "inputs",
                portable / "resource",
                portable / "output",
            )
        parser.assert_not_called()
        packet = portable / "output/annotation_batches" / self.batch["batch_id"]
        self.assertEqual(read_json(packet / "batch_result.json"), previous.record)
        self.assertEqual(result["batch_result_id"], previous.record["batch_result_id"])
        for accession, record in result["results"].items():
            self.assertEqual(record["action"], "reuse")
            self.assertNotEqual(
                record["result_id"], planned["previous_result_ids"][accession]
            )
            source = (
                self.root
                / "unavailable-first/samples"
                / accession
                / "annotation/eggnog/normalized"
            )
            destination = (
                portable / "output/samples" / accession / "annotation/eggnog/normalized"
            )
            self.assertEqual(inventory(source), inventory(destination))

    def test_interpretation_or_missing_member_renormalizes_the_whole_batch_once(self):
        previous, paths = self.first_run()
        changed = copy.deepcopy(self.entry)
        changed["interpretation"]["fixture_revision"] = 2
        changed["method_id"] = identity(
            dict(
                search=changed["search_method"],
                interpretation=changed["interpretation"],
            )
        )
        for name, entry, previous_results in (
            ("changed", changed, paths),
            ("missing", self.entry, {"A": paths["A"]}),
        ):
            with self.subTest(name=name):
                planned = self.plan(name + "-task", previous, previous_results, entry)
                self.assertEqual(planned["action"], "renormalize")
                self.assertNotIn("previous_result_ids", planned)
                self.assertFalse(
                    (self.root / (name + "-task/previous_results")).exists()
                )
                with patch.object(
                    egg, "normalize_batch", wraps=egg.normalize_batch
                ) as parser:
                    result = self.complete(
                        name,
                        task=name + "-task",
                        raw=self.root / (name + "-task/previous_batch/raw"),
                    )
                self.assertEqual(parser.call_count, 1)
                self.assertEqual(
                    result["batch_result_id"], previous.record["batch_result_id"]
                )
                self.assertEqual(
                    {record["action"] for record in result["results"].values()},
                    {"renormalize"},
                )

    def test_membership_method_and_individual_results_cannot_reuse_native_batch(self):
        previous, paths = self.first_run()
        for name, bundles, manifests, entry in (
            ("membership", self.fixture.bundles[:2], self.manifests[:2], self.entry),
            ("method", self.fixture.bundles, self.manifests, copy.deepcopy(self.entry)),
        ):
            if name == "method":
                entry["search_method"]["command"] = commands("eggnog", 8, 64)
                entry["method_id"] = identity(
                    dict(
                        search=entry["search_method"],
                        interpretation=entry["interpretation"],
                    )
                )
            path, batch = prepare_batches(
                bundles,
                entry["search_method"],
                self.root / (name + "-inputs"),
                target_fasta_bytes=100_000,
            )[0]
            members = [
                engine.member_task_identity(manifest, entry, batch)
                for manifest in manifests
            ]
            result = engine.plan_batch(
                path, entry, members, self.root / (name + "-task"), previous, paths
            )
            self.assertEqual(result["action"], "run")
            self.assertNotEqual(
                members[0]["search_fingerprint"], self.members[0]["search_fingerprint"]
            )
        individual = task_identity(self.manifests[0], "eggnog", self.entry)
        self.assertNotEqual(
            individual["search_fingerprint"], self.members[0]["search_fingerprint"]
        )
        self.assertEqual(self.plan("individual-task", None, paths)["action"], "run")

    def test_native_failure_and_native_parse_failure_publish_all_failed_members(self):
        self.plan("run-task")
        original = (self.raw / "eggnog.emapper.seed_orthologs").read_text()
        for name, exit_code in (("native-failed", 7), ("parse-failed", 0)):
            (self.raw / "exit_code.txt").write_text(f"{exit_code}\n")
            if exit_code == 0:
                (self.raw / "eggnog.emapper.seed_orthologs").write_text(
                    original.replace("## 2 queries scanned\n", "")
                )
            with patch.object(
                egg, "normalize_batch", wraps=egg.normalize_batch
            ) as parser:
                result = self.complete(name)
            self.assertEqual(parser.call_count, int(exit_code == 0))
            native, paths = self.published(name)
            for record in result["results"].values():
                self.assertEqual(record["status"], "failed")
                self.assertIsNone(record["reported_hits"])
                self.assertIsNone(record["accepted_proteins"])
                self.assertIsNone(record["evidence"])
                self.assertEqual(record["normalized_files"], {})
                self.assertEqual(record["raw_files"], native.record["raw_files"])
            self.assertTrue(
                all(not (path / "normalized").exists() for path in paths.values())
            )
            if exit_code != 0:
                self.assertEqual(
                    self.plan("retry-task", native, paths)["action"], "run"
                )

    def test_task_member_code_resource_and_copied_input_changes_fail_explicitly(self):
        task = self.plan("run-task")
        task_path = self.root / "run-task/task.json"
        tampered = copy.deepcopy(task)
        tampered["members"][0]["input_id"] = "0" * 64
        tampered["task_id"] = identity(
            {key: value for key, value in tampered.items() if key != "task_id"}
        )
        write_json(task_path, tampered)
        with self.assertRaisesRegex(AnnotationError, "Task member differs"):
            self.complete("member-changed")
        write_json(task_path, task)
        with patch.object(engine, "code_identity", return_value={}):
            with self.assertRaisesRegex(AnnotationError, "method/code"):
                self.complete("code-changed")
        resource_path = self.resource / RESOURCE_FILE
        resource_bytes = resource_path.read_bytes()
        resource_path.write_bytes(resource_bytes + b" ")
        with self.assertRaisesRegex(AnnotationError, "Task/resource identity"):
            self.complete("resource-changed")
        resource_path.write_bytes(resource_bytes)
        (self.root / "run-task/input.faa").write_text(">wrong\nMWA\n")
        with self.assertRaisesRegex(AnnotationError, "input or native command changed"):
            self.complete("input-changed")
        for name in (
            "member-changed",
            "code-changed",
            "resource-changed",
            "input-changed",
        ):
            self.assertFalse((self.root / name).exists())

    def test_staged_directory_roots_are_allowed_but_internal_raw_links_fail(self):
        self.plan("run-task")
        staged = self.root / "staged"
        staged.mkdir()
        for name, source in (
            ("task", self.root / "run-task"),
            ("raw", self.raw),
            ("batch", self.batchdir),
            ("resource", self.resource),
        ):
            (staged / name).symlink_to(source, target_is_directory=True)
        result = engine.complete_batch(
            staged / "task",
            staged / "raw",
            staged / "batch",
            staged / "resource",
            self.root / "staged-output",
        )
        self.assertEqual(
            {record["status"] for record in result["results"].values()}, {"success"}
        )
        (self.raw / "linked.log").symlink_to(self.raw / "exit_code.txt")
        with self.assertRaisesRegex(AnnotationError, "Linked evidence"):
            engine.complete_batch(
                staged / "task",
                staged / "raw",
                staged / "batch",
                staged / "resource",
                self.root / "linked-output",
            )
        self.assertFalse((self.root / "linked-output").exists())

    def test_raw_and_normalized_tampering_fail_at_planning_and_completion_boundaries(
        self,
    ):
        previous, paths = self.first_run()
        normalized = paths["A"] / "normalized/eggnog_annotations.tsv"
        original = normalized.read_bytes()
        normalized.write_bytes(original + b"\n")
        with self.assertRaisesRegex(
            AnnotationError, "Normalized annotation evidence changed"
        ):
            self.plan("bad-normalized-task", previous, paths)
        normalized.write_bytes(original)
        self.plan("reuse-task", previous, paths)
        archive = self.root / "reuse-task/previous_results/A"
        (archive / "normalized/eggnog_annotations.tsv").write_bytes(original + b"\n")
        old = read_json(archive / "result.json")
        old["normalized_files"] = inventory(archive / "normalized")
        old["result_id"] = identity(
            {key: value for key, value in old.items() if key != "result_id"}
        )
        write_json(archive / "result.json", old)
        with self.assertRaisesRegex(AnnotationError, "differs from the planned result"):
            self.complete(
                "substituted-member",
                task="reuse-task",
                raw=self.root / "reuse-task/previous_batch/raw",
            )
        self.plan("second-reuse-task", previous, paths)
        packet = self.root / "second-reuse-task/previous_batch"
        (packet / "raw/tool.log").write_text("substituted raw\n")
        record = read_json(packet / "batch_result.json")
        record["raw_files"] = inventory(packet / "raw")
        record["batch_result_id"] = identity(
            {key: value for key, value in record.items() if key != "batch_result_id"}
        )
        write_json(packet / "batch_result.json", record)
        with self.assertRaisesRegex(AnnotationError, "differs from the planned packet"):
            self.complete(
                "substituted-native", task="second-reuse-task", raw=packet / "raw"
            )
        (previous.root / "raw/exit_code.txt").write_text("1\n")
        with self.assertRaisesRegex(AnnotationError, "Native batch evidence changed"):
            self.plan("bad-raw-task", previous, paths)
        self.assertFalse((self.root / "substituted-member").exists())
        self.assertFalse((self.root / "substituted-native").exists())


if __name__ == "__main__":
    unittest.main()
