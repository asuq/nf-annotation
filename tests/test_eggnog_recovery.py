"""Recovery executes only unfinished proteomes, with unchanged search evidence."""

from __future__ import annotations

import shutil
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import test_eggnog_native_execution as native_fixture
from annotation_common import digest, identity, read_json, write_json
from annotation_resources import RESOURCE_FILE, file_records
from annotation_result import inventory
from eggnog_batches import prepare_batches
from recover_eggnog_batch import ANNOTATION_FILES


class EggnogRecoveryTests(unittest.TestCase):
    def setUp(self):
        self.native = native_fixture.EggnogNativeExecutionTests()
        self.native.setUp()
        self.addCleanup(self.native.doCleanups)
        self.root, self.fixture = self.native.root, self.native.fixture
        # The synthetic database needs the same six-file layout as the native mapper.
        for name in (*ANNOTATION_FILES, "eggnog_proteins.dmnd"):
            path = self.fixture.resource / name
            if not path.exists():
                path.write_bytes(b"Synthetic annotation database control\n")
        record = read_json(self.fixture.resource / RESOURCE_FILE)
        record["contract"]["files"] = file_records(self.fixture.resource)
        record["resource_id"] = identity(record["contract"])
        write_json(self.fixture.resource / RESOURCE_FILE, record)
        entry = self.fixture.entry
        entry["resource"] = record
        entry["resource_manifest_sha256"] = digest(
            self.fixture.resource / RESOURCE_FILE
        )
        entry["search_method"]["resource_id"] = record["resource_id"]
        entry["method_id"] = identity(
            {
                "search": entry["search_method"],
                "interpretation": entry["interpretation"],
            }
        )
        self.fixture.batchdir, self.fixture.batch = prepare_batches(
            self.fixture.fixture.bundles,
            entry["search_method"],
            self.root / "recovery-inputs",
            target_fasta_bytes=100_000,
        )[0]
        shutil.rmtree(self.native.task)
        self.fixture.plan("run-task")
        (self.native.work / "bundle").unlink()
        (self.native.work / "bundle").symlink_to(self.fixture.batchdir)
        executable = self.native.work / "bin/emapper.py"
        executable.write_text(
            executable.read_text().replace(
                f"if Path({str(self.native.failure)!r}).exists() and",
                f"if any(p['query'].startswith('B__') for p in proteins) and Path({str(self.native.failure)!r}).exists() and",
            )
        )
        self.native.failure.write_text("no_search\n")
        self.failed_raw = self.native.run_native(expected_exit=9)
        self.fixture.complete("failed", raw=self.failed_raw)
        self.packet = (
            self.root / "failed/annotation_batches" / self.fixture.batch["batch_id"]
        )
        self.members = [
            self.root / "failed/samples" / accession / "annotation/eggnog"
            for accession in ("A", "B", "C")
        ]
        self.before = inventory(self.root / "failed")
        self.recovery = self.root / "recovery"
        self.recovery.mkdir()
        for name, target in (
            ("task", self.native.task),
            ("bundle", self.fixture.batchdir),
            ("resource", self.fixture.resource),
        ):
            (self.recovery / name).symlink_to(target)
        self.native.failure.unlink()

    def recover(self, expected_exit=0):
        result = subprocess.run(
            [
                sys.executable,
                str(ROOT / "bin/recover_eggnog_batch.py"),
                "--task",
                "task",
                "--batch",
                "bundle",
                "--resource",
                "resource",
                "--previous",
                str(self.packet),
                "--members",
                *map(str, self.members),
            ],
            cwd=self.recovery,
            env=self.native.environment,
            text=True,
            capture_output=True,
            check=False,
            timeout=30,
        )
        self.assertEqual(
            result.returncode, expected_exit, result.stdout + result.stderr
        )
        return result

    def reseal(self):
        packet = read_json(self.packet / "batch_result.json")
        packet["raw_files"] = inventory(self.packet / "raw")
        packet["batch_result_id"] = identity(
            {key: value for key, value in packet.items() if key != "batch_result_id"}
        )
        write_json(self.packet / "batch_result.json", packet)
        for member in self.members:
            record = read_json(member / "result.json")
            record["raw_files"] = packet["raw_files"]
            record["native_batch"]["batch_result_id"] = packet["batch_result_id"]
            record["result_id"] = identity(
                {key: value for key, value in record.items() if key != "result_id"}
            )
            write_json(member / "result.json", record)

    def test_recovery_preserves_search_success_and_source_bytes(self):
        self.recover()
        calls = self.native.calls()
        self.assertEqual(
            [row["mode"] for row in calls],
            ["diamond", "no_search", "no_search", "no_search", "no_search"],
        )
        self.assertEqual(
            [[query.split("__")[0] for query in row["queries"]] for row in calls[1:]],
            [["A"], ["B"], ["B"], ["C"]],
        )
        raw = self.recovery / "raw"
        self.assertEqual(inventory(self.root / "failed"), self.before)
        for directory in ("search", "annotations/member00000000"):
            self.assertEqual(
                inventory(raw / directory), inventory(self.packet / "raw" / directory)
            )
        receipt = read_json(raw / "recovery.json")
        self.assertEqual(
            receipt["reused_stages"], ["search", "annotations/member00000000"]
        )
        self.assertEqual(
            receipt["executed_stages"],
            ["annotations/member00000001", "annotations/member00000002"],
        )
        self.assertEqual(set(receipt["database_files"]), set(ANNOTATION_FILES))
        self.assertTrue(
            (self.recovery / "eggnog_data/eggnog_proteins.dmnd").is_symlink()
        )
        archive = raw / "failed_attempts" / receipt["source_batch_result_id"]
        self.assertEqual(
            read_json(archive / "batch_result.json"),
            read_json(self.packet / "batch_result.json"),
        )
        self.assertEqual((archive / "exit_code.txt").read_text(), "9\n")
        self.assertEqual(
            (archive / "annotations/member00000001/exit_code.txt").read_text(), "9\n"
        )
        results = self.fixture.complete("recovered", raw=raw)["results"]
        self.assertEqual({record["status"] for record in results.values()}, {"success"})
        self.assertEqual(
            results["A"]["evidence"]["features"]["Preferred_name"],
            {"A::gene_1": ["hflB"]},
        )

    def test_second_failure_stays_failed_and_does_not_run_later_members(self):
        self.native.failure.write_text("no_search\n")
        self.recover(expected_exit=9)
        self.assertEqual(len(self.native.calls()), 4)
        self.assertEqual((self.recovery / "raw/exit_code.txt").read_text(), "9\n")
        results = self.fixture.complete("again-failed", raw=self.recovery / "raw")[
            "results"
        ]
        self.assertEqual({record["status"] for record in results.values()}, {"failed"})

    def test_incomplete_success_is_rejected_before_native_execution(self):
        path = self.packet / "raw/annotations/member00000000/eggnog.emapper.orthologs"
        path.write_text(path.read_text().split("##")[0])
        self.reseal()
        result = self.recover(expected_exit=1)
        self.assertIn("ortholog evidence", result.stderr)
        self.assertEqual(len(self.native.calls()), 3)

    def test_corrupt_search_is_rejected_before_native_execution(self):
        path = self.packet / "raw/search/eggnog.emapper.seed_orthologs"
        path.write_text(
            path.read_text().replace("## 3 queries scanned", "## 2 queries scanned")
        )
        self.reseal()
        result = self.recover(expected_exit=1)
        self.assertIn("completion count", result.stderr)
        self.assertEqual(len(self.native.calls()), 3)

    def test_coordinate_change_cannot_reuse_other_member_outputs(self):
        path = self.members[0] / "result.json"
        record = read_json(path)
        record["coordinate_id"] = "0" * 64
        record["result_id"] = identity(
            {key: value for key, value in record.items() if key != "result_id"}
        )
        write_json(path, record)
        result = self.recover(expected_exit=1)
        self.assertIn("member identity", result.stderr)
        self.assertEqual(len(self.native.calls()), 3)

    def test_resource_corruption_is_rejected_before_native_execution(self):
        path = self.fixture.resource / "eggnog.db.taxids.bin"
        path.write_bytes(b"X" * path.stat().st_size)
        result = self.recover(expected_exit=1)
        self.assertIn("database while staging", result.stderr)
        self.assertEqual(len(self.native.calls()), 3)

    def test_failed_sort_is_archived_and_rebuilt_from_original_seed_partition(self):
        sorted_seeds = self.packet / "raw/derived/member00000001/seeds.tsv.sorted"
        sorted_seeds.write_text("truncated failed sort\n")
        self.reseal()
        self.recover()
        raw = self.recovery / "raw"
        source_id = read_json(raw / "recovery.json")["source_batch_result_id"]
        self.assertEqual(
            (
                raw
                / "failed_attempts"
                / source_id
                / "derived/member00000001/seeds.tsv.sorted"
            ).read_text(),
            "truncated failed sort\n",
        )
        self.assertNotEqual(
            (raw / "derived/member00000001/seeds.tsv.sorted").read_text(),
            "truncated failed sort\n",
        )
        self.assertEqual(len(self.native.calls()), 5)

    def test_repeat_in_the_same_task_directory_preserves_the_finished_attempt(self):
        self.recover()
        before = inventory(self.recovery / "raw")
        result = self.recover(expected_exit=1)
        self.assertIn("fresh task directory", result.stderr)
        self.assertEqual(inventory(self.recovery / "raw"), before)
        self.assertEqual(len(self.native.calls()), 5)


if __name__ == "__main__":
    unittest.main()
