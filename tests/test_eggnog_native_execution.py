"""Execute the actual batch wrapper against context-sensitive synthetic native tools."""

from __future__ import annotations

import copy
import json
import os
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import test_annotation_batch_tasks as fixture
from annotation_common import AnnotationError, digest, read_json, write_json
from annotation_result import inventory
from eggnog_batches import validate_batch
from eggnog_native import validate_execution
from test_eggnog_batch_integration import EMAPPER


class EggnogNativeExecutionTests(unittest.TestCase):
    def setUp(self):
        self.fixture = fixture.AnnotationBatchTaskTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root
        self.fixture.plan("run-task")
        self.work = self.root / "work"
        self.work.mkdir()
        self.task = self.root / "run-task"
        for name, source in (
            ("task", self.task),
            ("bundle", self.fixture.batchdir),
            ("resource", self.fixture.resource),
        ):
            (self.work / name).symlink_to(source, target_is_directory=True)
        executable_directory = self.work / "bin"
        executable_directory.mkdir()
        self.ledger, self.failure, self.no_hits = (
            self.root / name for name in ("calls.jsonl", "fail", "no-hits.json")
        )
        executable = executable_directory / "emapper.py"
        executable.write_text(
            f"#!{sys.executable}\n"
            + EMAPPER.replace("__LEDGER__", repr(str(self.ledger)))
            .replace("__FAIL__", repr(str(self.failure)))
            .replace("__NOHITS__", repr(str(self.no_hits)))
        )
        executable.chmod(0o755)
        diamond = executable_directory / "diamond"
        diamond.write_text("#!/bin/sh\nprintf 'Synthetic DIAMOND control\\n'\n")
        diamond.chmod(0o755)
        self.environment = dict(
            os.environ,
            PATH=str(executable_directory) + os.pathsep + os.environ["PATH"],
            PYTHONPATH=str(ROOT / "bin"),
        )
        self.batch, self.proteins = validate_batch(self.fixture.batchdir)

    def run_native(self, expected_exit=0):
        result = subprocess.run(
            ["bash", "task/run.sh"],
            cwd=self.work,
            env=self.environment,
            capture_output=True,
            text=True,
            check=False,
            timeout=30,
        )
        self.assertEqual(
            result.returncode, expected_exit, result.stdout + result.stderr
        )
        return self.work / "raw"

    def calls(self):
        return [json.loads(line) for line in self.ledger.read_text().splitlines()]

    def complete(self, raw):
        return self.fixture.complete("completed", raw=raw)["results"]

    def test_one_search_and_separate_annotations_preserve_input_bytes_and_names(self):
        before = inventory(self.task)
        raw = self.run_native()
        self.assertEqual((raw / "exit_code.txt").read_text(), "0\n")
        self.assertEqual(inventory(self.task), before)
        self.assertFalse((self.task / "__pycache__").exists())
        calls = self.calls()
        self.assertEqual(
            [row["mode"] for row in calls],
            ["diamond", "no_search", "no_search", "no_search"],
        )
        self.assertEqual([len(row["queries"]) for row in calls], [3, 1, 1, 1])
        self.assertFalse((raw / "search/eggnog.emapper.annotations").exists())
        results = self.complete(raw)
        for protein in self.proteins:
            record = results[protein["accession"]]
            self.assertEqual(record["status"], "success")
            self.assertEqual(
                record["evidence"]["features"]["Preferred_name"],
                {protein["gene_id"]: ["hflB"]},
            )

    def test_empty_seed_member_runs_native_annotation_and_keeps_valid_zero(self):
        protein = next(row for row in self.proteins if row["accession"] == "C")
        write_json(self.no_hits, [protein["tool_id"]])
        raw = self.run_native()
        self.assertEqual([row["mode"] for row in self.calls()].count("no_search"), 3)
        self.assertEqual(
            (raw / "derived/member00000002/seeds.tsv.sorted").read_bytes(), b""
        )
        record = self.complete(raw)["C"]
        self.assertEqual(record["status"], "success")
        self.assertEqual(record["reported_hits"], 0)
        self.assertEqual(record["accepted_proteins"], 0)
        self.assertEqual(record["evidence"]["mapped_genes"], [])

    def test_annotation_failure_stops_later_phases_and_fails_the_whole_batch(self):
        self.failure.write_text("no_search\n")
        raw = self.run_native(expected_exit=9)
        self.assertEqual((raw / "exit_code.txt").read_text(), "9\n")
        stages = read_json(raw / "execution.json")["stages"]
        self.assertEqual([stage["exit_code"] for stage in stages], [0, 9])
        self.assertEqual(
            [row["mode"] for row in self.calls()], ["diamond", "no_search"]
        )
        records = self.complete(raw)
        self.assertEqual(set(records), {"A", "B", "C"})
        self.assertEqual({record["status"] for record in records.values()}, {"failed"})

    def test_changed_input_fails_before_search_or_annotation(self):
        path = self.task / "input.faa"
        path.write_bytes(path.read_bytes().replace(b"MWA", b"MWF"))
        raw = self.run_native(expected_exit=1)
        self.assertEqual((raw / "exit_code.txt").read_text(), "1\n")
        self.assertFalse(self.ledger.exists())
        self.assertIn("input differs", (raw / "tool.log").read_text())

    def test_native_phase_commands_and_status_cannot_be_relabelled(self):
        raw = self.run_native()
        original = read_json(raw / "execution.json")
        for change in (
            lambda record: record["stages"].pop(),
            lambda record: record["stages"][1].update(exit_code=False),
            lambda record: record["stages"][1]["command"].append("--no_annot"),
        ):
            record = copy.deepcopy(original)
            change(record)
            write_json(raw / "execution.json", record)
            with self.assertRaises(AnnotationError):
                validate_execution(raw, self.batch, self.proteins)

    def test_seed_partition_swap_fails_even_with_a_changed_local_checksum(self):
        raw = self.run_native()
        target = raw / "derived/member00000000"
        (target / "seeds.tsv").write_bytes(
            (raw / "derived/member00000001/seeds.tsv").read_bytes()
        )
        receipt = read_json(target / "partition.json")
        receipt["derived_seed_sha256"] = digest(target / "seeds.tsv")
        write_json(target / "partition.json", receipt)
        with self.assertRaisesRegex(AnnotationError, "partition"):
            validate_execution(raw, self.batch, self.proteins)


if __name__ == "__main__":
    unittest.main()
