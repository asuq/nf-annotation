"""Validate ordered annotation manifests and execute their bounded shell argv."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

from annotation_common import AnnotationError, read_json, write_json, write_tsv
from annotation_path_lists import read_path_list
from annotation_tasks import preflight


class AnnotationPathListTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(
            prefix="annotation-path-lists-", dir="/tmp"
        )
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.manifest = self.root / "input paths.json"

    def test_empty_and_ordered_literal_paths(self):
        for values in (
            [],
            ["z path", "a'quote", 'double"quote', "back\\slash", "line\nbreak", " "],
        ):
            with self.subTest(values=values):
                write_json(self.manifest, values)
                self.assertEqual(
                    read_path_list(self.manifest), [Path(p) for p in values]
                )

    def test_malformed_arrays_and_duplicate_paths_fail(self):
        for values in (
            None,
            "bundle",
            {"paths": ["bundle"]},
            [""],
            [None],
            [1],
            [True],
            [["bundle"]],
            [{}],
            ["a\0b"],
            ["bundle", "bundle"],
            ["bundle", "./bundle"],
        ):
            with self.subTest(values=values):
                write_json(self.manifest, values)
                with self.assertRaises(AnnotationError):
                    read_path_list(self.manifest)
        self.manifest.write_text('["truncated"')
        with self.assertRaisesRegex(AnnotationError, "Invalid JSON"):
            read_path_list(self.manifest)

    def planner_cli(self, *arguments):
        return subprocess.run(
            [
                sys.executable,
                str(ROOT / "bin/prepare_annotation_tasks.py"),
                "plan",
                *arguments,
            ],
            cwd=self.root,
            text=True,
            capture_output=True,
            timeout=30,
        )

    def test_planner_cli_accepts_explicit_empty_array_and_preserves_sample_order(self):
        samples, receipt = self.root / "samples.tsv", self.root / "receipt.json"
        write_tsv(
            samples, ("accession",), [{"accession": acc} for acc in ("B", "A")]
        )
        write_json(receipt, preflight({"annotation_tools": ""}))
        write_json(self.manifest, [])
        output = self.root / "planned"
        result = self.planner_cli(
            "--samples", str(samples), "--preflight", str(receipt),
            "--bundle-list", str(self.manifest), "--output", str(output),
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        plan = read_json(output / "annotation_plan.json")
        self.assertEqual(plan["accessions"], ["B", "A"])
        self.assertEqual(
            [row["accession"] for row in plan["tasks"]], ["B"] * 5 + ["A"] * 5
        )
        self.assertEqual(
            {row["status"] for row in plan["tasks"]}, {"skipped_disabled"}
        )

    def test_planner_cli_requires_new_flag_and_rejects_invalid_manifest(self):
        common = [
            "--samples", "samples.tsv", "--preflight", "receipt.json",
            "--output", "planned",
        ]
        for arguments in (common, [*common, "--bundle", "bundle"]):
            result = self.planner_cli(*arguments)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("--bundle-list", result.stderr)
        write_json(self.manifest, ["duplicate", "duplicate"])
        result = self.planner_cli(*common, "--bundle-list", str(self.manifest))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Duplicate path", result.stderr)
        self.assertFalse((self.root / "planned").exists())

    def render_shell(self, process: str, bundles: list[str], results: list[str]) -> str:
        """Interpolate the real module's shell template without a Nextflow runtime."""
        module = (ROOT / "modules/local/functional_annotation.nf").read_text()
        block = module.split(f"process {process} {{\n", 1)[1]
        block = block.split("\nprocess ", 1)[0]
        block = block.split("\n    script:\n", 1)[1]
        shell = block.split('    """\n', 1)[1].rsplit('\n    """', 1)[0]
        # Groovy removes the escaped dollar before interpolating the path JSON.
        shell = shell.replace("\\$", "$")
        replacements = {
            "bundleList": json.dumps(bundles),
            "resultList": json.dumps(results),
            "samples": "samples.tsv",
            "receipt": "receipt.json",
            "sourceArgs": "",
            "plan": "annotation_plan.json",
        }
        for name, value in replacements.items():
            shell = shell.replace("${" + name + "}", value)
        return "set -euo pipefail\n" + shell + "\n"

    def execute_shell(
        self, process: str, bundles: list[str], results: list[str], name: str
    ):
        work = self.root / name
        work.mkdir()
        executables = work / "bin"
        executables.mkdir()
        helper = executables / (
            "prepare_annotation_tasks.py"
            if process == "PLAN_ANNOTATIONS"
            else "aggregate_annotations.py"
        )
        helper.write_text(
            "import json, sys\nfrom pathlib import Path\n"
            "Path('argv.json').write_text(json.dumps(sys.argv[1:]))\n"
            "if sys.argv[1] == 'plan':\n"
            " Path('planned').mkdir()\n"
            " Path('planned/annotation_plan.json').write_text('{}')\n"
            " Path('planned/annotation_plan.tsv').write_text('accession\\n')\n"
        )
        helper.chmod(0o755)
        script = work / "run.sh"
        script.write_text(self.render_shell(process, bundles, results))
        result = subprocess.run(
            ["bash", str(script)],
            cwd=work,
            text=True,
            capture_output=True,
            timeout=30,
            env={
                **os.environ,
                "PATH": str(executables) + os.pathsep + os.environ["PATH"],
            },
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(read_json(work / "bundle_list.json"), bundles)
        if process == "AGGREGATE_ANNOTATIONS":
            self.assertEqual(read_json(work / "result_list.json"), results)
        self.assertFalse((work / "INJECTED_DOLLAR").exists())
        self.assertFalse((work / "INJECTED_BACKTICK").exists())
        return read_json(work / "argv.json")

    def test_module_heredocs_preserve_literals_and_keep_argv_bounded(self):
        literals = [
            "z path", "a'quote", 'double"quote', "back\\slash", "line\nbreak",
            "$(touch INJECTED_DOLLAR)", "`touch INJECTED_BACKTICK`",
        ]
        bundles = [f"bundles/bundle{i:02d}" for i in range(10000)]
        results = [f"results/result{i:03d}" for i in range(50000)]
        for process in ("PLAN_ANNOTATIONS", "AGGREGATE_ANNOTATIONS"):
            with self.subTest(process=process):
                empty = self.execute_shell(process, [], [], process + "-empty")
                small = self.execute_shell(
                    process, literals, list(reversed(literals)), process + "-literal"
                )
                large = self.execute_shell(process, bundles, results, process + "-large")
                self.assertEqual(empty, small)
                self.assertEqual(small, large)
                self.assertLess(sum(len(value.encode()) + 1 for value in large), 1024)
                self.assertIn("--bundle-list", large)
                self.assertNotIn("--bundle", large)
                self.assertNotIn("--result", large)


if __name__ == "__main__":
    unittest.main()
