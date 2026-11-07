"""Real Nextflow publication recovery with synthetic, counted native tools."""

from __future__ import annotations

import shutil
import subprocess
import unittest
from pathlib import Path

import test_annotation_bundle as bundle_fixture
import test_eggnog_batch_integration as batch_fixture
from annotation_common import read_json, read_tsv
from annotation_result import inventory
from annotation_source import validate_source

ROOT = Path(__file__).resolve().parents[1]
# Recover a failure produced by the original release's execution code.
BASE = "6da8a4fb520076e8f88757137a31cdaedd87f963"
ROUTING_FILES = (
    "modules/local/functional_annotation.nf",
    "subworkflows/local/functional_annotation.nf",
    "nextflow.config",
    "reannotate.nf",
)


@unittest.skipUnless(
    batch_fixture.NEXTFLOW, "Nextflow is required for recovery integration"
)
class EggnogRecoveryIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.batch = batch_fixture.EggnogBatchIntegrationTests()
        self.batch.setUp()
        self.addCleanup(self.batch.doCleanups)
        self.root, self.project = self.batch.root, self.batch.project
        fixture = bundle_fixture.AnnotationBundleTests()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        self.source = self.batch.input_source(
            "six-members",
            list("ABCDEF"),
            {name: fixture.build(name, f"extra-{name}") for name in "CDEF"},
        )
        # ABCD share the failed batch: A/B passed, C fails, D never started.
        # EF form a separate successful batch; KOfam succeeds for all six genomes.
        executable = self.project / "bin/emapper.py"
        text = executable.read_text().replace(
            f"if Path({str(self.batch.fail_native)!r}).exists() and",
            f"if any(p['query'].startswith('C__') for p in proteins) and Path({str(self.batch.fail_native)!r}).exists() and",
        )
        executable.write_text(text)
        # This fake container executes on macOS, where Linux stat syntax is absent.
        stat = self.project / "bin/stat"
        stat.write_text(
            "#!/bin/sh\nif [ \"$*\" = '-f -c %T .' ]; then printf 'overlayfs\\n'; else exec /usr/bin/stat \"$@\"; fi\n"
        )
        stat.chmod(0o755)
        self.config = self.root / "recovery.config"
        self.config.write_text(
            """params {
 annotation_tools = 'eggnog,kofam'
 eggnog_container = 'sha256:EGG'
 python_container = 'sha256:HELPER'
 kofam_container = 'sha256:KOFAM'
 annotation_cpus = 1
 annotation_memory = 18.GB
 max_cpus = 1
 max_memory = 18.GB
 eggnog_group_samples = 4
 eggnog_db = 'EGGDB'
 kofam_db = 'KODB'
}
""".replace("EGGDB", str(self.batch.database))
            .replace("KODB", str(self.batch.fixture.database))
            .replace("EGG", "1" * 64)
            .replace("HELPER", "2" * 64)
            .replace("KOFAM", "3" * 64)
        )
        self.old = self.root / "failed"
        self.work = self.root / "work"

    def run_workflow(self, name, *, recovery=False, succeeds=True, extra=()):
        destination = self.root / name
        command = [
            batch_fixture.NEXTFLOW,
            "run",
            str(self.project / "reannotate.nf"),
            "-profile",
            "test,docker",
            "-c",
            str(self.config),
            "--annotation_from",
            str(self.source),
            "--outdir",
            str(destination),
            "-work-dir",
            str(self.work),
        ]
        if recovery:
            command += ["--eggnog_recovery_from", str(self.old)]
        command += list(extra)
        result = subprocess.run(
            command,
            cwd=self.root,
            env=self.batch.fixture.environment,
            capture_output=True,
            check=False,
            text=True,
            timeout=180,
        )
        (self.root / f"{name}.log").write_text(result.stdout + result.stderr)
        self.assertEqual(
            result.returncode == 0, succeeds, (result.stdout + result.stderr)[-9000:]
        )
        return destination, result.stdout + result.stderr

    def original_failure(self):
        # Start with the unmodified release's scripts, then apply only this fix.
        for name in ROUTING_FILES:
            original = subprocess.check_output(
                ["git", "show", f"{BASE}:{name}"], cwd=ROOT
            )
            (self.project / name).write_bytes(original)
        self.batch.fail_native.write_text("no_search\n")
        self.run_workflow("failed", succeeds=False)
        self.batch.fail_native.unlink()
        upgrade = self.root / "pipeline-upgrade"
        shutil.copytree(self.project, upgrade)
        self.project = upgrade
        for name in ROUTING_FILES:
            shutil.copyfile(ROOT / name, self.project / name)
        self.protected = inventory(self.old)
        self.previous_calls = self.batch.calls("no_search")
        self.searches = self.batch.calls()
        # No original task outputs/cache are needed to reuse published successes.
        shutil.rmtree(self.work)
        with self.config.open("a") as handle:
            handle.write(
                "process { withName: ANNOTATION_SEARCH { beforeScript = 'echo UNEXPECTED_NATIVE_SEARCH >&2; exit 97' } }\n"
            )

    def test_released_run_recovers_without_rerunning_successful_native_tasks(self):
        self.original_failure()
        output, log = self.run_workflow("recovered", recovery=True)
        manifest, batches = validate_source(output)
        self.assertTrue(manifest["complete"])
        self.assertEqual(inventory(self.old), self.protected)
        self.assertEqual(self.batch.calls(), self.searches)
        new = self.batch.calls("no_search")[len(self.previous_calls) :]
        self.assertEqual(
            [[query.split("__")[0] for query in queries] for queries in new],
            [["C"], ["D"]],
        )
        traces = read_tsv(output / "pipeline_info/trace.tsv")
        native = [row for row in traces if ":ANNOTATION_SEARCH " in row["name"]]
        self.assertEqual(native, [])
        plans = read_tsv(output / "pipeline_info/annotation_plan.tsv")
        self.assertEqual(
            {row["action"] for row in plans if row["tool"] == "kofam"}, {"reuse"}
        )
        self.assertEqual(
            {
                row["action"]
                for row in plans
                if row["tool"] == "eggnog" and row["accession"] in ("E", "F")
            },
            {"reuse"},
        )
        recovered = [row for row in traces if ":RECOVER_EGGNOG_BATCH " in row["name"]]
        self.assertEqual(len(recovered), 1)
        self.assertEqual(recovered[0]["status"], "COMPLETED")
        receipts = [
            read_json(batch.root / "raw/recovery.json")
            for batch in batches.values()
            if (batch.root / "raw/recovery.json").exists()
        ]
        self.assertEqual(len(receipts), 1)
        self.assertEqual(
            receipts[0]["reused_stages"],
            ["search", "annotations/member00000000", "annotations/member00000001"],
        )
        self.assertIn("ANNOTATION_ACCEPTANCE", log)

        # The healed publication reuses normally and then resumes genuinely cached.
        self.source = output
        before_annotations = self.batch.calls("no_search")
        self.run_workflow("reused")
        resumed, _ = self.run_workflow("reused", extra=("-resume",))
        self.assertTrue(validate_source(resumed)[0]["complete"])
        completions = [
            row
            for row in read_tsv(resumed / "pipeline_info/trace.tsv")
            if ":COMPLETE_EGGNOG_BATCH " in row["name"]
        ]
        self.assertEqual({row["status"] for row in completions}, {"CACHED"})
        self.assertEqual(self.batch.calls(), self.searches)
        self.assertEqual(self.batch.calls("no_search"), before_annotations)

    def test_changed_plan_stops_before_any_native_task(self):
        self.original_failure()
        _, log = self.run_workflow(
            "changed",
            recovery=True,
            succeeds=False,
            extra=("--annotation_cpus", "2", "--max_cpus", "2"),
        )
        self.assertIn("identical cohort", log)
        self.assertEqual(self.batch.calls(), self.searches)
        self.assertEqual(self.batch.calls("no_search"), self.previous_calls)
        traces = read_tsv(self.root / "changed/pipeline_info/trace.tsv")
        self.assertFalse(
            any(
                ":ANNOTATION_SEARCH " in row["name"]
                or ":RECOVER_EGGNOG_BATCH " in row["name"]
                for row in traces
            )
        )


if __name__ == "__main__":
    unittest.main()
