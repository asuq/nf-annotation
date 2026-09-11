"""Integration checks for stub-run resume safety after output pruning."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
from annotation_common import bundle_proteins

NEXTFLOW = shutil.which("nextflow")
TARGET_PROCESSES = (
    "BARRNAP",
    "CHECKM2",
    "BUSCO",
    "PROKKA",
    "CCFINDER",
    "CODETTA",
    "BUILD_COHORT_16S",
    "SUMMARISE_BUSCO",
    "SUMMARISE_CODETTA",
)


class StubResumeStorageTestCase(unittest.TestCase):
    """Verify that the pruned task outputs still support `-resume`."""

    def run_pipeline(
        self,
        *,
        run_dir: Path,
        outdir: Path,
        workdir: Path,
        resume: bool = False,
    ) -> subprocess.CompletedProcess[str]:
        """Run the stub pipeline with isolated Nextflow state."""
        command = [
            NEXTFLOW,
            "run",
            ".",
            "-profile",
            "test",
            "-stub-run",
            "--outdir",
            str(outdir),
            "-work-dir",
            str(workdir),
        ]
        if resume:
            command.append("-resume")

        environment = os.environ.copy()
        environment["NXF_ANSI_LOG"] = "false"
        environment["NXF_DISABLE_CHECK_LATEST"] = "true"
        environment["NXF_HOME"] = str(run_dir / ".nxf")
        return subprocess.run(
            command,
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
            env=environment,
        )

    def work_leaf_dirs(self, workdir: Path) -> set[str]:
        """Return the hashed leaf work directories created by Nextflow."""
        if not workdir.exists():
            return set()
        leaf_dirs: set[str] = set()
        for path in workdir.glob("*/*"):
            if path.is_dir():
                leaf_dirs.add(path.relative_to(workdir).as_posix())
        return leaf_dirs

    def task_work_leaf_dirs(self, workdir: Path) -> set[str]:
        """Return only real task work directories, excluding `work/tmp/*`."""
        return {path for path in self.work_leaf_dirs(workdir) if not path.startswith("tmp/")}

    @unittest.skipUnless(NEXTFLOW is not None, "nextflow is required for resume tests.")
    def test_stub_run_pruned_outputs_remain_resume_safe(self) -> None:
        """Reuse the same work cache only for the `-resume` run."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            outdir = tmpdir / "results"
            workdir = tmpdir / "work"

            first = self.run_pipeline(run_dir=tmpdir, outdir=outdir, workdir=workdir)
            self.assertEqual(first.returncode, 0, first.stderr)
            first_work_dirs = self.task_work_leaf_dirs(workdir)
            self.assertTrue(first_work_dirs)

            second = self.run_pipeline(run_dir=tmpdir, outdir=outdir, workdir=workdir)
            self.assertEqual(second.returncode, 0, second.stderr)
            second_work_dirs = self.task_work_leaf_dirs(workdir)
            self.assertGreater(len(second_work_dirs), len(first_work_dirs))

            resumed = self.run_pipeline(run_dir=tmpdir, outdir=outdir, workdir=workdir, resume=True)
            self.assertEqual(resumed.returncode, 0, resumed.stderr)
            resumed_work_dirs = self.task_work_leaf_dirs(workdir)
            self.assertEqual(resumed_work_dirs, second_work_dirs)
            self.assertTrue((outdir / "tables" / "master_table.tsv").is_file())
            bundle = outdir / "samples/TEST_ACC/annotation/bundle"
            manifest, proteins = bundle_proteins(bundle)
            self.assertEqual(manifest["genetic_code"], 4)
            self.assertEqual(len(proteins), 1)
            self.assertEqual((bundle / "proteins.faa").read_text().splitlines()[1], "MWA")

            resume_output = "\n".join((resumed.stdout, resumed.stderr))
            self.assertIn("cached", resume_output.lower())
            self.assertTrue(any(process in resume_output for process in TARGET_PROCESSES))
