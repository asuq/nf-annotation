"""Exercise real Nextflow Codetta staging with input-driven fake callers."""

import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")


@unittest.skipUnless(NEXTFLOW, "Nextflow is required")
class CodettaStagingTests(unittest.TestCase):
    def exercise(self, symlink=False):
        with tempfile.TemporaryDirectory(prefix="codetta-staging-", dir="/tmp") as name:
            root = Path(name)
            project = root / "project"
            (project / "bin").mkdir(parents=True)
            shutil.copyfile(ROOT / "modules/local/codetta.nf", project / "codetta.nf")
            database = root / "codetta"
            database.mkdir()
            (database / "Pfam-A_enone.hmm").write_text("reference\n")
            home = root / "runtime"
            (home / "resources").mkdir(parents=True)
            (home / "resources/bad_pfams.txt").write_text("resource\n")
            genome = root / "genome.fasta"
            genome.write_text(">contig\nATGAAATAA\n")
            outside = root / "outside"
            outside.mkdir()
            caller = (
                "#!"
                + sys.executable
                + "\n"
                + r"""
import json
import sys
from pathlib import Path
args = sys.argv[1:]
resource = Path(args[args.index('--resource_directory') + 1])
assert (resource / 'Pfam-A_enone.hmm').read_text() == 'reference\n'
assert (resource / 'bad_pfams.txt').read_text() == 'resource\n'
with Path('codetta/calls.jsonl').open('a') as handle:
    handle.write(json.dumps([Path(sys.argv[0]).name, args]) + '\n')
if '--inference_output' in args:
    Path(args[args.index('--inference_output') + 1]).write_text('inference fixture\n')
    Path(args[args.index('--results_summary') + 1]).write_text('summary fixture\n')
"""
            )
            for script in (
                "codetta_align.py",
                "codetta_summary.py",
                "codetta_infer.py",
            ):
                path = project / "bin" / script
                path.write_text(caller)
                path.chmod(0o755)
            (project / "main.nf").write_text(
                "include { CODETTA } from './codetta'\n"
                "workflow { CODETTA(Channel.of(tuple([accession:'S',internal_id:'S'], "
                "file(params.genome), file(params.database)))) }\n"
            )
            (project / "nextflow.config").write_text(
                "params { soft_fail_attempts = 1; codetta_extra_args = ''; "
                f"outdir = '{root / 'results'}'; genome = '{genome}'; database = '{database}' }}\n"
                "process.executor = 'local'\n"
                + (
                    f"process.beforeScript = 'ln -s {outside} codetta'\n"
                    if symlink
                    else ""
                )
            )
            before = {
                p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                for p in database.iterdir()
            }
            result = subprocess.run(
                [
                    NEXTFLOW,
                    "run",
                    str(project / "main.nf"),
                    "-work-dir",
                    str(root / "work"),
                ],
                cwd=project,
                env=dict(
                    os.environ,
                    CODETTA_HOME=str(home),
                    NXF_ANSI_LOG="false",
                    NXF_DISABLE_CHECK_LATEST="true",
                ),
                capture_output=True,
                text=True,
                timeout=180,
                check=False,
            )
            after = {
                p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                for p in database.iterdir()
            }
            self.assertEqual(before, after, result.stdout + result.stderr)
            self.assertEqual(list(outside.iterdir()), [])
            if symlink:
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(
                    "Codetta output must be an ordinary directory",
                    result.stdout + result.stderr,
                )
            else:
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                output = root / "results/samples/S/codetta"
                self.assertFalse(output.is_symlink())
                self.assertEqual(
                    (output / "codetta_inference.txt").read_text(),
                    "inference fixture\n",
                )
                self.assertEqual(
                    len((output / "calls.jsonl").read_text().splitlines()), 3
                )
                self.assertIn("exit_code=0", (output / "codetta.log").read_text())

    def test_database_named_codetta_is_separate_from_output(self):
        self.exercise()

    def test_output_symlink_is_rejected_before_native_writes(self):
        self.exercise(symlink=True)


if __name__ == "__main__":
    unittest.main()
