"""Distinguish a successful 16S-negative Barrnap run from native failure."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")


@unittest.skipUnless(NEXTFLOW, "Nextflow is required for Barrnap integration")
class BarrnapRuntimeTests(unittest.TestCase):
    def test_node_scratch_and_native_exit_are_explicit(self) -> None:
        with tempfile.TemporaryDirectory(prefix="nf-barrnap-", dir="/tmp") as temp:
            root = Path(temp)
            project = root / "project"
            module = project / "modules/local"
            module.mkdir(parents=True)
            shutil.copy2(ROOT / "modules/local/barrnap.nf", module / "barrnap.nf")
            (project / "main.nf").write_text(
                "nextflow.enable.dsl = 2\n"
                "include { BARRNAP } from './modules/local/barrnap'\n"
                "workflow { BARRNAP(Channel.of(tuple([accession: 'A'], "
                "file(params.genome, checkIfExists: true)))) }\n"
            )
            (project / "nextflow.config").write_text(
                "params.soft_fail_attempts = 1\n"
                "params.barrnap_container = null\n"
                "params.barrnap_kingdom = 'bac'\n"
                f"process.scratch = '{root / 'node_scratch'}'\n"
                "process.stageOutMode = 'copy'\n"
                "process.executor = 'local'\n"
                "process.errorStrategy = 'finish'\n"
                "process.maxRetries = 0\n"
            )
            (root / "node_scratch").mkdir()
            genome = root / "genome.fasta"
            genome.write_text(">A\nACGTACGTACGT\n")
            bindir = root / "bin"
            bindir.mkdir()
            caller = bindir / "barrnap"
            caller.write_text(
                "#!/usr/bin/env python3\n"
                "import os, sys\n"
                "from pathlib import Path\n"
                "if '--version' in sys.argv:\n"
                "    print('barrnap synthetic')\n"
                "    raise SystemExit(0)\n"
                "temp = Path(os.environ['TMPDIR'])\n"
                "Path(os.environ['FAKE_BARRNAP_DIAGNOSTIC']).write_text(\n"
                "    f'temp={temp} cwd={Path.cwd()}\\n')\n"
                "if temp.parent.resolve() != Path.cwd().resolve():\n"
                "    raise SystemExit('Barrnap did not receive task-owned TMPDIR')\n"
                "(temp / 'native-temp').write_text('created')\n"
                "mode = os.environ['FAKE_BARRNAP_MODE']\n"
                "if mode == 'failure':\n"
                "    print('Synthetic Barrnap failure', file=sys.stderr)\n"
                "    raise SystemExit(30)\n"
                "Path(sys.argv[sys.argv.index('--outseq') + 1]).write_text(\n"
                "    '>16S\\nACGT\\n' if mode == 'positive' else '')\n"
                "if mode == 'positive':\n"
                "    print('A\\tbarrnap\\trRNA\\t1\\t4\\t.\\t+\\t.\\tName=16S_rRNA')\n"
            )
            caller.chmod(0o755)

            for mode in ("positive", "negative", "failure"):
                output = root / f"results-{mode}"
                environment = dict(
                    os.environ,
                    PATH=str(bindir) + os.pathsep + os.environ["PATH"],
                    FAKE_BARRNAP_MODE=mode,
                    FAKE_BARRNAP_DIAGNOSTIC=str(root / f'diagnostic-{mode}.txt'),
                    NXF_ANSI_LOG="false",
                    NXF_DISABLE_CHECK_LATEST="true",
                    NXF_OFFLINE="true",
                )
                result = subprocess.run(
                    [NEXTFLOW, "run", str(project / "main.nf"),
                     "--genome", str(genome), "--outdir", str(output),
                     "-work-dir", str(root / f"work-{mode}")],
                    cwd=root, env=environment, text=True, capture_output=True,
                    timeout=120,
                )
                published = output / "samples/A/barrnap"
                if mode == "failure":
                    self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
                    self.assertFalse((published / "rrna.gff").exists())
                    self.assertFalse((published / "rrna.fa").exists())
                else:
                    diagnostics = '\n'.join(
                        path.read_text(errors='replace')
                        for path in (root / f'work-{mode}').glob('*/*/barrnap.log')
                    )
                    diagnostic = root / f'diagnostic-{mode}.txt'
                    diagnostics += diagnostic.read_text() if diagnostic.exists() else 'fake caller did not start'
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr + diagnostics)
                    self.assertEqual((published / "rrna.gff").stat().st_size > 0,
                                     mode == "positive")
                    self.assertEqual((published / "rrna.fa").stat().st_size > 0,
                                     mode == "positive")
                    self.assertIn("exit_code=0", (published / "barrnap.log").read_text())


if __name__ == "__main__":
    unittest.main()
