"""Verify native wrapper arguments independently of the synthetic cohort stubs."""

import json
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
class DomainToolCommandTests(unittest.TestCase):
    def test_both_domains_keep_explicit_codes_and_barrnap_models(self):
        with tempfile.TemporaryDirectory(prefix="domain-commands-", dir="/tmp") as name:
            root = Path(name)
            (root / "bin").mkdir()
            for tool in ("prokka", "barrnap"):
                shutil.copyfile(ROOT / f"modules/local/{tool}.nf", root / f"{tool}.nf")
            shutil.copy2(
                ROOT / "bin/annotation_storage.sh", root / "bin/annotation_storage.sh"
            )
            code = (
                "#!"
                + sys.executable
                + "\n"
                + r"""
import json
import sys
from pathlib import Path
if sys.argv[1:] == ['--version']:
    print('synthetic wrapper control')
    sys.exit(0)
args = sys.argv[1:]
tool = Path(sys.argv[0]).name
identifier = Path(args[0] if tool == 'prokka' else args[-1]).stem
expected = {'B4':('Bacteria','bac','4'), 'B11':('Bacteria','bac','11'),
            'A4':('Archaea','arc','4'), 'A11':('Archaea','arc','11')}[identifier]
kingdom = args[args.index('--kingdom') + 1]
if tool == 'prokka':
    genetic_code = args[args.index('--gcode') + 1]
    assert (kingdom, genetic_code) == (expected[0], expected[2])
    assert '--rfam' in args and '--compliant' in args
    target = Path(args[args.index('--outdir') + 1]); target.mkdir()
    for suffix in ('gff','faa','gbk'):
        (target / (identifier + '.' + suffix)).write_text(json.dumps([kingdom,genetic_code]) + '\n')
else:
    assert kingdom == expected[1]
    Path(args[args.index('--outseq') + 1]).write_text('>control\nACGT\n')
    print(kingdom)
"""
            )
            for tool in ("prokka", "barrnap"):
                path = root / "bin" / tool
                path.write_text(code)
                path.chmod(0o755)
            jobs = []
            cases = {
                "B4": ("Bacteria", "bac", 4),
                "B11": ("Bacteria", "bac", 11),
                "A4": ("Archaea", "arc", 4),
                "A11": ("Archaea", "arc", 11),
            }
            for accession, (kingdom, barrnap, gcode) in cases.items():
                path = root / f"{accession}.fasta"
                path.write_text(">contig\n" + "ACGT" * 100 + "\n")
                prokka = "" if kingdom == "Bacteria" else f",prokka_kingdom:'{kingdom}'"
                jobs.append(
                    f"tuple([accession:'{accession}',internal_id:'{accession}',barrnap_kingdom:'{barrnap}'{prokka}],file('{path}'),{gcode})"
                )
            (root / "main.nf").write_text(
                "include { PROKKA } from './prokka'\ninclude { BARRNAP } from './barrnap'\nworkflow {\n samples=Channel.of("
                + ",".join(jobs)
                + ")\n PROKKA(samples)\n BARRNAP(samples.map { item -> tuple(item[0],item[1]) })\n}\n"
            )
            (root / "nextflow.config").write_text(
                f"params {{ outdir='{root / 'results'}'; soft_fail_attempts=1; annotation_cleanup='success'; domain='metadata' }}\nprocess.executor='local'\n"
            )
            result = subprocess.run(
                [NEXTFLOW, "run", "main.nf"],
                cwd=root,
                env=dict(
                    os.environ, NXF_ANSI_LOG="false", NXF_DISABLE_CHECK_LATEST="true"
                ),
                capture_output=True,
                text=True,
                timeout=180,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            for accession, (kingdom, barrnap, gcode) in cases.items():
                sample = root / "results/samples" / accession
                self.assertEqual(
                    json.loads((sample / "prokka/prokka.gff").read_text()),
                    [kingdom, str(gcode)],
                )
                self.assertEqual(
                    (sample / "barrnap/rrna.gff").read_text().strip(), barrnap
                )
