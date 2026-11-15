"""Observe released tools without changing their arguments or output stream."""

import hashlib
import json
import os
import signal
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "bin/run_prokka.py"


class ProkkaExecutionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir="/tmp")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.tools = self.root / "tools"
        self.tools.mkdir()
        self.prodigal = self.tools / "prodigal"
        self.prodigal.write_text(
            f"#!{sys.executable}\n"
            + """import json,os,signal,sys
from pathlib import Path
if '-v' in sys.argv:
    print('Prodigal V2.6.3',file=sys.stderr)
    sys.exit(0)
Path('child_arguments.json').write_text(json.dumps(sys.argv[1:]))
print('unchanged prediction stream',flush=True)
code=int(os.environ['CHILD_EXIT'])
if code<0:os.kill(os.getpid(),-code)
sys.exit(code)
"""
        )
        self.prokka = self.tools / "prokka"
        self.prokka.write_text(
            f"#!{sys.executable}\n"
            + """import json,os,subprocess,sys
from pathlib import Path
Path('parent_arguments.json').write_text(json.dumps(sys.argv[1:]))
subprocess.run(['prodigal','-v'],check=False)
if not os.environ.get('OMIT_PREDICTION'):
    result=subprocess.run(['prodigal','-i','input.fna','-c','-m','-g','11','-p','single','-f','sco','-q'],capture_output=True,check=False)
    Path('received_stream').write_bytes(result.stdout)
# Mimic the released parent's masked failure and zero-CDS message.
print('[12:00:00] Found 0 CDS',flush=True)
sys.exit(int(os.environ.get('PARENT_EXIT','0')))
"""
        )
        for path in (self.prokka, self.prodigal):
            path.chmod(0o755)
        self.original = {
            p: hashlib.md5(p.read_bytes()).hexdigest()
            for p in (self.prokka, self.prodigal)
        }

    def run_case(self, child_exit=0, **environment):
        arguments = [
            "input.fna",
            "--outdir",
            "prokka",
            "--compliant",
            "--kingdom",
            "Bacteria",
            "--gcode",
            "11",
            "--cpus",
            "2",
            "--rfam",
        ]
        result = subprocess.run(
            [sys.executable, str(SCRIPT), *arguments],
            cwd=self.root,
            env=dict(
                os.environ,
                PATH=str(self.tools) + os.pathsep + os.environ["PATH"],
                CHILD_EXIT=str(child_exit),
                **environment,
            ),
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(
            json.loads((self.root / "parent_arguments.json").read_text()), arguments
        )
        self.assertEqual(
            self.original,
            {p: hashlib.md5(p.read_bytes()).hexdigest() for p in self.original},
        )
        return result, json.loads((self.root / "prokka_execution.json").read_text())

    def test_successful_zero_keeps_arguments_stream_and_released_binaries(self):
        result, record = self.run_case()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(
            (self.root / "received_stream").read_bytes(),
            b"unchanged prediction stream\n",
        )
        self.assertEqual(
            json.loads((self.root / "child_arguments.json").read_text()),
            [
                "-i",
                "input.fna",
                "-c",
                "-m",
                "-g",
                "11",
                "-p",
                "single",
                "-f",
                "sco",
                "-q",
            ],
        )
        self.assertEqual(record["child"]["returncode"], 0)
        self.assertEqual(record["status"], "success")
        self.assertEqual(record["cds_count"], 0)
        self.assertIn("NF_PROKKA_PREDICTION_V1\tsingle\t0", result.stdout)

    def test_masked_nonzero_and_signalled_children_remain_failures(self):
        for code in (7, -signal.SIGTERM, -signal.SIGKILL):
            with self.subTest(code=code):
                result, record = self.run_case(code)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(record["parent_returncode"], 0)
                self.assertEqual(record["child"]["returncode"], code)
                self.assertEqual(record["status"], "failed")
                self.assertIsNone(record["cds_count"])
                self.assertNotIn("NF_PROKKA_PREDICTION_V1", result.stdout)

    def test_parent_failure_or_missing_child_evidence_is_not_zero(self):
        for settings in ({"PARENT_EXIT": "5"}, {"OMIT_PREDICTION": "1"}):
            with self.subTest(settings=settings):
                result, record = self.run_case(**settings)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(record["status"], "failed")
                self.assertIsNone(record["cds_count"])
