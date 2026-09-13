"""Fail-closed policy selection and the retired Nextflow flag migration."""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

from ani_common import derive_sixteen_s_ani_exclusion_reason


class Ani16sPolicyTests(unittest.TestCase):
    def test_invalid_policy_or_status_does_not_become_eligible(self):
        for policy in (True, False, None, "", "all", "IGNORE"):
            with (
                self.subTest(policy=policy),
                self.assertRaisesRegex(ValueError, "policy"),
            ):
                derive_sixteen_s_ani_exclusion_reason("Yes", policy=policy)
        with self.assertRaisesRegex(ValueError, "status"):
            derive_sixteen_s_ani_exclusion_reason("corrupt", policy="ignore")
        for absent in (None, "", "NA"):
            self.assertIsNone(
                derive_sixteen_s_ani_exclusion_reason(absent, policy="ignore")
            )
            self.assertEqual(derive_sixteen_s_ani_exclusion_reason(absent), "16s_na")

    @unittest.skipUnless(shutil.which("nextflow"), "Nextflow is required")
    def test_removed_flag_and_unknown_policy_fail_before_analysis(self):
        cases = (
            (["--ani_allow_incomplete_16s"], "has been removed; use --ani_16s_policy"),
            (
                ["--ani_16s_policy", "everything"],
                "must be complete, allow_incomplete, or ignore",
            ),
        )
        for arguments, message in cases:
            with (
                self.subTest(arguments=arguments),
                tempfile.TemporaryDirectory(dir="/tmp") as name,
            ):
                result = subprocess.run(
                    [
                        shutil.which("nextflow"),
                        "run",
                        str(ROOT / "main.nf"),
                        "-profile",
                        "test",
                        *arguments,
                    ],
                    cwd=name,
                    check=False,
                    capture_output=True,
                    text=True,
                    timeout=120,
                    env=dict(
                        os.environ,
                        NXF_ANSI_LOG="false",
                        NXF_DISABLE_CHECK_LATEST="true",
                    ),
                )
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(message, result.stdout + result.stderr)
                self.assertNotIn("Submitted process", result.stdout + result.stderr)
