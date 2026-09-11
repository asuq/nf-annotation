"""Container contract tests for the custom PADLOC runtime image."""

from __future__ import annotations

import os
import subprocess
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCKERFILE = ROOT / "docker" / "padloc" / "Dockerfile"
IMAGE = "codex-padloc:test"
RUN_DOCKER_TESTS = os.environ.get("RUN_DOCKER_TESTS") == "1"


def run_command(
    args: list[str], *, timeout: int | None = None
) -> subprocess.CompletedProcess[str]:
    """Run one subprocess and return the completed process."""
    return subprocess.run(
        args,
        check=False,
        capture_output=True,
        text=True,
        cwd=ROOT,
        timeout=timeout,
    )


class PadlocContainerContractTestCase(unittest.TestCase):
    """Lock the custom PADLOC image contract."""

    def test_dockerfile_patches_the_launcher_at_build_time(self) -> None:
        """Require the PADLOC Dockerfile to bundle one fixed launcher data directory."""
        dockerfile_text = DOCKERFILE.read_text(encoding="utf-8")

        self.assertIn(
            "FROM quay.io/biocontainers/padloc:2.0.0--hdfd78af_1", dockerfile_text
        )
        self.assertIn("ENV PADLOC_DATA_DIR=/opt/padloc-data", dockerfile_text)
        self.assertIn(
            'sed -i "s#mkdir -p \\"\\${SRC_DIR}/../data\\"#mkdir -p \\"${PADLOC_DATA_DIR}\\"#"',
            dockerfile_text,
        )
        self.assertIn(
            'sed -i "s#DATA=\\$(normpath \\"\\${SRC_DIR}/../data\\")#DATA=\\$(normpath \\"${PADLOC_DATA_DIR}\\")#"',
            dockerfile_text,
        )
        self.assertNotIn("--db-update", dockerfile_text)
        self.assertIn("@sha256:", dockerfile_text)
        self.assertIn("ADD --checksum=sha256:", dockerfile_text)
        self.assertIn("7f99b47b75e232b111c18626badb9ac32e8e0b5a", dockerfile_text)
        self.assertIn("ENV LC_ALL=C", dockerfile_text)

    @unittest.skipUnless(
        RUN_DOCKER_TESTS,
        "Set RUN_DOCKER_TESTS=1 to run the PADLOC container image contract test.",
    )
    def test_custom_image_bundles_the_fixed_padloc_database(self) -> None:
        """Require the custom image to bundle one fixed PADLOC database."""
        build_result = run_command(
            [
                "docker",
                "build",
                "--platform",
                "linux/amd64",
                "-f",
                str(DOCKERFILE),
                "-t",
                IMAGE,
                str(DOCKERFILE.parent),
            ]
        )
        self.assertEqual(
            build_result.returncode,
            0,
            msg=f"Docker build failed.\nSTDOUT:\n{build_result.stdout}\nSTDERR:\n{build_result.stderr}",
        )

        result = run_command(
            [
                "docker",
                "run",
                "--rm",
                "--platform",
                "linux/amd64",
                IMAGE,
                "bash",
                "-lc",
                (
                    "set -euo pipefail; "
                    "PADLOC_BIN=$(command -v padloc); "
                    'grep -F \'mkdir -p "/opt/padloc-data"\' "$PADLOC_BIN" >/dev/null; '
                    'grep -F \'DATA=$(normpath "/opt/padloc-data")\' "$PADLOC_BIN" >/dev/null; '
                    "test -f /opt/padloc-data/hmm/padlocdb.hmm; "
                    "padloc --version"
                ),
            ],
            timeout=120,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        output = result.stdout + result.stderr

        self.assertIn("v2.0.0", output)
        self.assertNotIn('DATA=$(normpath "${SRC_DIR}/../data")', output)


if __name__ == "__main__":
    unittest.main()
