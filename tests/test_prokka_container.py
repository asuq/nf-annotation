"""Container contract tests for the custom Prokka runtime image."""

from __future__ import annotations

import os
import subprocess
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCKERFILE = ROOT / "docker" / "prokka" / "Dockerfile"
WRAPPER = ROOT / "docker" / "prokka" / "minced-wrapper.sh"
IMAGE = "codex-prokka:test"
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


class ProkkaContainerContractTestCase(unittest.TestCase):
    """Lock the custom Prokka image contract."""

    def test_dockerfile_wraps_minced_at_build_time(self) -> None:
        """Require the Dockerfile to install a deterministic minced wrapper."""
        dockerfile_text = DOCKERFILE.read_text(encoding="utf-8")
        wrapper_text = WRAPPER.read_text(encoding="utf-8")

        self.assertIn(
            "FROM quay.io/biocontainers/prokka:1.15.6--pl5321hdfd78af_0",
            dockerfile_text,
        )
        self.assertIn('MINCED_BIN="$(command -v minced)"', dockerfile_text)
        self.assertIn('mv "${MINCED_BIN}" "${MINCED_BIN}.real"', dockerfile_text)
        self.assertIn(
            "COPY docker/prokka/minced-wrapper.sh /usr/local/bin/minced",
            dockerfile_text,
        )
        self.assertIn('test "${first_line}" = "minced 0.4.2"', dockerfile_text)
        self.assertIn("-XX:+PerfDisableSharedMem", wrapper_text)
        self.assertIn(
            'real_path="$(readlink -f "${script_dir}/minced.real")"', wrapper_text
        )
        self.assertIn('jar_dir="$(dirname "${real_path}")"', wrapper_text)
        self.assertIn(
            'exec "${java_bin}" -XX:+PerfDisableSharedMem -jar "${jar_dir}/minced.jar" "$@"',
            wrapper_text,
        )

    @unittest.skipUnless(
        RUN_DOCKER_TESTS,
        "Set RUN_DOCKER_TESTS=1 to run the Prokka container image contract test.",
    )
    def test_custom_image_keeps_minced_version_probe_clean(self) -> None:
        """Build the custom image and require a clean minced version probe."""
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
                ".",
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
                    "first_line=\"$(minced --version 2>&1 | sed -n '1p')\"; "
                    'test "$first_line" = "minced 0.4.2"; '
                    'case "$first_line" in *hsperfdata*) exit 1 ;; esac; '
                    "test -L /usr/local/bin/minced.real; "
                    "prokka --version"
                ),
            ],
            timeout=30,
        )
        self.assertEqual(
            result.returncode,
            0,
            msg=f"Custom Prokka image contract failed.\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}",
        )
        self.assertIn("prokka 1.15.6", result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
