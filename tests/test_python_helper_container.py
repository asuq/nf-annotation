"""Container contract tests for the shared Python helper image."""

from __future__ import annotations

import os
import subprocess
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCKERFILE = ROOT / "docker" / "python_helper" / "Dockerfile"
IMAGE = "codex-python-helper:test"
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


class PythonHelperContainerContractTestCase(unittest.TestCase):
    """Lock the shared Python helper image contract."""

    def test_dockerfile_installs_biopython_numpy_and_scipy(self) -> None:
        """Require the helper Dockerfile to install the Python reporting stack."""
        dockerfile_text = DOCKERFILE.read_text(encoding="utf-8")

        self.assertIn("FROM python:3.12", dockerfile_text)
        self.assertIn("biopython==1.84", dockerfile_text)
        self.assertIn("numpy==2.1.1", dockerfile_text)
        self.assertIn("scipy==1.14.1", dockerfile_text)

    @unittest.skipUnless(
        RUN_DOCKER_TESTS,
        "Set RUN_DOCKER_TESTS=1 to run the shared Python helper image contract test.",
    )
    def test_custom_image_imports_biopython_numpy_and_scipy_on_linux_amd64(
        self,
    ) -> None:
        """Build the helper image and require imports to work on amd64."""
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
                "python",
                "-c",
                (
                    "import Bio, numpy, scipy; "
                    "print(Bio.__version__); "
                    "print(numpy.__version__); "
                    "print(scipy.__version__)"
                ),
            ],
            timeout=30,
        )
        self.assertEqual(
            result.returncode,
            0,
            msg=f"Python helper image contract failed.\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}",
        )
        self.assertIn("1.84", result.stdout)
        self.assertIn("2.1.1", result.stdout)
        self.assertIn("1.14.1", result.stdout)


if __name__ == "__main__":
    unittest.main()
