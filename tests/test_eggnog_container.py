"""Container contract tests for the custom eggNOG runtime image."""

from __future__ import annotations

import os
import subprocess
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCKERFILE = ROOT / "docker" / "eggnog" / "Dockerfile"
IMAGE = "codex-eggnog:test"
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


class EggnogContainerContractTestCase(unittest.TestCase):
    """Lock the custom eggNOG image contract."""

    def test_dockerfile_uses_locked_v3_with_namespace_export(self) -> None:
        """Require the v3 source contract and qualified namespace export."""
        dockerfile_text = DOCKERFILE.read_text(encoding="utf-8")
        manifest = DOCKERFILE.with_name("pixi.toml").read_text(encoding="utf-8")
        self.assertIn("@sha256:", dockerfile_text)
        self.assertIn("pixi install --locked", dockerfile_text)
        self.assertIn("b3757a6d226047527a729546e58ff530d76a5d7d", manifest)
        self.assertIn("python install_go_export.py", dockerfile_text)
        self.assertIn("python verify_go_export.py", dockerfile_text)
        self.assertNotIn("2.1.13", dockerfile_text)

    @unittest.skipUnless(
        RUN_DOCKER_TESTS,
        "Set RUN_DOCKER_TESTS=1 to run the eggNOG container image contract test.",
    )
    def test_custom_image_qualifies_namespace_export(self) -> None:
        """Build v3 and compare its namespace export with the original engine."""
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
                "python /opt/nf-annotation/eggnog/verify_go_export.py && emapper.py --version",
            ],
            timeout=120,
        )
        self.assertEqual(
            result.returncode,
            0,
            msg=f"Custom eggNOG image contract failed.\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}",
        )
        self.assertIn("emapper-3.0.0-beta6", result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
