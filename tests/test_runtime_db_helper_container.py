"""Container contract tests for the runtime database helper image."""

from __future__ import annotations

import os
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCKERFILE = ROOT / "docker" / "runtime_db_helper" / "Dockerfile"
IMAGE_TAG = "codex-runtime-db-helper:test"
RUN_DOCKER_TESTS = os.environ.get("RUN_DOCKER_TESTS") == "1"


def run_command(args: list[str]) -> subprocess.CompletedProcess[str]:
    """Run one subprocess and return the completed process."""
    return subprocess.run(
        args,
        check=False,
        capture_output=True,
        text=True,
        cwd=ROOT,
    )


class RuntimeDbHelperContainerContractTestCase(unittest.TestCase):
    """Lock the dedicated runtime database helper image contract."""

    def test_dockerfile_installs_aria2_and_helper_scripts(self) -> None:
        """Require the helper image to bundle aria2 and the repo CLIs."""
        dockerfile_text = DOCKERFILE.read_text(encoding="utf-8")

        self.assertIn("FROM python:3.12-slim", dockerfile_text)
        self.assertIn("aria2", dockerfile_text)
        self.assertIn("procps", dockerfile_text)
        self.assertIn("COPY bin/prepare_runtime_databases.py", dockerfile_text)
        self.assertIn("COPY bin/finalise_runtime_database.py", dockerfile_text)
        self.assertIn("COPY bin/merge_runtime_database_reports.py", dockerfile_text)
        self.assertIn(
            "COPY assets/runtime/runtime_database_sources.json", dockerfile_text
        )
        self.assertIn(
            "ENV NF_MYCO_RUNTIME_DB_SOURCE_MANIFEST=/opt/nf-myco_annotation/runtime_database_sources.json",
            dockerfile_text,
        )
        helper_text = (ROOT / "bin" / "prepare_runtime_databases.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("NF_MYCO_RUNTIME_DB_SOURCE_MANIFEST", helper_text)
        self.assertIn("--codetta-source", helper_text)
        self.assertIn("--codetta-dest", helper_text)
        self.assertIn("--codetta-version", helper_text)

    @unittest.skipUnless(
        RUN_DOCKER_TESTS,
        "Set RUN_DOCKER_TESTS=1 to build the runtime database helper image.",
    )
    def test_built_helper_image_exposes_aria2_and_the_manifest(self) -> None:
        """Build the helper image and require the bundled tools to be present."""
        build_result = run_command(
            [
                "docker",
                "build",
                "--platform",
                "linux/amd64",
                "-f",
                str(DOCKERFILE),
                "-t",
                IMAGE_TAG,
                ".",
            ]
        )
        self.assertEqual(
            build_result.returncode,
            0,
            msg=f"Docker build failed.\nSTDOUT:\n{build_result.stdout}\nSTDERR:\n{build_result.stderr}",
        )

        run_result = run_command(
            [
                "docker",
                "run",
                "--rm",
                "--platform",
                "linux/amd64",
                IMAGE_TAG,
                "bash",
                "-lc",
                (
                    "set -euo pipefail; "
                    "command -v ps >/dev/null; "
                    "command -v aria2c >/dev/null; "
                    "command -v prepare_runtime_databases.py >/dev/null; "
                    "command -v finalise_runtime_database.py >/dev/null; "
                    "command -v merge_runtime_database_reports.py >/dev/null; "
                    "prepare_runtime_databases.py --help 2>&1 | grep -F -- '--codetta-source' >/dev/null; "
                    "prepare_runtime_databases.py --help 2>&1 | grep -F -- '--codetta-dest' >/dev/null; "
                    "prepare_runtime_databases.py --help 2>&1 | grep -F -- '--codetta-version' >/dev/null; "
                    "test -f /opt/nf-myco_annotation/runtime_database_sources.json; "
                    "python3 -c \"import json; data = json.load(open('/opt/nf-myco_annotation/runtime_database_sources.json', encoding='utf-8')); assert 'codetta' in data\""
                ),
            ]
        )
        self.assertEqual(
            run_result.returncode,
            0,
            msg=f"Runtime helper contract failed.\nSTDOUT:\n{run_result.stdout}\nSTDERR:\n{run_result.stderr}",
        )


if __name__ == "__main__":
    unittest.main()
