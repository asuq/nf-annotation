"""Container contract tests for the custom Codetta runtime image."""

from __future__ import annotations

import os
import subprocess
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCKERFILE = ROOT / "docker" / "codetta" / "Dockerfile"
IMAGE = "quay.io/asuq1617/codetta:2.0"
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


class CodettaContainerContractTestCase(unittest.TestCase):
    """Lock the custom Codetta image contract."""

    def test_dockerfile_pins_the_upstream_commit_and_runtime_layout(self) -> None:
        """Require the Dockerfile to pin Codetta and expose its resources."""
        dockerfile_text = DOCKERFILE.read_text(encoding="utf-8")

        self.assertIn("FROM python:3.12-slim", dockerfile_text)
        self.assertIn(
            "ARG CODETTA_REF=863359ed326276602d44e48227b6003ac6ffd266", dockerfile_text
        )
        self.assertIn("ARG CODETTA_VERSION=v2.0", dockerfile_text)
        self.assertIn("ARG CODETTA_SOURCE_BRANCH=main", dockerfile_text)
        self.assertIn("procps", dockerfile_text)
        self.assertIn(
            'python3 -m pip install --no-cache-dir numpy scipy "setuptools<82"',
            dockerfile_text,
        )
        self.assertIn('python3 -c "import pkg_resources"', dockerfile_text)
        self.assertIn("bash setup.sh", dockerfile_text)
        self.assertIn(
            "ln -s /opt/codetta/codetta_align.py /usr/local/bin/codetta_align.py",
            dockerfile_text,
        )
        self.assertIn(
            "ln -s /opt/codetta/codetta_summary.py /usr/local/bin/codetta_summary.py",
            dockerfile_text,
        )
        self.assertIn(
            "ln -s /opt/codetta/codetta_infer.py /usr/local/bin/codetta_infer.py",
            dockerfile_text,
        )
        self.assertIn(
            'LABEL org.opencontainers.image.version="${CODETTA_VERSION}"',
            dockerfile_text,
        )
        self.assertIn("ENV CODETTA_HOME=/opt/codetta", dockerfile_text)

    @unittest.skipUnless(
        RUN_DOCKER_TESTS,
        "Set RUN_DOCKER_TESTS=1 to run the Codetta container image contract test.",
    )
    def test_custom_image_exposes_codetta_scripts_and_resources(self) -> None:
        """Require the built image to expose the Codetta scripts and resources."""
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
                    "command -v ps >/dev/null; "
                    "command -v codetta_align.py >/dev/null; "
                    "command -v codetta_summary.py >/dev/null; "
                    "command -v codetta_infer.py >/dev/null; "
                    'python3 -c "import pkg_resources" >/dev/null; '
                    "test -d /opt/codetta/resources; "
                    "test -f /opt/codetta/.codetta_version; "
                    "test -f /opt/codetta/.codetta_source_commit; "
                    "test -f /opt/codetta/.codetta_source_branch; "
                    "python3 - <<'PY'\n"
                    "from pathlib import Path\n"
                    "meta_root = Path('/opt/codetta')\n"
                    "root = Path('/opt/codetta/resources')\n"
                    "required = [\n"
                    "    'bad_pfams.txt',\n"
                    "    'mito_pfams.txt',\n"
                    "    'pyrrolysine_pfams.txt',\n"
                    "    'selenocysteine_pfams.txt',\n"
                    "    'transposon_pfams.txt',\n"
                    "    'viral_pfams.txt',\n"
                    "]\n"
                    "missing = [name for name in required if not (root / name).is_file()]\n"
                    "if missing:\n"
                    "    raise SystemExit('missing:' + ','.join(missing))\n"
                    "if meta_root.joinpath('.codetta_version').read_text().strip() != 'v2.0':\n"
                    "    raise SystemExit('bad_version')\n"
                    "if meta_root.joinpath('.codetta_source_commit').read_text().strip() != '863359ed326276602d44e48227b6003ac6ffd266':\n"
                    "    raise SystemExit('bad_commit')\n"
                    "if meta_root.joinpath('.codetta_source_branch').read_text().strip() != 'main':\n"
                    "    raise SystemExit('bad_branch')\n"
                    "PY"
                ),
            ],
            timeout=120,
        )
        self.assertEqual(
            result.returncode,
            0,
            msg=f"Codetta container contract failed.\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}",
        )


if __name__ == "__main__":
    unittest.main()
