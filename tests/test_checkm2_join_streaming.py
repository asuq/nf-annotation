"""Behaviour tests for keyed CheckM2 report joins."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path

NEXTFLOW = shutil.which("nextflow") or ""


class Checkm2JoinStreamingTestCase(unittest.TestCase):
    """Verify that keyed joins stream paired CheckM2 reports per accession."""

    def write_workflow(self, path: Path) -> Path:
        """Write one temporary Nextflow workflow that mirrors the QC join pattern."""
        path.write_text(
            textwrap.dedent(
                """\
                nextflow.enable.dsl=2

                workflow {
                    long startMillis = System.currentTimeMillis()

                    checkm2Gcode4 = Channel.of(
                        [accession: 'ACC_EARLY', delay_ms: 0L],
                        [accession: 'ACC_LATE', delay_ms: 300L],
                    ).map { row ->
                        Thread.sleep(row.delay_ms as long)
                        tuple(row.accession, [accession: row.accession], "${row.accession}_gcode4.tsv")
                    }

                    checkm2Gcode11 = Channel.of(
                        [accession: 'ACC_EARLY', delay_ms: 100L],
                        [accession: 'ACC_LATE', delay_ms: 1600L],
                    ).map { row ->
                        Thread.sleep(row.delay_ms as long)
                        tuple(row.accession, "${row.accession}_gcode11.tsv")
                    }

                    checkm2Gcode4
                        .join(checkm2Gcode11)
                        .map { accession, meta, gcode4Report, gcode11Report ->
                            "${accession}\\t${System.currentTimeMillis() - startMillis}"
                        }
                        .view { line -> "EMIT\\t${line}" }
                }
                """
            ),
            encoding="ascii",
        )
        return path

    def run_workflow(
        self, workflow_path: Path, run_dir: Path
    ) -> subprocess.CompletedProcess[str]:
        """Run one temporary Nextflow workflow with local state isolated to the tempdir."""
        environment = os.environ.copy()
        environment["NXF_ANSI_LOG"] = "false"
        environment["NXF_DISABLE_CHECK_LATEST"] = "true"
        environment["NXF_HOME"] = str(run_dir / ".nxf")
        return subprocess.run(
            [NEXTFLOW, "run", str(workflow_path)],
            cwd=run_dir,
            text=True,
            capture_output=True,
            check=False,
            env=environment,
        )

    def parse_emissions(self, stdout: str) -> list[tuple[str, int]]:
        """Parse timestamped emission lines from the workflow stdout."""
        emissions: list[tuple[str, int]] = []
        for line in stdout.splitlines():
            if not line.startswith("EMIT\t"):
                continue
            _, accession, elapsed_millis = line.split("\t")
            emissions.append((accession, int(elapsed_millis)))
        return emissions

    @unittest.skipUnless(
        bool(NEXTFLOW), "nextflow is required for join-streaming tests."
    )
    def test_join_emits_early_pair_before_late_pair_finishes(self) -> None:
        """Emit the first accession well before the slow second pair completes."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            workflow_path = self.write_workflow(tmpdir / "main.nf")

            result = self.run_workflow(workflow_path, tmpdir)

            self.assertEqual(result.returncode, 0, result.stderr)
            emissions = self.parse_emissions(result.stdout)
            self.assertEqual(
                [accession for accession, _ in emissions],
                ["ACC_EARLY", "ACC_LATE"],
                result.stdout,
            )
            self.assertGreaterEqual(len(emissions), 2)
            self.assertGreater(emissions[1][1] - emissions[0][1], 900, emissions)


if __name__ == "__main__":
    unittest.main()
