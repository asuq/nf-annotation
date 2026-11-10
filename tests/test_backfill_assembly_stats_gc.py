"""Tests for backfilling gc_content into previous-run assembly stats."""

from __future__ import annotations

import csv
import os
import stat
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "bin" / "backfill_assembly_stats_gc.sh"


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read one TSV file into a header and row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames is not None
        return reader.fieldnames, list(reader)


class BackfillAssemblyStatsGcTestCase(unittest.TestCase):
    """Cover the one-off previous-run GC backfill helper."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write one UTF-8 text fixture file."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def write_tsv_rows(
        self,
        path: Path,
        header: tuple[str, ...],
        rows: list[dict[str, str]],
    ) -> Path:
        """Write one TSV fixture with a fixed header."""
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(header), delimiter="\t")
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
        return path

    def install_fake_seqtk(self, directory: Path) -> Path:
        """Install a fake seqtk executable that supports `seqtk comp`."""
        seqtk_path = directory / "seqtk"
        seqtk_path.write_text(
            textwrap.dedent(
                """\
                #!/usr/bin/env python3
                import sys
                from pathlib import Path

                def read_fasta(path: Path):
                    header = None
                    chunks = []
                    for raw_line in path.read_text(encoding="utf-8").splitlines():
                        line = raw_line.strip()
                        if not line:
                            continue
                        if line.startswith(">"):
                            if header is not None:
                                yield header, "".join(chunks)
                            header = line[1:]
                            chunks = []
                            continue
                        if header is None:
                            print("invalid FASTA", file=sys.stderr)
                            return
                        chunks.append(line)
                    if header is not None:
                        yield header, "".join(chunks)

                if len(sys.argv) != 3 or sys.argv[1] != "comp":
                    print("unsupported arguments", file=sys.stderr)
                    sys.exit(1)

                for name, sequence in read_fasta(Path(sys.argv[2])):
                    upper = sequence.upper()
                    print(
                        f"{name}\\t{len(sequence)}\\t"
                        f"{upper.count('A')}\\t{upper.count('C')}\\t"
                        f"{upper.count('G')}\\t{upper.count('T')}"
                    )
                """
            ),
            encoding="utf-8",
        )
        seqtk_path.chmod(seqtk_path.stat().st_mode | stat.S_IXUSR)
        return seqtk_path

    def write_previous_run_fixture(
        self,
        root: Path,
        *,
        assembly_rows: list[dict[str, str]],
    ) -> Path:
        """Create one minimal previous-run results tree."""
        results_dir = root / "results"
        self.write_tsv_rows(
            results_dir / "cohort" / "assembly_stats" / "assembly_stats.tsv",
            ("accession", "n50", "scaffolds", "genome_size"),
            assembly_rows,
        )
        return results_dir

    def write_staged_fasta(
        self, results_dir: Path, accession: str, filename: str, sequence: str
    ) -> Path:
        """Write one staged FASTA fixture."""
        return self.write_text_file(
            results_dir / "samples" / accession / "staged" / filename,
            f">{filename}\n{sequence}\n",
        )

    def run_script(
        self,
        *,
        results_dir: Path,
        path_prefix: Path,
        extra_args: list[str] | None = None,
    ) -> subprocess.CompletedProcess[str]:
        """Run the backfill helper with a controlled PATH."""
        environment = os.environ.copy()
        environment["PATH"] = f"{path_prefix}:{environment['PATH']}"
        return subprocess.run(
            [
                "bash",
                str(SCRIPT),
                "--results-dir",
                str(results_dir),
                *(extra_args or []),
            ],
            check=False,
            capture_output=True,
            text=True,
            env=environment,
        )

    def test_backfills_gc_content_to_default_sibling_output(self) -> None:
        """Append gc_content while preserving row order and prior metrics."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            results_dir = self.write_previous_run_fixture(
                tmpdir,
                assembly_rows=[
                    {
                        "accession": "ACC2",
                        "n50": "99999",
                        "scaffolds": "9",
                        "genome_size": "999999",
                    },
                    {
                        "accession": "ACC1",
                        "n50": "50000",
                        "scaffolds": "2",
                        "genome_size": "800000",
                    },
                ],
            )
            self.write_staged_fasta(results_dir, "ACC2", "ID2.fasta", "GGGGAAAA")
            self.write_text_file(
                results_dir / "samples" / "ACC2" / "staged" / "ID2.fasta.fai", ""
            )
            self.write_staged_fasta(results_dir, "ACC1", "ID1.fasta", "AACCGGTT")

            result = self.run_script(results_dir=results_dir, path_prefix=tmpdir)

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("INFO: Backfill configuration:", result.stderr)
            self.assertIn("INFO: Queued accession ACC2: row=1", result.stderr)
            self.assertIn("INFO: Queued accession ACC1: row=2", result.stderr)
            self.assertIn(
                "INFO: Starting GC backfill workers: rows=2 jobs=1", result.stderr
            )
            self.assertIn(
                "INFO: Backfilled accession ACC2: row=1 gc_content=50", result.stderr
            )
            self.assertIn(
                "INFO: Backfilled accession ACC1: row=2 gc_content=50", result.stderr
            )
            self.assertIn("INFO: Writing merged GC backfill output to", result.stderr)
            self.assertIn("INFO: GC backfill complete: rows=2 output=", result.stderr)
            output_path = (
                results_dir / "cohort" / "assembly_stats" / "assembly_stats.with_gc.tsv"
            )
            header, rows = read_tsv(output_path)
            self.assertEqual(
                header,
                ["accession", "n50", "scaffolds", "genome_size", "gc_content"],
            )
            self.assertEqual([row["accession"] for row in rows], ["ACC2", "ACC1"])
            self.assertEqual(rows[0]["n50"], "99999")
            self.assertEqual(rows[0]["scaffolds"], "9")
            self.assertEqual(rows[0]["genome_size"], "999999")
            self.assertEqual(rows[0]["gc_content"], "50")
            self.assertEqual(rows[1]["gc_content"], "50")

    def test_parallel_jobs_preserve_input_row_order(self) -> None:
        """Keep output rows in the same order when multiple workers run."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            results_dir = self.write_previous_run_fixture(
                tmpdir,
                assembly_rows=[
                    {
                        "accession": "ACC3",
                        "n50": "3",
                        "scaffolds": "3",
                        "genome_size": "30",
                    },
                    {
                        "accession": "ACC1",
                        "n50": "1",
                        "scaffolds": "1",
                        "genome_size": "10",
                    },
                    {
                        "accession": "ACC2",
                        "n50": "2",
                        "scaffolds": "2",
                        "genome_size": "20",
                    },
                ],
            )
            self.write_staged_fasta(results_dir, "ACC3", "ID3.fasta", "GGGG")
            self.write_staged_fasta(results_dir, "ACC1", "ID1.fasta", "AAAA")
            self.write_staged_fasta(results_dir, "ACC2", "ID2.fasta", "CCCC")

            result = self.run_script(
                results_dir=results_dir,
                path_prefix=tmpdir,
                extra_args=["--jobs", "2"],
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn(
                "INFO: Detected 3 assembly-stats row(s) for GC backfill.", result.stderr
            )
            self.assertIn(
                "INFO: Starting GC backfill workers: rows=3 jobs=2", result.stderr
            )
            self.assertIn(
                "INFO: Backfilled accession ACC3: row=1 gc_content=100", result.stderr
            )
            self.assertIn(
                "INFO: Backfilled accession ACC1: row=2 gc_content=0", result.stderr
            )
            self.assertIn(
                "INFO: Backfilled accession ACC2: row=3 gc_content=100", result.stderr
            )
            self.assertIn("INFO: Merging GC backfill results into", result.stderr)
            self.assertIn("INFO: GC backfill complete: rows=3 output=", result.stderr)
            _header, rows = read_tsv(
                results_dir / "cohort" / "assembly_stats" / "assembly_stats.with_gc.tsv"
            )
            self.assertEqual(
                [row["accession"] for row in rows], ["ACC3", "ACC1", "ACC2"]
            )
            self.assertEqual([row["gc_content"] for row in rows], ["100", "0", "100"])

    def test_parallel_jobs_accept_staged_paths_with_spaces(self) -> None:
        """Handle staged FASTA paths containing spaces during parallel execution."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            seqtk_path = self.install_fake_seqtk(tmpdir)
            spaced_root = tmpdir / "with spaces"
            results_dir = self.write_previous_run_fixture(
                spaced_root,
                assembly_rows=[
                    {
                        "accession": "ACC5",
                        "n50": "5",
                        "scaffolds": "1",
                        "genome_size": "8",
                    }
                ],
            )
            self.write_staged_fasta(results_dir, "ACC5", "ID5.fasta", "AACCGGTT")

            result = self.run_script(
                results_dir=results_dir,
                path_prefix=tmpdir,
                extra_args=["--jobs", "2", "--seqtk-binary", str(seqtk_path)],
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("INFO: Backfill configuration:", result.stderr)
            self.assertIn(
                "INFO: Backfilled accession ACC5: row=1 gc_content=50", result.stderr
            )
            _header, rows = read_tsv(
                results_dir / "cohort" / "assembly_stats" / "assembly_stats.with_gc.tsv"
            )
            self.assertEqual(rows[0]["accession"], "ACC5")
            self.assertEqual(rows[0]["gc_content"], "50")

    def test_excludes_ambiguous_bases_from_gc_denominator(self) -> None:
        """Ignore ambiguous bases when computing gc_content."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            results_dir = self.write_previous_run_fixture(
                tmpdir,
                assembly_rows=[
                    {
                        "accession": "ACC1",
                        "n50": "1",
                        "scaffolds": "1",
                        "genome_size": "8",
                    }
                ],
            )
            self.write_staged_fasta(results_dir, "ACC1", "ID1.fasta", "ACGTNNNN")

            result = self.run_script(results_dir=results_dir, path_prefix=tmpdir)

            self.assertEqual(result.returncode, 0, result.stderr)
            _header, rows = read_tsv(
                results_dir / "cohort" / "assembly_stats" / "assembly_stats.with_gc.tsv"
            )
            self.assertEqual(rows[0]["gc_content"], "50")

    def test_fails_when_staged_directory_is_missing(self) -> None:
        """Report a missing staged directory for one accession."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            results_dir = self.write_previous_run_fixture(
                tmpdir,
                assembly_rows=[
                    {
                        "accession": "ACC1",
                        "n50": "1",
                        "scaffolds": "1",
                        "genome_size": "8",
                    }
                ],
            )

            result = self.run_script(results_dir=results_dir, path_prefix=tmpdir)

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("INFO: Backfill configuration:", result.stderr)
            self.assertIn("Missing staged directory for ACC1", result.stderr)

    def test_fails_when_multiple_staged_fastas_are_present(self) -> None:
        """Reject staged directories with more than one FASTA candidate."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            results_dir = self.write_previous_run_fixture(
                tmpdir,
                assembly_rows=[
                    {
                        "accession": "ACC1",
                        "n50": "1",
                        "scaffolds": "1",
                        "genome_size": "8",
                    }
                ],
            )
            self.write_staged_fasta(results_dir, "ACC1", "ID1.fasta", "ACGT")
            self.write_staged_fasta(results_dir, "ACC1", "ID1_copy.fasta", "ACGT")

            result = self.run_script(results_dir=results_dir, path_prefix=tmpdir)

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("INFO: Backfill configuration:", result.stderr)
            self.assertIn(
                "Expected exactly one staged FASTA for ACC1, found multiple",
                result.stderr,
            )

    def test_fails_when_input_lacks_required_column(self) -> None:
        """Reject assembly stats tables without the required metrics."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            results_dir = tmpdir / "results"
            self.write_tsv_rows(
                results_dir / "cohort" / "assembly_stats" / "assembly_stats.tsv",
                ("accession", "n50", "genome_size"),
                [{"accession": "ACC1", "n50": "1", "genome_size": "8"}],
            )

            result = self.run_script(results_dir=results_dir, path_prefix=tmpdir)

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("INFO: Backfill configuration:", result.stderr)
            self.assertIn("missing the scaffolds column", result.stderr)

    def test_fails_when_staged_fasta_has_no_canonical_bases(self) -> None:
        """Reject staged genomes whose canonical base count is zero."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            results_dir = self.write_previous_run_fixture(
                tmpdir,
                assembly_rows=[
                    {
                        "accession": "ACC1",
                        "n50": "1",
                        "scaffolds": "1",
                        "genome_size": "8",
                    }
                ],
            )
            self.write_staged_fasta(results_dir, "ACC1", "ID1.fasta", "NNNNNNNN")

            result = self.run_script(results_dir=results_dir, path_prefix=tmpdir)

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("INFO: Starting accession ACC1: row=1", result.stderr)
            self.assertIn(
                "INFO: Backfill failed for accession ACC1: row=1", result.stderr
            )
            self.assertIn("Failed to calculate gc_content for ACC1", result.stderr)

    def test_parallel_jobs_report_failing_accession(self) -> None:
        """Surface the accession-specific worker error when one parallel job fails."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            results_dir = self.write_previous_run_fixture(
                tmpdir,
                assembly_rows=[
                    {
                        "accession": "ACC1",
                        "n50": "1",
                        "scaffolds": "1",
                        "genome_size": "8",
                    },
                    {
                        "accession": "ACC2",
                        "n50": "2",
                        "scaffolds": "1",
                        "genome_size": "8",
                    },
                ],
            )
            self.write_staged_fasta(results_dir, "ACC1", "ID1.fasta", "AACCGGTT")
            self.write_staged_fasta(results_dir, "ACC2", "ID2.fasta", "NNNNNNNN")

            result = self.run_script(
                results_dir=results_dir,
                path_prefix=tmpdir,
                extra_args=["--jobs", "2"],
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "INFO: Starting GC backfill workers: rows=2 jobs=2", result.stderr
            )
            self.assertIn(
                "INFO: Backfilled accession ACC1: row=1 gc_content=50", result.stderr
            )
            self.assertIn(
                "INFO: Backfill failed for accession ACC2: row=2", result.stderr
            )
            self.assertIn("Failed to calculate gc_content for ACC2", result.stderr)

    def test_rejects_invalid_jobs_argument(self) -> None:
        """Fail early when --jobs is not a positive integer."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            self.install_fake_seqtk(tmpdir)
            results_dir = self.write_previous_run_fixture(
                tmpdir,
                assembly_rows=[
                    {
                        "accession": "ACC1",
                        "n50": "1",
                        "scaffolds": "1",
                        "genome_size": "8",
                    }
                ],
            )
            self.write_staged_fasta(results_dir, "ACC1", "ID1.fasta", "AACCGGTT")

            result = self.run_script(
                results_dir=results_dir,
                path_prefix=tmpdir,
                extra_args=["--jobs", "0"],
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertNotIn("INFO: Backfill configuration:", result.stderr)
            self.assertIn("--jobs must be a positive integer", result.stderr)


if __name__ == "__main__":
    unittest.main()
