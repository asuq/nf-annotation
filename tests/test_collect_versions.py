"""Tests for the final tool and database provenance collector."""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import collect_versions  # noqa: E402


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a TSV file into row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class CollectVersionsTestCase(unittest.TestCase):
    """Cover version-file parsing and runtime/context reporting."""

    def test_inherited_versions_keep_their_origin_and_do_not_replace_current_versions(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory(dir="/tmp") as directory:
            root = Path(directory)
            inherited = root / "inherited"
            inherited.mkdir()
            source = inherited / "report.tsv"
            source.write_text(
                "component\tkind\tversion\timage_or_path\tnotes\n"
                "checkm2\ttool\told-version\tNA\treused source fixture\n"
            )
            current = root / "versions.yml"
            current.write_text('"CHECKM2":\n  checkm2: "new-version"\n')
            output = root / "versions.tsv"
            self.assertEqual(
                collect_versions.main(
                    [
                        "--version-file",
                        str(current),
                        "--inherited-version-dir",
                        str(inherited),
                        "--output",
                        str(output),
                    ]
                ),
                0,
            )
            rows = [
                row
                for row in read_tsv(output)
                if row["component"] == "checkm2" and row["kind"] == "tool"
            ]
            self.assertEqual(
                {row["version"] for row in rows}, {"old-version", "new-version"}
            )
            self.assertEqual(
                next(row["notes"] for row in rows if row["version"] == "old-version"),
                "reused source fixture",
            )

    def test_malformed_or_missing_inherited_provenance_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(dir="/tmp") as directory:
            root = Path(directory)
            inherited = root / "inherited"
            inherited.mkdir()
            command = [
                "--inherited-version-dir",
                str(inherited),
                "--output",
                str(root / "out.tsv"),
            ]
            self.assertEqual(collect_versions.main(command), 1)
            path = inherited / "report.tsv"
            path.write_text("component\tversion\ncheckm2\t1\n")
            self.assertEqual(collect_versions.main(command), 1)
            path.write_text(
                "component\tkind\tversion\timage_or_path\tnotes\ncheckm2\ttool\n"
            )
            self.assertEqual(collect_versions.main(command), 1)

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write test text to a file and return its path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def test_collects_versions_resources_and_containers(self) -> None:
        """Write a stable TSV from process versions plus config context."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            versions_a = self.write_text_file(
                tmpdir / "validate_versions.yml",
                '\n'.join(
                    [
                        '"VALIDATE_INPUTS":',
                        '  python: "3.12.1"',
                        '  script: "bin/validate_inputs.py"',
                    ]
                )
                + "\n",
            )
            versions_b = self.write_text_file(
                tmpdir / "checkm2_versions.yml",
                '\n'.join(
                    [
                        '"CHECKM2":',
                        '  checkm2: "1.0.2"',
                    ]
                )
                + "\n",
            )
            versions_c = self.write_text_file(
                tmpdir / "codetta_versions.yml",
                '\n'.join(
                    [
                        '"CODETTA":',
                        '  codetta: "v2.0"',
                        '  codetta_source_commit: "863359ed326276602d44e48227b6003ac6ffd266"',
                    ]
                )
                + "\n",
            )
            output = tmpdir / "tool_and_db_versions.tsv"

            exit_code = collect_versions.main(
                [
                    "--version-file",
                    str(versions_a),
                    "--version-file",
                    str(versions_b),
                    "--version-file",
                    str(versions_c),
                    "--nextflow-version",
                    "25.04.8",
                    "--pipeline-version",
                    "0.1.0",
                    "--git-commit",
                    "abc1234",
                    "--container-engine",
                    "singularity",
                    "--checkm2-db",
                    "/db/checkm2",
                    "--codetta-db",
                    "/db/codetta",
                    "--codetta-db-label",
                    "Pfam35_Shulgina21",
                    "--taxdump",
                    "/db/taxdump",
                    "--busco-lineage",
                    "bacillota_odb12",
                    "--busco-lineage",
                    "mycoplasmatota_odb12",
                    "--busco-db",
                    "/db/busco",
                    "--eggnog-db",
                    "/db/eggnog",
                    "--container-ref",
                    "python=python:3.12",
                    "--container-ref",
                    "checkm2=quay.io/biocontainers/checkm2:1.0.2",
                    "--container-ref",
                    "codetta=quay.io/asuq1617/codetta:2.0",
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            rows = read_tsv(output)
            self.assertEqual(list(rows[0]), list(collect_versions.OUTPUT_COLUMNS))

            row_map = {
                (row["component"], row["kind"], row["notes"]): row
                for row in rows
            }
            self.assertEqual(
                row_map[("nextflow", "runtime", "workflow")]["version"],
                "25.04.8",
            )
            self.assertEqual(
                row_map[("busco_datasets", "database", "Configured BUSCO lineages in order")][
                    "version"
                ],
                "bacillota_odb12;mycoplasmatota_odb12",
            )
            self.assertEqual(
                row_map[("bacillota_odb12", "database", "BUSCO lineage dataset")][
                    "image_or_path"
                ],
                "/db/busco/bacillota_odb12",
            )
            self.assertEqual(
                row_map[("python", "container", "params container reference")][
                    "image_or_path"
                ],
                "python:3.12",
            )
            self.assertEqual(
                row_map[("codetta_db", "database", "Codetta profile database")]["version"],
                "Pfam35_Shulgina21",
            )
            self.assertEqual(
                row_map[("codetta_db", "database", "Codetta profile database")][
                    "image_or_path"
                ],
                "/db/codetta",
            )
            self.assertEqual(
                row_map[
                    (
                        "codetta_source_commit",
                        "tool",
                        "reported by CODETTA",
                    )
                ]["version"],
                "863359ed326276602d44e48227b6003ac6ffd266",
            )
            self.assertEqual(
                row_map[("codetta", "tool", "reported by CODETTA")]["version"],
                "v2.0",
            )
            self.assertEqual(
                row_map[("codetta", "container", "params container reference")]["image_or_path"],
                "quay.io/asuq1617/codetta:2.0",
            )
            self.assertEqual(
                row_map[("python", "runtime", "reported by VALIDATE_INPUTS")]["version"],
                "3.12.1",
            )
            self.assertEqual(
                row_map[("VALIDATE_INPUTS", "pipeline", "script")]["image_or_path"],
                "bin/validate_inputs.py",
            )
            self.assertEqual(
                row_map[("checkm2", "tool", "reported by CHECKM2")]["version"],
                "1.0.2",
            )
            self.assertNotIn(("use_biocontainers", "pipeline", "config flag"), row_map)

    def test_missing_version_file_fails(self) -> None:
        """Fail cleanly when an input versions file is missing."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            output = tmpdir / "tool_and_db_versions.tsv"

            exit_code = collect_versions.main(
                [
                    "--version-file",
                    str(tmpdir / "missing.yml"),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse(output.exists())

    def test_ignores_trailing_heredoc_delimiters_in_versions_files(self) -> None:
        """Tolerate cached versions files that include a literal heredoc terminator."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            versions = self.write_text_file(
                tmpdir / "versions.yml",
                '\n'.join(
                    [
                        '"DOWNLOAD_BUSCO_DATASET":',
                        '  busco: "BUSCO 6.0.0"',
                        "EOF",
                    ]
                )
                + "\n",
            )
            output = tmpdir / "tool_and_db_versions.tsv"

            exit_code = collect_versions.main(
                [
                    "--version-file",
                    str(versions),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            rows = read_tsv(output)
            busco_rows = [row for row in rows if row["component"] == "busco"]
            self.assertEqual(len(busco_rows), 1)
            self.assertEqual(busco_rows[0]["version"], "BUSCO 6.0.0")

    def test_rejects_multiline_quoted_values_in_versions_files(self) -> None:
        """Fail when a versions file contains a split quoted value."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            versions = self.write_text_file(
                tmpdir / "versions.yml",
                '\n'.join(
                    [
                        '"STAGE_INPUTS":',
                        '  seqtk: "',
                        'NA"',
                        '  samtools: "NA"',
                    ]
                )
                + "\n",
            )
            output = tmpdir / "tool_and_db_versions.tsv"

            exit_code = collect_versions.main(
                [
                    "--version-file",
                    str(versions),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse(output.exists())

    def test_rejects_indented_header_lines_in_versions_files(self) -> None:
        """Fail when a header line is indented like a malformed entry."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            versions = self.write_text_file(
                tmpdir / "versions.yml",
                '\n'.join(
                    [
                        '    "PER_SAMPLE_ANNOTATION:EGGNOG":',
                        '      eggnog_mapper: "2.1.13"',
                    ]
                )
                + "\n",
            )
            output = tmpdir / "tool_and_db_versions.tsv"

            exit_code = collect_versions.main(
                [
                    "--version-file",
                    str(versions),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse(output.exists())

    def test_collects_single_line_eggnog_version_entry(self) -> None:
        """Keep a single-line eggNOG version entry intact."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            versions = self.write_text_file(
                tmpdir / "versions.yml",
                '\n'.join(
                    [
                        '"PER_SAMPLE_ANNOTATION:EGGNOG":',
                        '  eggnog_mapper: "2.1.13"',
                    ]
                )
                + "\n",
            )
            output = tmpdir / "tool_and_db_versions.tsv"

            exit_code = collect_versions.main(
                [
                    "--version-file",
                    str(versions),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            rows = read_tsv(output)
            eggnog_rows = [
                row for row in rows if row["component"] == "eggnog_mapper"
            ]
            self.assertEqual(len(eggnog_rows), 1)
            self.assertEqual(eggnog_rows[0]["version"], "2.1.13")

    def test_collects_canonical_process_rows_from_version_dir(self) -> None:
        """Keep one canonical row set when a process reports identical versions repeatedly."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            version_dir = tmpdir / "version_files"
            self.write_text_file(
                version_dir / "versions01.yml",
                '\n'.join(
                    [
                        '"PER_SAMPLE_QC:BARRNAP":',
                        '  barrnap: "0.9"',
                    ]
                )
                + "\n",
            )
            self.write_text_file(
                version_dir / "versions02.yml",
                '\n'.join(
                    [
                        '"PER_SAMPLE_QC:BARRNAP":',
                        '  barrnap: "0.9"',
                    ]
                )
                + "\n",
            )
            output = tmpdir / "tool_and_db_versions.tsv"

            exit_code = collect_versions.main(
                [
                    "--version-dir",
                    str(version_dir),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            rows = read_tsv(output)
            barrnap_rows = [
                row
                for row in rows
                if row["component"] == "barrnap"
                and row["notes"] == "reported by PER_SAMPLE_QC:BARRNAP"
            ]
            self.assertEqual(len(barrnap_rows), 1)

    def test_version_dir_conflicts_fail_for_one_process(self) -> None:
        """Fail when the same process reports conflicting version rows."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            version_dir = tmpdir / "version_files"
            self.write_text_file(
                version_dir / "versions01.yml",
                '\n'.join(
                    [
                        '"PER_SAMPLE_QC:BARRNAP":',
                        '  barrnap: "0.9"',
                    ]
                )
                + "\n",
            )
            self.write_text_file(
                version_dir / "versions02.yml",
                '\n'.join(
                    [
                        '"PER_SAMPLE_QC:BARRNAP":',
                        '  barrnap: "0.10"',
                    ]
                )
                + "\n",
            )
            output = tmpdir / "tool_and_db_versions.tsv"

            exit_code = collect_versions.main(
                [
                    "--version-dir",
                    str(version_dir),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse(output.exists())

    def test_version_dir_malformed_file_still_fails(self) -> None:
        """Keep malformed staged versions files as hard failures."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            version_dir = tmpdir / "version_files"
            self.write_text_file(
                version_dir / "versions01.yml",
                '\n'.join(
                    [
                        '    "COHORT_16S:BUILD_COHORT_16S":',
                        '      python: "3.12.13"',
                        '    EOF',
                    ]
                )
                + "\n",
            )
            output = tmpdir / "tool_and_db_versions.tsv"

            exit_code = collect_versions.main(
                [
                    "--version-dir",
                    str(version_dir),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse(output.exists())

    def test_supports_mixed_version_file_and_version_dir_inputs(self) -> None:
        """Allow one-off version files to be combined with a staged versions directory."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            version_dir = tmpdir / "version_files"
            direct_file = self.write_text_file(
                tmpdir / "direct_versions.yml",
                '\n'.join(
                    [
                        '"INPUT_VALIDATION_AND_STAGING:VALIDATE_INPUTS":',
                        '  python: "3.12.13"',
                    ]
                )
                + "\n",
            )
            self.write_text_file(
                version_dir / "versions01.yml",
                '\n'.join(
                    [
                        '"PER_SAMPLE_QC:BARRNAP":',
                        '  barrnap: "0.9"',
                    ]
                )
                + "\n",
            )
            output = tmpdir / "tool_and_db_versions.tsv"

            exit_code = collect_versions.main(
                [
                    "--version-file",
                    str(direct_file),
                    "--version-dir",
                    str(version_dir),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            rows = read_tsv(output)
            components = {row["component"] for row in rows}
            self.assertIn("python", components)
            self.assertIn("barrnap", components)

    def test_invalid_container_ref_fails(self) -> None:
        """Fail cleanly when a container reference is not name=value."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            output = tmpdir / "tool_and_db_versions.tsv"

            exit_code = collect_versions.main(
                [
                    "--container-ref",
                    "python",
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse(output.exists())

    def test_malformed_versions_file_fails(self) -> None:
        """Fail cleanly when a versions file is not simple YAML-like key/value text."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            versions = self.write_text_file(
                tmpdir / "broken_versions.yml",
                '"VALIDATE_INPUTS"\n  python "3.12.1"\n',
            )
            output = tmpdir / "tool_and_db_versions.tsv"

            exit_code = collect_versions.main(
                [
                    "--version-file",
                    str(versions),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
