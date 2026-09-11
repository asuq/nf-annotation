"""Tests for the runtime database preparation helper."""

from __future__ import annotations

import csv
import gzip
import io
import json
import logging
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest
import zipfile
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import prepare_runtime_databases  # noqa: E402


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read one TSV file into row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class PrepareRuntimeDatabasesTestCase(unittest.TestCase):
    """Cover in-place runtime database preparation and reuse."""

    def reset_logging(self) -> None:
        """Reset root logging handlers so stderr capture stays predictable."""
        root_logger = logging.getLogger()
        for handler in list(root_logger.handlers):
            root_logger.removeHandler(handler)

    def run_cli(self, args: list[str]) -> tuple[int, str, str]:
        """Run the helper CLI and capture stdout plus stderr."""
        self.reset_logging()
        stdout = io.StringIO()
        stderr = io.StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            exit_code = prepare_runtime_databases.main(args)
        return exit_code, stdout.getvalue(), stderr.getvalue()

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write one UTF-8 test file and return its path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def read_marker(self, path: Path) -> dict[str, object]:
        """Read one ready marker JSON object."""
        return json.loads(
            (path / prepare_runtime_databases.MARKER_FILE_NAME).read_text(
                encoding="utf-8"
            )
        )

    def create_taxdump_dir(self, path: Path) -> Path:
        """Create one valid taxdump source directory."""
        self.write_text_file(path / "names.dmp", "1\t|\troot\t|\t\t|\tscientific name\t|\n")
        self.write_text_file(path / "nodes.dmp", "1\t|\t1\t|\tno rank\t|\n")
        return path

    def create_checkm2_dir(self, path: Path) -> Path:
        """Create one valid CheckM2 source directory."""
        self.write_text_file(path / "CheckM2_database.dmnd", "fake-dmnd\n")
        return path

    def create_busco_lineage_dir(self, path: Path) -> Path:
        """Create one valid BUSCO lineage source directory."""
        self.write_text_file(path / "dataset.cfg", "name=fake\n")
        return path

    def create_eggnog_dir(self, path: Path) -> Path:
        """Create a synthetic, already prepared v0.4 eggNOG resource."""
        import annotation_resources

        for name in ("eggnog.db", "eggnog_proteins.dmnd", "eggnog.db.taxids.bin", "eggnog.db.fieldpresence.bin", "eggnog.taxa.db", "eggnog.taxa.db.traverse.pkl", "go-basic.obo"):
            self.write_text_file(path / name, "synthetic-prepared-resource\n")
        contract = {"component": "eggnog", "version": "eggNOG7-emapper3.0", "files": annotation_resources.file_records(path), "settings": {"synthetic_fixture": True}}
        payload = {"schema_version": 1, "component": "eggnog", "resource_id": annotation_resources.identity(contract), "contract": contract, "acquisition": {"files": []}}
        self.write_text_file(path / annotation_resources.RESOURCE_FILE, json.dumps(payload))
        return path

    def create_codetta_dir(self, path: Path) -> Path:
        """Create one valid Codetta source directory."""
        self.write_text_file(path / "Pfam-A_enone.hmm", "profile-placeholder\n")
        self.write_text_file(path / "Pfam-A_enone.hmm.h3f", "index-placeholder\n")
        self.write_text_file(path / "Pfam-A_enone.hmm.h3i", "index-placeholder\n")
        self.write_text_file(path / "Pfam-A_enone.hmm.h3m", "index-placeholder\n")
        self.write_text_file(path / "Pfam-A_enone.hmm.h3p", "index-placeholder\n")
        return path

    def create_padloc_dir(self, path: Path) -> Path:
        """Create one valid PADLOC source directory."""
        self.write_text_file(path / "hmm" / "padlocdb.hmm", "hmm-placeholder\n")
        return path

    def create_tar_archive(self, source: Path, archive_path: Path, root_name: str) -> Path:
        """Create one tar archive containing the source directory under one root."""
        archive_path.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(archive_path, "w:gz") as handle:
            handle.add(source, arcname=root_name)
        return archive_path

    def create_zip_archive(self, source: Path, archive_path: Path, root_name: str) -> Path:
        """Create one zip archive containing the source directory under one root."""
        archive_path.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(archive_path, "w") as handle:
            for child in source.rglob("*"):
                if child.is_dir():
                    continue
                relative = child.relative_to(source)
                handle.write(child, Path(root_name) / relative)
        return archive_path

    def create_gzip_file(self, path: Path, content: str) -> Path:
        """Create one gzip-compressed text file."""
        path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(content)
        return path

    def checksum_md5(self, path: Path) -> str:
        """Return the MD5 checksum for one file."""
        return prepare_runtime_databases.compute_checksum(path, "md5")

    def write_remote_manifest(self, path: Path, payload: dict[str, object]) -> Path:
        """Write one remote source manifest JSON file."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        return path

    def build_remote_manifest(
        self,
        *,
        taxdump_url: str,
        checkm2_url: str,
        busco_template: str,
        codetta_url: str | None = None,
        padloc_url: str,
        taxdump_checksum_url: str | None = None,
        taxdump_checksum_value: str | None = None,
        checkm2_checksum_value: str | None = None,
        taxdump_versions: dict[str, dict[str, object]] | None = None,
    ) -> dict[str, object]:
        """Build one minimal curated remote manifest for testing."""
        return {
            "taxdump": {
                "default_version": "current",
                "versions": taxdump_versions
                or {
                    "current": {
                        "kind": "archive",
                        "url": taxdump_url,
                        "archive_name": "taxdump.zip",
                        "checksum": (
                            {
                                "type": "md5",
                                "url": taxdump_checksum_url,
                                "expected_name": "taxdump.zip",
                            }
                            if taxdump_checksum_url
                            else {
                                "type": "md5",
                                "value": taxdump_checksum_value,
                            }
                        ),
                    }
                },
            },
            "checkm2": {
                "default_version": "current",
                "versions": {
                    "current": {
                        "kind": "archive",
                        "url": checkm2_url,
                        "archive_name": "checkm2.tar.gz",
                        "checksum": {
                            "type": "md5",
                            "value": checkm2_checksum_value,
                        },
                    }
                },
            },
            "codetta": {
                "default_version": "current",
                "versions": {
                    "current": {
                        "kind": "archive",
                        "url": codetta_url or "https://example.invalid/codetta.tar.gz",
                        "archive_name": "Pfam-A_enone.tar.gz",
                    }
                },
            },
            "busco": {
                "default_version": "current",
                "versions": {
                    "current": {
                        "kind": "lineage_archives",
                        "lineage_url_template": busco_template,
                        "archive_name_template": "{lineage}.tar.gz",
                    }
                },
            },
            "padloc": {
                "default_version": "current",
                "versions": {
                    "current": {
                        "kind": "archive",
                        "url": padloc_url,
                        "archive_name": "padloc.zip",
                    }
                },
            },
        }

    def make_fake_aria2(
        self,
        source_map: dict[str, Path],
        recorded_calls: list[list[str]],
        failures: dict[str, str] | None = None,
    ) -> mock.Mock:
        """Build one fake aria2 runner that copies local fixtures."""
        failure_map = failures or {}

        def fake_run(
            args: list[str],
            *,
            check: bool = False,
            capture_output: bool = True,
            text: bool = True,
        ) -> subprocess.CompletedProcess[str]:
            self.assertEqual(args[0], "aria2c")
            self.assertFalse(check)
            self.assertTrue(capture_output)
            self.assertTrue(text)
            recorded_calls.append(args)
            url = args[-1]
            if url in failure_map:
                return subprocess.CompletedProcess(args, 1, "", failure_map[url])
            source = source_map[url]
            destination_dir = Path(args[args.index("--dir") + 1])
            destination_name = args[args.index("--out") + 1]
            destination_dir.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination_dir / destination_name)
            return subprocess.CompletedProcess(args, 0, "", "")

        return mock.Mock(side_effect=fake_run)

    def test_main_prepares_all_directory_sources_and_writes_report(self) -> None:
        """Prepare all runtime database classes from directory sources."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            taxdump_source = self.create_taxdump_dir(tmpdir / "sources" / "taxdump")
            checkm2_source = self.create_checkm2_dir(tmpdir / "sources" / "checkm2")
            bacillota_source = self.create_busco_lineage_dir(
                tmpdir / "sources" / "busco" / "bacillota_odb12"
            )
            myco_source = self.create_busco_lineage_dir(
                tmpdir / "sources" / "busco" / "mycoplasmatota_odb12"
            )
            codetta_source = self.create_codetta_dir(tmpdir / "sources" / "codetta")
            eggnog_source = self.create_eggnog_dir(tmpdir / "sources" / "eggnog")
            padloc_source = self.create_padloc_dir(tmpdir / "sources" / "padloc")
            report_path = tmpdir / "report.tsv"

            exit_code, stdout, _stderr = self.run_cli(
                [
                    "--taxdump-source",
                    str(taxdump_source),
                    "--taxdump-dest",
                    str(tmpdir / "prepared" / "taxdump"),
                    "--checkm2-source",
                    str(checkm2_source),
                    "--checkm2-dest",
                    str(tmpdir / "prepared" / "checkm2"),
                    "--busco-lineage-source",
                    f"bacillota_odb12={bacillota_source}",
                    "--busco-lineage-source",
                    f"mycoplasmatota_odb12={myco_source}",
                    "--busco-dest-root",
                    str(tmpdir / "prepared" / "busco"),
                    "--codetta-source",
                    str(codetta_source),
                    "--codetta-dest",
                    str(tmpdir / "prepared" / "codetta"),
                    "--eggnog-source",
                    str(eggnog_source),
                    "--eggnog-dest",
                    str(tmpdir / "prepared" / "eggnog"),
                    "--padloc-source",
                    str(padloc_source),
                    "--padloc-dest",
                    str(tmpdir / "prepared" / "padloc"),
                    "--report",
                    str(report_path),
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertIn("--taxdump", stdout)
            self.assertIn("--checkm2_db", stdout)
            self.assertIn("--codetta_db", stdout)
            self.assertIn("--busco_db", stdout)
            self.assertIn("--eggnog_db", stdout)
            self.assertNotIn("--padloc_db", stdout)

            rows = read_tsv(report_path)
            self.assertEqual(len(rows), 8)
            row_map = {row["component"]: row for row in rows}
            self.assertEqual(row_map["taxdump"]["status"], "prepared")
            self.assertEqual(row_map["checkm2"]["status"], "prepared")
            self.assertEqual(row_map["codetta"]["status"], "prepared")
            self.assertEqual(row_map["busco:bacillota_odb12"]["status"], "prepared")
            self.assertEqual(row_map["busco:mycoplasmatota_odb12"]["status"], "prepared")
            self.assertEqual(
                row_map["busco_root"]["destination"],
                str((tmpdir / "prepared" / "busco").resolve()),
            )

            taxdump_dest = tmpdir / "prepared" / "taxdump"
            checkm2_dest = tmpdir / "prepared" / "checkm2"
            bacillota_dest = tmpdir / "prepared" / "busco" / "bacillota_odb12"
            myco_dest = tmpdir / "prepared" / "busco" / "mycoplasmatota_odb12"
            codetta_dest = tmpdir / "prepared" / "codetta"
            eggnog_dest = tmpdir / "prepared" / "eggnog"
            padloc_dest = tmpdir / "prepared" / "padloc"

            self.assertTrue((taxdump_dest / "names.dmp").is_file())
            self.assertTrue((checkm2_dest / "CheckM2_database.dmnd").is_file())
            self.assertTrue((bacillota_dest / "dataset.cfg").is_file())
            self.assertTrue((myco_dest / "dataset.cfg").is_file())
            self.assertTrue((codetta_dest / "Pfam-A_enone.hmm").is_file())
            self.assertTrue((codetta_dest / "Pfam-A_enone.hmm.h3f").is_file())
            self.assertTrue((eggnog_dest / "eggnog.db").is_file())
            self.assertTrue((padloc_dest / "hmm" / "padlocdb.hmm").is_file())

            self.assertEqual(self.read_marker(taxdump_dest)["component"], "taxdump")
            self.assertEqual(self.read_marker(checkm2_dest)["component"], "checkm2")
            self.assertEqual(
                self.read_marker(bacillota_dest)["component"],
                "busco:bacillota_odb12",
            )
            self.assertEqual(self.read_marker(codetta_dest)["component"], "codetta")
            self.assertEqual(self.read_marker(eggnog_dest)["component"], "eggnog")
            self.assertEqual(self.read_marker(padloc_dest)["component"], "padloc")

    def test_main_accepts_checkm2_raw_dmnd_source(self) -> None:
        """Materialise one raw CheckM2 dmnd file into a prepared directory."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            dmnd_source = self.write_text_file(tmpdir / "sources" / "checkm2.dmnd", "fake\n")
            destination = tmpdir / "prepared" / "checkm2"

            exit_code, _stdout, _stderr = self.run_cli(
                [
                    "--checkm2-source",
                    str(dmnd_source),
                    "--checkm2-dest",
                    str(destination),
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertTrue((destination / "checkm2.dmnd").is_file())
            self.assertEqual(self.read_marker(destination)["component"], "checkm2")

    def test_main_skips_valid_existing_destination_with_marker(self) -> None:
        """Reuse one ready destination rather than rebuilding it."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            source = self.create_taxdump_dir(tmpdir / "sources" / "taxdump")
            destination = tmpdir / "prepared" / "taxdump"

            first_exit_code, _stdout, _stderr = self.run_cli(
                [
                    "--taxdump-source",
                    str(source),
                    "--taxdump-dest",
                    str(destination),
                ]
            )
            second_exit_code, _stdout, _stderr = self.run_cli(
                [
                    "--taxdump-source",
                    str(source),
                    "--taxdump-dest",
                    str(destination),
                    "--report",
                    str(tmpdir / "report.tsv"),
                ]
            )

            self.assertEqual(first_exit_code, 0)
            self.assertEqual(second_exit_code, 0)
            rows = read_tsv(tmpdir / "report.tsv")
            self.assertEqual(rows[0]["status"], "present")

    def test_main_reports_incomplete_destination_without_marker(self) -> None:
        """Fail when one destination is valid on disk but lacks the ready marker."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            source = self.create_taxdump_dir(tmpdir / "sources" / "taxdump")
            destination = self.create_taxdump_dir(tmpdir / "prepared" / "taxdump")

            exit_code, _stdout, stderr = self.run_cli(
                [
                    "--taxdump-source",
                    str(source),
                    "--taxdump-dest",
                    str(destination),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertIn(prepare_runtime_databases.MARKER_FILE_NAME, stderr)

    def test_codetta_destination_without_ready_marker_is_incomplete(self) -> None:
        """Reject one valid Codetta directory that lacks the ready marker."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            source = self.create_codetta_dir(tmpdir / "sources" / "codetta")
            destination = self.create_codetta_dir(tmpdir / "prepared" / "codetta")

            exit_code, _stdout, stderr = self.run_cli(
                [
                    "--codetta-source",
                    str(source),
                    "--codetta-dest",
                    str(destination),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertIn("Destination is valid but incomplete for codetta", stderr)
            self.assertIn(prepare_runtime_databases.MARKER_FILE_NAME, stderr)

    def test_codetta_force_rebuilds_destination_missing_ready_marker(self) -> None:
        """Rebuild one valid-but-incomplete Codetta destination when forced."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            source = self.create_codetta_dir(tmpdir / "sources" / "codetta")
            destination = self.create_codetta_dir(tmpdir / "prepared" / "codetta")
            report_path = tmpdir / "report.tsv"

            exit_code, _stdout, _stderr = self.run_cli(
                [
                    "--codetta-source",
                    str(source),
                    "--codetta-dest",
                    str(destination),
                    "--force",
                    "--report",
                    str(report_path),
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertEqual(self.read_marker(destination)["component"], "codetta")
            self.assertEqual(read_tsv(report_path)[0]["status"], "prepared")

    def test_main_fails_on_invalid_existing_destination_without_force(self) -> None:
        """Fail cleanly when one destination exists but does not validate."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            source = self.create_padloc_dir(tmpdir / "sources" / "padloc")
            self.write_text_file(tmpdir / "prepared" / "padloc" / "broken.txt", "broken\n")

            exit_code, _stdout, stderr = self.run_cli(
                [
                    "--padloc-source",
                    str(source),
                    "--padloc-dest",
                    str(tmpdir / "prepared" / "padloc"),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertIn("invalid or partially prepared", stderr)

    def test_main_force_rebuilds_invalid_destination(self) -> None:
        """Replace one invalid destination when force is requested."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            source = self.create_padloc_dir(tmpdir / "sources" / "padloc")
            self.write_text_file(tmpdir / "prepared" / "padloc" / "broken.txt", "broken\n")

            exit_code, _stdout, _stderr = self.run_cli(
                [
                    "--padloc-source",
                    str(source),
                    "--padloc-dest",
                    str(tmpdir / "prepared" / "padloc"),
                    "--force",
                ]
            )

            self.assertEqual(exit_code, 0)
            self.assertTrue(
                (tmpdir / "prepared" / "padloc" / "hmm" / "padlocdb.hmm").is_file()
            )
            self.assertFalse((tmpdir / "prepared" / "padloc" / "broken.txt").exists())

    def test_main_honours_preparation_lock(self) -> None:
        """Fail fast when another preparation lock is already present."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            source = self.create_taxdump_dir(tmpdir / "sources" / "taxdump")
            destination = tmpdir / "prepared" / "taxdump"
            lock_path = prepare_runtime_databases.lock_path(destination)
            self.write_text_file(lock_path, "{}\n")

            exit_code, _stdout, stderr = self.run_cli(
                [
                    "--taxdump-source",
                    str(source),
                    "--taxdump-dest",
                    str(destination),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertIn("Preparation lock already exists", stderr)

    def test_main_prepares_archive_sources_without_using_system_tmp(self) -> None:
        """Extract archive sources beside the destination parent rather than via system tmp."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)

            taxdump_archive = self.create_zip_archive(
                self.create_taxdump_dir(tmpdir / "sources" / "taxdump_dir"),
                tmpdir / "archives" / "taxdump.zip",
                "taxdump_payload",
            )
            busco_archive = self.create_tar_archive(
                self.create_busco_lineage_dir(tmpdir / "sources" / "busco_dir"),
                tmpdir / "archives" / "busco.tar.gz",
                "busco_payload",
            )
            padloc_archive = self.create_zip_archive(
                self.create_padloc_dir(tmpdir / "sources" / "padloc_dir"),
                tmpdir / "archives" / "padloc.zip",
                "padloc_payload",
            )

            with mock.patch.dict(os.environ, {"TMPDIR": str(tmpdir / "blocked-tmp")}):
                exit_code, _stdout, _stderr = self.run_cli(
                    [
                        "--taxdump-source",
                        str(taxdump_archive),
                        "--taxdump-dest",
                        str(tmpdir / "prepared" / "taxdump"),
                        "--busco-lineage-source",
                        f"bacillota_odb12={busco_archive}",
                        "--busco-dest-root",
                        str(tmpdir / "prepared" / "busco"),
                        "--padloc-source",
                        str(padloc_archive),
                        "--padloc-dest",
                        str(tmpdir / "prepared" / "padloc"),
                    ]
                )

            self.assertEqual(exit_code, 0)
            self.assertEqual(list((tmpdir / "prepared").glob(".prepare-*")), [])
            self.assertTrue((tmpdir / "prepared" / "taxdump" / "names.dmp").is_file())
            self.assertTrue(
                (tmpdir / "prepared" / "busco" / "bacillota_odb12" / "dataset.cfg").is_file()
            )
            self.assertTrue(
                (tmpdir / "prepared" / "padloc" / "hmm" / "padlocdb.hmm").is_file()
            )

    def test_main_downloads_remote_default_sources_with_aria2(self) -> None:
        """Prepare all runtime databases from curated remote defaults."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            taxdump_archive = self.create_zip_archive(
                self.create_taxdump_dir(tmpdir / "fixtures" / "taxdump"),
                tmpdir / "fixtures" / "taxdump.zip",
                "taxdump_payload",
            )
            taxdump_checksum = self.write_text_file(
                tmpdir / "fixtures" / "taxdump.zip.md5",
                f"{self.checksum_md5(taxdump_archive)}  taxdump.zip\n",
            )
            checkm2_archive = self.create_tar_archive(
                self.create_checkm2_dir(tmpdir / "fixtures" / "checkm2"),
                tmpdir / "fixtures" / "checkm2.tar.gz",
                "checkm2_payload",
            )
            bacillota_archive = self.create_tar_archive(
                self.create_busco_lineage_dir(tmpdir / "fixtures" / "bacillota_odb12"),
                tmpdir / "fixtures" / "bacillota_odb12.tar.gz",
                "bacillota_payload",
            )
            myco_archive = self.create_tar_archive(
                self.create_busco_lineage_dir(tmpdir / "fixtures" / "mycoplasmatota_odb12"),
                tmpdir / "fixtures" / "mycoplasmatota_odb12.tar.gz",
                "myco_payload",
            )
            codetta_archive = self.create_tar_archive(
                self.create_codetta_dir(tmpdir / "fixtures" / "codetta"),
                tmpdir / "fixtures" / "codetta.tar.gz",
                "codetta_payload",
            )
            padloc_archive = self.create_zip_archive(
                self.create_padloc_dir(tmpdir / "fixtures" / "padloc"),
                tmpdir / "fixtures" / "padloc.zip",
                "padloc_payload",
            )
            manifest_path = self.write_remote_manifest(
                tmpdir / "remote_sources.json",
                self.build_remote_manifest(
                    taxdump_url="https://example.invalid/taxdump.zip",
                    checkm2_url="https://example.invalid/checkm2.tar.gz",
                    busco_template="https://example.invalid/{lineage}.tar.gz",
                    codetta_url="https://example.invalid/codetta.tar.gz",
                    padloc_url="https://example.invalid/padloc.zip",
                    taxdump_checksum_url="https://example.invalid/taxdump.zip.md5",
                    checkm2_checksum_value=self.checksum_md5(checkm2_archive),
                ),
            )
            source_map = {
                "https://example.invalid/taxdump.zip": taxdump_archive,
                "https://example.invalid/taxdump.zip.md5": taxdump_checksum,
                "https://example.invalid/checkm2.tar.gz": checkm2_archive,
                "https://example.invalid/bacillota_odb12.tar.gz": bacillota_archive,
                "https://example.invalid/mycoplasmatota_odb12.tar.gz": myco_archive,
                "https://example.invalid/codetta.tar.gz": codetta_archive,
                "https://example.invalid/padloc.zip": padloc_archive,
            }
            recorded_calls: list[list[str]] = []

            with (
                mock.patch.object(prepare_runtime_databases.shutil, "which", return_value="/usr/bin/aria2c"),
                mock.patch.object(
                    prepare_runtime_databases.subprocess,
                    "run",
                    self.make_fake_aria2(source_map, recorded_calls),
                ),
            ):
                exit_code, _stdout, _stderr = self.run_cli(
                    [
                        "--taxdump-dest",
                        str(tmpdir / "prepared" / "taxdump"),
                        "--checkm2-dest",
                        str(tmpdir / "prepared" / "checkm2"),
                        "--busco-dest-root",
                        str(tmpdir / "prepared" / "busco"),
                        "--codetta-dest",
                        str(tmpdir / "prepared" / "codetta"),
                        "--padloc-dest",
                        str(tmpdir / "prepared" / "padloc"),
                        "--download",
                        "--remote-source-manifest",
                        str(manifest_path),
                        "--report",
                        str(tmpdir / "report.tsv"),
                    ]
                )

            self.assertEqual(exit_code, 0)
            self.assertTrue((tmpdir / "prepared" / "taxdump" / "names.dmp").is_file())
            self.assertTrue(
                (tmpdir / "prepared" / "checkm2" / "CheckM2_database.dmnd").is_file()
            )
            self.assertTrue(
                (tmpdir / "prepared" / "busco" / "bacillota_odb12" / "dataset.cfg").is_file()
            )
            self.assertTrue(
                (tmpdir / "prepared" / "busco" / "mycoplasmatota_odb12" / "dataset.cfg").is_file()
            )
            self.assertTrue((tmpdir / "prepared" / "codetta" / "Pfam-A_enone.hmm").is_file())
            self.assertTrue(
                (tmpdir / "prepared" / "padloc" / "hmm" / "padlocdb.hmm").is_file()
            )

            taxdump_marker = self.read_marker(tmpdir / "prepared" / "taxdump")
            self.assertEqual(taxdump_marker["source_metadata"]["transport"], "aria2")
            self.assertEqual(taxdump_marker["source_metadata"]["version"], "current")
            self.assertEqual(
                taxdump_marker["source"],
                "https://example.invalid/taxdump.zip",
            )

            rows = read_tsv(tmpdir / "report.tsv")
            row_map = {row["component"]: row for row in rows}
            self.assertIn("transport=aria2", row_map["taxdump"]["details"])
            self.assertIn("source_mode=remote", row_map["taxdump"]["details"])
            self.assertEqual(row_map["codetta"]["status"], "prepared")

            first_call = recorded_calls[0]
            self.assertIn("--allow-overwrite=true", first_call)
            self.assertIn("--auto-file-renaming=false", first_call)
            self.assertIn("--dir", first_call)
            self.assertIn("--out", first_call)
            self.assertEqual(first_call[-1], "https://example.invalid/taxdump.zip")

    def test_busco_lineage_checksums_are_verified_before_preparation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            lineage = "bacillota_odb12"
            archive = self.create_tar_archive(
                self.create_busco_lineage_dir(root / "fixture"),
                root / "lineage.tar.gz", "payload",
            )
            config = {
                "kind": "lineage_archives",
                "lineage_url_template": "https://example.invalid/{lineage}.tar.gz",
                "checksums": {lineage: {"type": "md5", "value": self.checksum_md5(archive)}},
            }
            manifest = {"busco": {"default_version": "pinned", "versions": {"pinned": config}}}
            arguments = dict(
                manifest=manifest, remote_component="busco", component_label=f"busco:{lineage}",
                version="pinned", validator=prepare_runtime_databases.build_busco_lineage_validator(lineage),
                scratch_root=None, lineage=lineage,
            )
            with (
                mock.patch.object(prepare_runtime_databases.shutil, "which", return_value="aria2c"),
                mock.patch.object(prepare_runtime_databases.subprocess, "run", self.make_fake_aria2(
                    {f"https://example.invalid/{lineage}.tar.gz": archive}, [],
                )),
            ):
                _, _, metadata = prepare_runtime_databases.prepare_remote_component(
                    destination=root / "valid", **arguments,
                )
                self.assertEqual(metadata["checksum"], f"md5:{self.checksum_md5(archive)}")
                config["checksums"][lineage]["value"] = "0" * 32
                with self.assertRaises(prepare_runtime_databases.ChecksumMismatchError):
                    prepare_runtime_databases.prepare_remote_component(destination=root / "bad", **arguments)
                self.assertFalse((root / "bad").exists())
                config["checksums"] = {}
                with self.assertRaisesRegex(prepare_runtime_databases.PrepareRuntimeDatabasesError, "no checksum"):
                    prepare_runtime_databases.prepare_remote_component(destination=root / "missing", **arguments)

    def test_main_remote_version_pin_selects_the_requested_release(self) -> None:
        """Use the requested remote version rather than the default remote version."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            release_one = self.create_zip_archive(
                self.create_taxdump_dir(tmpdir / "fixtures" / "taxdump_v1"),
                tmpdir / "fixtures" / "taxdump_v1.zip",
                "taxdump_payload_v1",
            )
            release_two_dir = self.create_taxdump_dir(tmpdir / "fixtures" / "taxdump_v2")
            self.write_text_file(
                release_two_dir / "names.dmp",
                "2\t|\trelease-two\t|\t\t|\tscientific name\t|\n",
            )
            release_two = self.create_zip_archive(
                release_two_dir,
                tmpdir / "fixtures" / "taxdump_v2.zip",
                "taxdump_payload_v2",
            )
            manifest_path = self.write_remote_manifest(
                tmpdir / "remote_sources.json",
                self.build_remote_manifest(
                    taxdump_url="https://example.invalid/taxdump-v1.zip",
                    checkm2_url="https://example.invalid/checkm2.tar.gz",
                    busco_template="https://example.invalid/{lineage}.tar.gz",
                    padloc_url="https://example.invalid/padloc.zip",
                    taxdump_checksum_value=self.checksum_md5(release_one),
                    checkm2_checksum_value="unused",
                    taxdump_versions={
                        "v1": {
                            "kind": "archive",
                            "url": "https://example.invalid/taxdump-v1.zip",
                            "archive_name": "taxdump-v1.zip",
                            "checksum": {
                                "type": "md5",
                                "value": self.checksum_md5(release_one),
                            },
                        },
                        "v2": {
                            "kind": "archive",
                            "url": "https://example.invalid/taxdump-v2.zip",
                            "archive_name": "taxdump-v2.zip",
                            "checksum": {
                                "type": "md5",
                                "value": self.checksum_md5(release_two),
                            },
                        },
                    },
                ),
            )
            source_map = {
                "https://example.invalid/taxdump-v1.zip": release_one,
                "https://example.invalid/taxdump-v2.zip": release_two,
            }

            with (
                mock.patch.object(prepare_runtime_databases.shutil, "which", return_value="/usr/bin/aria2c"),
                mock.patch.object(
                    prepare_runtime_databases.subprocess,
                    "run",
                    self.make_fake_aria2(source_map, []),
                ),
            ):
                exit_code, _stdout, _stderr = self.run_cli(
                    [
                        "--taxdump-dest",
                        str(tmpdir / "prepared" / "taxdump"),
                        "--taxdump-version",
                        "v2",
                        "--download",
                        "--remote-source-manifest",
                        str(manifest_path),
                    ]
                )

            self.assertEqual(exit_code, 0)
            names_content = (tmpdir / "prepared" / "taxdump" / "names.dmp").read_text(
                encoding="utf-8"
            )
            self.assertIn("release-two", names_content)

    def test_main_prefers_local_source_over_remote_download(self) -> None:
        """Use a valid local source even when remote download is enabled."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            local_source = self.create_taxdump_dir(tmpdir / "sources" / "taxdump")
            manifest_path = self.write_remote_manifest(
                tmpdir / "remote_sources.json",
                self.build_remote_manifest(
                    taxdump_url="https://example.invalid/taxdump.zip",
                    checkm2_url="https://example.invalid/checkm2.tar.gz",
                    busco_template="https://example.invalid/{lineage}.tar.gz",
                    padloc_url="https://example.invalid/padloc.zip",
                    taxdump_checksum_value="unused",
                    checkm2_checksum_value="unused",
                ),
            )

            with (
                mock.patch.object(prepare_runtime_databases.shutil, "which", return_value="/usr/bin/aria2c"),
                mock.patch.object(prepare_runtime_databases.subprocess, "run") as mocked_run,
            ):
                exit_code, _stdout, _stderr = self.run_cli(
                    [
                        "--taxdump-source",
                        str(local_source),
                        "--taxdump-dest",
                        str(tmpdir / "prepared" / "taxdump"),
                        "--download",
                        "--remote-source-manifest",
                        str(manifest_path),
                    ]
                )

            self.assertEqual(exit_code, 0)
            mocked_run.assert_not_called()

    def test_main_checksum_mismatch_fails_remote_download(self) -> None:
        """Reject one downloaded archive when its checksum does not match."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            taxdump_archive = self.create_zip_archive(
                self.create_taxdump_dir(tmpdir / "fixtures" / "taxdump"),
                tmpdir / "fixtures" / "taxdump.zip",
                "taxdump_payload",
            )
            manifest_path = self.write_remote_manifest(
                tmpdir / "remote_sources.json",
                self.build_remote_manifest(
                    taxdump_url="https://example.invalid/taxdump.zip",
                    checkm2_url="https://example.invalid/checkm2.tar.gz",
                    busco_template="https://example.invalid/{lineage}.tar.gz",
                    padloc_url="https://example.invalid/padloc.zip",
                    taxdump_checksum_value="badchecksum",
                    checkm2_checksum_value="unused",
                ),
            )

            recorded_calls: list[list[str]] = []
            with (
                mock.patch.object(prepare_runtime_databases.shutil, "which", return_value="/usr/bin/aria2c"),
                mock.patch.object(
                    prepare_runtime_databases.subprocess,
                    "run",
                    self.make_fake_aria2(
                        {"https://example.invalid/taxdump.zip": taxdump_archive},
                        recorded_calls,
                    ),
                ),
            ):
                exit_code, _stdout, stderr = self.run_cli(
                    [
                        "--taxdump-dest",
                        str(tmpdir / "prepared" / "taxdump"),
                        "--download",
                        "--remote-source-manifest",
                        str(manifest_path),
                    ]
                )

            self.assertEqual(exit_code, 1)
            self.assertIn("Checksum mismatch", stderr)
            self.assertFalse((tmpdir / "prepared" / "taxdump").exists())
            self.assertEqual(
                sum(
                    1
                    for call in recorded_calls
                    if call[-1] == "https://example.invalid/taxdump.zip"
                ),
                prepare_runtime_databases.MAX_MD5_DOWNLOAD_ATTEMPTS,
            )

    def test_main_checksum_mismatch_retries_and_recovers(self) -> None:
        """Accept one remote archive when a later MD5 retry succeeds."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            bad_taxdump_dir = tmpdir / "fixtures" / "taxdump_bad"
            self.write_text_file(
                bad_taxdump_dir / "names.dmp",
                "1\t|\tbad root\t|\t\t|\tscientific name\t|\n",
            )
            self.write_text_file(
                bad_taxdump_dir / "nodes.dmp",
                "1\t|\t1\t|\tno rank\t|\n",
            )
            bad_archive = self.create_zip_archive(
                bad_taxdump_dir,
                tmpdir / "fixtures" / "taxdump_bad.zip",
                "taxdump_payload",
            )
            good_archive = self.create_zip_archive(
                self.create_taxdump_dir(tmpdir / "fixtures" / "taxdump_good"),
                tmpdir / "fixtures" / "taxdump_good.zip",
                "taxdump_payload",
            )
            manifest_path = self.write_remote_manifest(
                tmpdir / "remote_sources.json",
                self.build_remote_manifest(
                    taxdump_url="https://example.invalid/taxdump.zip",
                    checkm2_url="https://example.invalid/checkm2.tar.gz",
                    busco_template="https://example.invalid/{lineage}.tar.gz",
                    padloc_url="https://example.invalid/padloc.zip",
                    taxdump_checksum_value=self.checksum_md5(good_archive),
                    checkm2_checksum_value="unused",
                ),
            )

            recorded_calls: list[list[str]] = []
            download_attempts = 0

            def fake_run(
                args: list[str],
                *,
                check: bool = False,
                capture_output: bool = True,
                text: bool = True,
            ) -> subprocess.CompletedProcess[str]:
                """Return corrupt downloads first, then one valid archive."""
                nonlocal download_attempts
                self.assertEqual(args[0], "aria2c")
                self.assertFalse(check)
                self.assertTrue(capture_output)
                self.assertTrue(text)
                recorded_calls.append(args)
                url = args[-1]
                destination_dir = Path(args[args.index("--dir") + 1])
                destination_name = args[args.index("--out") + 1]
                destination_dir.mkdir(parents=True, exist_ok=True)
                if url == "https://example.invalid/taxdump.zip":
                    download_attempts += 1
                    source = bad_archive if download_attempts < 3 else good_archive
                    shutil.copy2(source, destination_dir / destination_name)
                    return subprocess.CompletedProcess(args, 0, "", "")
                raise AssertionError(f"Unexpected URL: {url}")

            with (
                mock.patch.object(prepare_runtime_databases.shutil, "which", return_value="/usr/bin/aria2c"),
                mock.patch.object(
                    prepare_runtime_databases.subprocess,
                    "run",
                    mock.Mock(side_effect=fake_run),
                ),
            ):
                exit_code, stdout, stderr = self.run_cli(
                    [
                        "--taxdump-dest",
                        str(tmpdir / "prepared" / "taxdump"),
                        "--download",
                        "--remote-source-manifest",
                        str(manifest_path),
                    ]
                )

            self.assertEqual(exit_code, 0)
            self.assertIn("--taxdump", stdout)
            self.assertIn("Checksum mismatch", stderr)
            self.assertEqual(download_attempts, 3)
            self.assertTrue(
                (
                    tmpdir
                    / "prepared"
                    / "taxdump"
                    / prepare_runtime_databases.MARKER_FILE_NAME
                ).is_file()
            )

    def test_main_remote_download_failure_propagates_aria2_error(self) -> None:
        """Surface aria2 failures when one remote download cannot be retrieved."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            manifest_path = self.write_remote_manifest(
                tmpdir / "remote_sources.json",
                self.build_remote_manifest(
                    taxdump_url="https://example.invalid/taxdump.zip",
                    checkm2_url="https://example.invalid/checkm2.tar.gz",
                    busco_template="https://example.invalid/{lineage}.tar.gz",
                    padloc_url="https://example.invalid/padloc.zip",
                    taxdump_checksum_value="unused",
                    checkm2_checksum_value="unused",
                ),
            )

            with (
                mock.patch.object(prepare_runtime_databases.shutil, "which", return_value="/usr/bin/aria2c"),
                mock.patch.object(
                    prepare_runtime_databases.subprocess,
                    "run",
                    self.make_fake_aria2(
                        {},
                        [],
                        failures={"https://example.invalid/taxdump.zip": "404 Not Found"},
                    ),
                ),
            ):
                exit_code, _stdout, stderr = self.run_cli(
                    [
                        "--taxdump-dest",
                        str(tmpdir / "prepared" / "taxdump"),
                        "--download",
                        "--remote-source-manifest",
                        str(manifest_path),
                    ]
                )

            self.assertEqual(exit_code, 1)
            self.assertIn("404 Not Found", stderr)


if __name__ == "__main__":
    unittest.main()
