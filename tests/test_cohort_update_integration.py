"""Run the real cohort workflow with input-driven external-tool substitutes."""

from __future__ import annotations

import csv
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")
TOOLS = {
    "barrnap",
    "checkm2",
    "busco",
    "codetta",
    "prokka",
    "ccfinder",
    "padloc",
    "eggnog",
    "fastani",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fingerprints(root: Path) -> dict[str, str]:
    return {
        str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in root.rglob("*")
        if path.is_file()
    }


@unittest.skipUnless(
    NEXTFLOW, "nextflow is required for cohort-update integration tests"
)
class CohortUpdateIntegrationTestCase(unittest.TestCase):
    """Prove reuse, cohort membership, scientific joins and repeated updates."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.root = Path(
            tempfile.mkdtemp(prefix="nf-annotation-cohort-update-", dir="/tmp")
        )
        cls.project = cls.root / "pipeline"
        cls.project.mkdir()
        for directory in ("bin", "modules", "subworkflows", "conf", "assets"):
            shutil.copytree(ROOT / directory, cls.project / directory)
        for name in ("main.nf", "nextflow.config"):
            shutil.copy2(ROOT / name, cls.project / name)
        shutil.copytree(
            ROOT / "external/nf-helper/conf", cls.project / "external/nf-helper/conf"
        )
        (cls.project / "tests").mkdir()
        shutil.copy2(
            ROOT / "tests/cohort_update_fixture.py",
            cls.project / "tests/cohort_update_fixture.py",
        )
        # Keep production scripts intact. Replace expensive stubs and execute
        # the real lightweight scripts in stub mode rather than canned tables.
        for path in (cls.project / "modules/local").glob("*.nf"):
            content = path.read_text()
            if "\n    stub:" not in content:
                continue
            content = content.split("\n    stub:", 1)[0]
            tool = path.stem
            if tool in TOOLS - {"fastani"}:
                arguments = "--process-name '${task.process}'"
                arguments += ' --accession "${meta.accession}" --internal-id "${meta.internal_id}"'
                if tool in {
                    "barrnap",
                    "checkm2",
                    "busco",
                    "codetta",
                    "prokka",
                    "ccfinder",
                }:
                    arguments += ' --genome "${genome}"'
                if tool == "checkm2":
                    arguments += ' --translation-table "${translation_table}"'
                if tool == "busco":
                    arguments += ' --lineage "${lineage}"'
                content += (
                    '\n    stub:\n    """\n    python3 "${projectDir}/tests/cohort_update_fixture.py" '
                    + tool
                    + " "
                    + arguments
                    + '\n    """\n'
                )
            path.write_text(content + "}\n")
        for executable, function in (
            ("seqtk", "fake_seqtk"),
            ("fastANI", "fake_fastani"),
        ):
            wrapper = cls.project / "bin" / executable
            wrapper.write_text(
                f"#!{sys.executable}\nimport sys\nsys.path.insert(0, {str(cls.project / 'tests')!r})\nfrom cohort_update_fixture import {function}\nsys.exit({function}(sys.argv[1:]))\n"
            )
            wrapper.chmod(0o755)
        cls.inputs = cls.root / "inputs"
        cls.inputs.mkdir()
        cls.genomes = {}
        for index, accession in enumerate(
            ("A", "B", "C", "D", "LOW", "NOGCODE", "FAILED", "ID-X", "ID.X")
        ):
            genome = cls.inputs / f"{accession}.fasta"
            genome.write_text(
                f">contig_{accession}\n" + "T" * index + ("ACGT" * 500)[index:] + "\n"
            )
            cls.genomes[accession] = genome
        cls.metadata = cls.inputs / "metadata.tsv"
        with cls.metadata.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "accession",
                    "Tax_ID",
                    "Organism_Name",
                    "Assembly_Level",
                    "Atypical_Warnings",
                ]
            )
            for accession in cls.genomes:
                writer.writerow(
                    [accession, "8", "Synthetic bacterium", "Scaffold", "NA"]
                )
        cls.input_fingerprints = {
            accession: hashlib.sha256(path.read_bytes()).hexdigest()
            for accession, path in cls.genomes.items()
        }

    @classmethod
    def tearDownClass(cls) -> None:
        if os.environ.get("NF_ANNOTATION_KEEP_UPDATE_TESTS"):
            print(f"Cohort-update test artefacts retained at {cls.root}")
        else:
            shutil.rmtree(cls.root)

    def manifest(self, name: str, accessions: list[str]) -> Path:
        path = self.inputs / f"{name}.csv"
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["accession", "is_new", "assembly_level", "genome_fasta"])
            for accession in accessions:
                writer.writerow([accession, "false", "NA", self.genomes[accession]])
        return path

    def run_pipeline(
        self,
        name: str,
        accessions: list[str],
        *,
        source: Path | None = None,
        extra: list[str] = (),
        resume: bool = False,
        expected_code: int = 0,
    ) -> Path:
        launch = self.root / f"launch-{name}"
        launch.mkdir(exist_ok=True)
        output = self.root / f"results-{name}"
        work = self.root / f"work-{name}"
        command = [
            NEXTFLOW,
            "run",
            str(self.project),
            "-profile",
            "test",
            "-stub-run",
            "-work-dir",
            str(work),
            "--outdir",
            str(output),
            "--sample_csv",
            str(self.manifest(name, accessions)),
            "--metadata",
            str(self.metadata),
            "--task_attempts",
            "1",
            *extra,
        ]
        if source:
            command.extend(["--update_from", str(source)])
        if resume:
            command.append("-resume")
        environment = os.environ.copy()
        environment.update(NXF_ANSI_LOG="false", NXF_DISABLE_CHECK_LATEST="true")
        result = subprocess.run(
            command,
            cwd=launch,
            env=environment,
            text=True,
            capture_output=True,
            timeout=180,
        )
        (launch / "captured.log").write_text(result.stdout + result.stderr)
        if expected_code == 0:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertTrue(
                (output / "tables/tool_and_db_versions.tsv").is_file(),
                result.stdout + result.stderr,
            )
            self.assertEqual(
                {
                    accession: hashlib.sha256(path.read_bytes()).hexdigest()
                    for accession, path in self.genomes.items()
                },
                self.input_fingerprints,
            )
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return output

    def assert_members(self, output: Path, accessions: set[str]) -> None:
        for name in ("validated_samples.tsv", "sample_status.tsv", "master_table.tsv"):
            rows = read_rows(output / "tables" / name)
            self.assertEqual({row["accession"] for row in rows}, accessions, name)
            self.assertEqual(len(rows), len(accessions), name)
        self.assertEqual(
            {path.name for path in (output / "samples").iterdir()}, accessions
        )

    def heavy_tasks(self, output: Path) -> list[dict[str, str]]:
        return [
            row
            for row in read_rows(output / "pipeline_info/trace.tsv")
            if row["name"].split(":")[-1].split(" (")[0].lower().split("_gcode")[0]
            in TOOLS - {"fastani"}
        ]

    def test_add_remove_and_successive_updates_without_original_work(self) -> None:
        source = self.run_pipeline("ab", ["A", "B"])
        self.assert_members(source, {"A", "B"})
        before = fingerprints(source)
        self.run_pipeline(
            "invalid-threshold",
            ["B", "C"],
            source=source,
            extra=["--ani_threshold", "0"],
            expected_code=1,
        )
        invalid_log = (self.root / "launch-invalid-threshold/captured.log").read_text()
        self.assertIn("ani_threshold", invalid_log)
        self.assertNotIn("Submitted process", invalid_log)
        incomplete_source = self.root / "incomplete-source"
        shutil.copytree(source, incomplete_source)
        (incomplete_source / "samples/B/checkm2/checkm2_summary.tsv").unlink()
        rejected = self.run_pipeline(
            "incomplete", ["B", "C"], source=incomplete_source, expected_code=1
        )
        self.assertEqual(self.heavy_tasks(rejected), [])
        self.assertIn(
            "Missing published artefact",
            (self.root / "launch-incomplete/captured.log").read_text(),
        )
        shutil.rmtree(self.root / "work-ab")
        shutil.rmtree(self.root / "launch-ab/.nextflow")
        updated = self.run_pipeline("bc", ["B", "C"], source=source)
        self.assert_members(updated, {"B", "C"})
        tasks = self.heavy_tasks(updated)
        self.assertEqual(len(tasks), 10)
        self.assertTrue(all("(C" in row["name"] for row in tasks), tasks)
        self.assertEqual(fingerprints(source), before)
        self.assertEqual(
            fingerprints(source / "samples/B"), fingerprints(updated / "samples/B")
        )
        actions = {
            row["accession"]: row["action"]
            for row in read_rows(updated / "tables/cohort_update.tsv")
        }
        self.assertEqual(actions, {"A": "removed", "B": "reused", "C": "added"})
        ani = read_rows(updated / "cohort/ani_clusters/ani_summary.tsv")
        self.assertEqual({row["Accession"] for row in ani}, {"B", "C"})
        self.assertEqual(
            {
                row["accession"]
                for row in read_rows(updated / "cohort/16s/all_best_16S_manifest.tsv")
            },
            {"B", "C"},
        )
        versions = read_rows(updated / "tables/tool_and_db_versions.tsv")
        self.assertTrue(
            any(row["notes"].startswith("reused source ") for row in versions)
        )
        fresh = self.run_pipeline("bc-fresh", ["B", "C"])
        self.assertEqual(
            read_rows(updated / "tables/master_table.tsv"),
            read_rows(fresh / "tables/master_table.tsv"),
        )
        self.assertEqual(
            read_rows(updated / "tables/sample_status.tsv"),
            read_rows(fresh / "tables/sample_status.tsv"),
        )
        self.assertEqual(ani, read_rows(fresh / "cohort/ani_clusters/ani_summary.tsv"))
        reclustered = self.run_pipeline(
            "high-threshold",
            ["B", "C"],
            source=updated,
            extra=[
                "--ani_threshold",
                "0.9999",
                "--ani_score_profile",
                "mag",
                "--busco_lineages",
                "mycoplasmatota_odb12",
            ],
        )
        self.assertEqual(self.heavy_tasks(reclustered), [])
        self.assertEqual(
            len(read_rows(reclustered / "cohort/ani_clusters/ani_representatives.tsv")),
            2,
        )
        shutil.rmtree(self.root / "work-bc")
        shutil.rmtree(self.root / "launch-bc/.nextflow")
        # Remove the source results too: copied B must no longer depend on it.
        shutil.rmtree(source)
        later = self.run_pipeline("bcd", ["B", "C", "D"], source=updated)
        self.assert_members(later, {"B", "C", "D"})
        self.assertTrue(all("(D" in row["name"] for row in self.heavy_tasks(later)))
        removed = self.run_pipeline(
            "b-only",
            ["B"],
            source=later,
            extra=[
                "--checkm2_db",
                "/missing/checkm2",
                "--codetta_db",
                "/missing/codetta",
                "--eggnog_db",
                "/missing/eggnog",
                "--busco_db",
                "/missing/busco",
                "--prepare_busco_datasets",
                "true",
            ],
        )
        self.assert_members(removed, {"B"})
        self.assertEqual(self.heavy_tasks(removed), [])
        self.assertEqual(
            len(read_rows(removed / "cohort/ani_clusters/ani_representatives.tsv")), 1
        )

    def test_stored_failures_skips_empty_ani_and_collision_ids(self) -> None:
        source = self.run_pipeline(
            "states",
            ["LOW", "NOGCODE", "FAILED", "ID-X"],
            extra=["--eggnog_only_accessions", "ID-X"],
        )
        initial = {
            row["accession"]: row
            for row in read_rows(source / "tables/sample_status.tsv")
        }
        self.assertEqual(initial["FAILED"]["prokka_status"], "failed")
        self.assertEqual(initial["FAILED"]["eggnog_status"], "skipped")
        updated = self.run_pipeline(
            "states-update", ["LOW", "NOGCODE", "FAILED", "ID-X", "ID.X"], source=source
        )
        final = {
            row["accession"]: row
            for row in read_rows(updated / "tables/sample_status.tsv")
        }
        for accession in initial:
            self.assertEqual(
                final[accession]["internal_id"], initial[accession]["internal_id"]
            )
            for column in initial[accession]:
                if column.endswith("_status"):
                    self.assertEqual(
                        final[accession][column],
                        initial[accession][column],
                        (accession, column),
                    )
        self.assertNotEqual(final["ID-X"]["internal_id"], final["ID.X"]["internal_id"])
        empty = self.run_pipeline("empty-ani", ["LOW", "NOGCODE"], source=updated)
        self.assert_members(empty, {"LOW", "NOGCODE"})
        self.assertEqual(read_rows(empty / "cohort/ani_clusters/ani_summary.tsv"), [])
        self.assertEqual(self.heavy_tasks(empty), [])
        same = self.run_pipeline("unchanged", ["LOW", "NOGCODE"], source=empty)
        self.assertEqual(self.heavy_tasks(same), [])
        resumed = self.run_pipeline(
            "unchanged", ["LOW", "NOGCODE"], source=empty, resume=True
        )
        self.assert_members(resumed, {"LOW", "NOGCODE"})
        self.assertTrue(
            any(
                row["status"] == "CACHED"
                for row in read_rows(resumed / "pipeline_info/trace.tsv")
            )
        )


if __name__ == "__main__":
    unittest.main()
