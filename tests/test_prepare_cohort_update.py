"""Fail-closed validation and identity tests for published cohort reuse."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import build_master_table  # noqa: E402
import collect_versions  # noqa: E402
import master_table_contract  # noqa: E402
import prepare_cohort_update as update  # noqa: E402
import validate_inputs  # noqa: E402

LINEAGE = "bacillota_odb12"
BUSCO = "C:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200"


def write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


def table(
    path: Path, header: list[str] | tuple[str, ...], rows: list[dict[str, str]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    update.write_table(path, header, rows)


class PrepareCohortUpdateTestCase(unittest.TestCase):
    """Exercise source corruption, sample identity and deterministic snapshots."""

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory(
            prefix="cohort-preflight-", dir="/tmp"
        )
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.source = self.root / "source"
        self.genomes: dict[str, Path] = {}

    def genome(self, accession: str) -> Path:
        if accession not in self.genomes:
            path = self.root / "inputs" / f"{accession}.fasta"
            write(path, ">contig1 description\n" + "ACGT" * 400 + "\n")
            self.genomes[accession] = path
        return self.genomes[accession]

    def published_cohort(
        self, accessions: list[str], *, failed: bool = False, skipped: bool = False
    ) -> None:
        internal_ids = validate_inputs.add_collision_suffixes(accessions)
        manifest = []
        masters = []
        statuses = []
        for accession in accessions:
            internal_id = internal_ids[accession]
            row = {
                "accession": accession,
                "internal_id": internal_id,
                "is_new": "false",
                "assembly_level": "NA",
                "genome_fasta": str(self.genome(accession)),
            }
            manifest.append(row)
            sample = self.source / "samples" / accession
            write(
                sample / f"staged/{internal_id}.fasta",
                self.genome(accession).read_text(),
            )
            for relative in (
                "barrnap/rrna.gff",
                "barrnap/rrna.fa",
                "barrnap/barrnap.log",
                "codetta/codetta.log",
                "checkm2_gcode4/quality_report.tsv",
                "checkm2_gcode4/checkm2.log",
                "checkm2_gcode11/quality_report.tsv",
                "checkm2_gcode11/checkm2.log",
            ):
                write(
                    sample / relative,
                    "exit_code=0\n" if relative.endswith(".log") else "fixture\n",
                )
            write(
                sample / "16s/best_16S.fna",
                ">16S_rRNA::contig1:0-1500(+)\n" + "ACGT" * 375 + "\n",
            )
            sixteen_s = {
                "accession": accession,
                "16S": "Yes",
                "best_16S_header": "16S_rRNA::contig1:0-1500(+)",
                "best_16S_length": "1500",
                "warnings": "",
            }
            table(sample / "16s/16S_status.tsv", list(sixteen_s), [sixteen_s])
            qc = {column: "1" for column in build_master_table.CHECKM2_COLUMNS}
            qc.update(accession=accession, Gcode="4", Low_quality="false", warnings="")
            table(sample / "checkm2/checkm2_summary.tsv", list(qc), [qc])
            for code in (4, 11):
                report = {
                    column.removesuffix(f"_gcode{code}"): value
                    for column, value in qc.items()
                    if column.endswith(f"_gcode{code}")
                }
                report["Name"] = internal_id
                table(
                    sample / f"checkm2_gcode{code}/quality_report.tsv",
                    list(report),
                    [report],
                )
            codetta = {
                "accession": accession,
                "Codetta_Genetic_Code": "F" * 64,
                "Codetta_NCBI_Table_Candidates": "NA",
                "codetta_status": "done",
                "warnings": "",
            }
            table(sample / "codetta/codetta_summary.tsv", list(codetta), [codetta])
            busco = {
                "accession": accession,
                "lineage": LINEAGE,
                f"BUSCO_{LINEAGE}": BUSCO,
                "busco_status": "done",
                "warnings": "",
            }
            table(
                sample / f"busco/{LINEAGE}/busco_summary_{LINEAGE}.tsv",
                list(busco),
                [busco],
            )
            write(
                sample / f"busco/{LINEAGE}/short_summary.json",
                '{"results": {"C": 98, "S": 98, "D": 0, "F": 1, "M": 1, "n": 200}}\n',
            )
            write(sample / f"busco/{LINEAGE}/busco.log", "exit_code=0\n")
            ccfinder = {
                "accession": accession,
                "CRISPRS": "0",
                "SPACERS_SUM": "0",
                "CRISPR_FRAC": "0",
                "ccfinder_status": "done",
                "warnings": "",
            }
            table(sample / "ccfinder/ccfinder_strains.tsv", list(ccfinder), [ccfinder])
            for relative in (
                "ccfinder/result.json",
                "ccfinder/ccfinder.log",
                "ccfinder/ccfinder_contigs.tsv",
                "ccfinder/ccfinder_crisprs.tsv",
                "padloc/padloc/results.tsv",
            ):
                write(sample / relative, "fixture\n")
            write(
                sample / "ccfinder/result.json",
                '{"Sequences": [{"Id": "contig1", "Length": 1600, "Crisprs": []}]}\n',
            )
            for extension in ("gff", "faa", "gbk"):
                write(
                    sample / f"prokka/prokka.{extension}", "" if failed else "fixture\n"
                )
            write(sample / "prokka/prokka.log", f"exit_code={1 if failed else 0}\n")
            write(sample / "padloc/padloc.log", "exit_code=0\n")
            if not skipped:
                write(sample / "eggnog/eggnog_annotations.tsv", "fixture\n")
                write(sample / "eggnog/eggnog.log", "exit_code=0\n")
                (sample / "eggnog/eggnog").mkdir()
            status = {
                column: "done" if column.endswith("_status") else ""
                for column in master_table_contract.build_sample_status_columns(
                    [LINEAGE]
                )
            }
            status.update(
                accession=accession,
                internal_id=internal_id,
                is_new="false",
                gcode="4",
                low_quality="false",
                prokka_status="failed" if failed else "done",
                eggnog_status="skipped" if skipped else "done",
                ani_included="true",
            )
            statuses.append(status)
            master = {
                column: "NA"
                for column in master_table_contract.build_append_columns([LINEAGE])
            }
            master.update(
                {
                    key: value
                    for summary in (qc, sixteen_s, codetta, ccfinder, busco)
                    for key, value in summary.items()
                    if key in master
                }
            )
            master.update(accession=accession, is_new="false")
            masters.append(master)
        table(
            self.source / "tables/validated_samples.tsv",
            [*validate_inputs.REQUIRED_SAMPLE_COLUMNS, "internal_id"],
            manifest,
        )
        table(
            self.source / "tables/sample_status.tsv",
            master_table_contract.build_sample_status_columns([LINEAGE]),
            statuses,
        )
        table(
            self.source / "tables/master_table.tsv",
            ["accession", *master_table_contract.build_append_columns([LINEAGE])],
            masters,
        )
        table(
            self.source / "tables/tool_and_db_versions.tsv",
            collect_versions.OUTPUT_COLUMNS,
            [
                {
                    "component": "checkm2",
                    "kind": "tool",
                    "version": "test",
                    "image_or_path": "NA",
                    "notes": "published source",
                }
            ],
        )

    def arguments(
        self, accessions: list[str], *, lineages: list[str] | None = None
    ) -> argparse.Namespace:
        lineages = [LINEAGE] if lineages is None else lineages
        manifest = self.root / "current.csv"
        with manifest.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(validate_inputs.REQUIRED_SAMPLE_COLUMNS)
            for accession in accessions:
                writer.writerow([accession, "false", "NA", self.genome(accession)])
        metadata = self.root / "metadata.tsv"
        table(
            metadata,
            ["accession", "Assembly_Level"],
            [
                {"accession": accession, "Assembly_Level": "Scaffold"}
                for accession in accessions
            ],
        )
        validated = self.root / "validated"
        validate_inputs.run_validation(
            manifest, metadata, validated, None, busco_lineages=lineages
        )
        inputs = self.root / "candidate_genomes"
        inputs.mkdir(exist_ok=True)
        for index, accession in enumerate(accessions, start=1):
            shutil.copy2(self.genome(accession), inputs / f"genome{index:06d}")
        settings = self.root / "settings.json"
        settings.write_text('{"ani_threshold": 0.95}\n')
        return argparse.Namespace(
            source_results=self.source,
            destination=self.root / "results",
            validated_samples=validated / "validated_samples.tsv",
            accession_map=validated / "accession_map.tsv",
            initial_status=validated / "sample_status.tsv",
            validation_warnings=validated / "validation_warnings.tsv",
            genome_inputs=inputs,
            metadata=metadata,
            settings=settings,
            busco_lineage=lineages,
            previous_update=None,
            outdir=self.root / "prepared",
        )

    def test_add_remove_uses_membership_instead_of_is_new(self) -> None:
        self.published_cohort(["A", "B"])
        args = self.arguments(["B", "C"])
        update.run_prepare(args)
        _, audit = update.read_table(args.outdir / "cohort_update.tsv")
        self.assertEqual(
            {row["accession"]: row["action"] for row in audit},
            {"A": "removed", "B": "reused", "C": "added"},
        )
        _, new = update.read_table(args.outdir / "new_samples.tsv")
        self.assertEqual(
            [(row["accession"], row["is_new"]) for row in new], [("C", "false")]
        )
        self.assertTrue(
            all(
                len(row["genome_sha256"]) == 64
                for row in audit
                if row["action"] != "removed"
            )
        )

    def test_replacement_cohort_does_not_require_old_busco_lineages(self) -> None:
        self.published_cohort(["A"])
        args = self.arguments(["B"], lineages=["another_odb12"])
        update.run_prepare(args)
        _, new = update.read_table(args.outdir / "new_samples.tsv")
        self.assertEqual([row["accession"] for row in new], ["B"])

    def test_supplemental_columns_cannot_override_internal_reuse_fields(self) -> None:
        self.published_cohort(["A"])
        args = self.arguments(["A"])
        header, rows = update.read_table(args.validated_samples)
        rows[0]["source_gcode"] = "supplemental metadata"
        table(args.validated_samples, [*header, "source_gcode"], rows)
        update.run_prepare(args)
        _, validated = update.read_table(args.outdir / "validated_samples.tsv")
        self.assertEqual(validated[0]["source_gcode"], "supplemental metadata")
        header, reused = update.read_table(args.outdir / "reused_samples.tsv")
        self.assertEqual(len(header), len(set(header)))
        self.assertEqual(reused[0]["source_gcode"], "4")

    def test_preserves_published_ids_when_adding_a_sanitisation_collision(self) -> None:
        self.published_cohort(["ID-X"])
        args = self.arguments(["ID-X", "ID.X"])
        update.run_prepare(args)
        _, rows = update.read_table(args.outdir / "validated_samples.tsv")
        ids = {row["accession"]: row["internal_id"] for row in rows}
        self.assertEqual(ids["ID-X"], "ID_X")
        self.assertTrue(ids["ID.X"].startswith("ID_X_"))
        _, initial = update.read_table(args.outdir / "sample_status.tsv")
        self.assertEqual({row["accession"]: row["internal_id"] for row in initial}, ids)

    def test_missing_artefact_blocks_all_additions(self) -> None:
        self.published_cohort(["A"])
        args = self.arguments(["A", "B"])
        (self.source / "samples/A/checkm2/checkm2_summary.tsv").unlink()
        with self.assertRaisesRegex(
            update.CohortUpdateError, "Missing published artefact"
        ):
            update.run_prepare(args)
        self.assertFalse(args.outdir.exists())

    def test_changed_retained_sequence_fails_before_output(self) -> None:
        self.published_cohort(["A"])
        args = self.arguments(["A"])
        write(args.genome_inputs / "genome000001", ">contig1\nTTTT\n")
        with self.assertRaisesRegex(update.CohortUpdateError, "Genome content changed"):
            update.run_prepare(args)
        self.assertFalse(args.outdir.exists())

    def test_relocated_rewrapped_gzip_input_matches_published_sequence(self) -> None:
        self.published_cohort(["A"])
        args = self.arguments(["A"])
        with gzip.open(args.genome_inputs / "genome000001", "wt") as handle:
            handle.write(">contig1 another description\n" + "ACGT\n" * 400)
        update.run_prepare(args)

    def test_soft_masking_change_is_not_silently_reused(self) -> None:
        self.published_cohort(["A"])
        args = self.arguments(["A"])
        write(args.genome_inputs / "genome000001", ">contig1\n" + "acgt" * 400 + "\n")
        with self.assertRaisesRegex(update.CohortUpdateError, "Genome content changed"):
            update.run_prepare(args)

    def test_documented_failed_and_skipped_annotations_are_reusable(self) -> None:
        self.published_cohort(["A"], failed=True, skipped=True)
        args = self.arguments(["A"])
        update.run_prepare(args)
        _, rows = update.read_table(args.outdir / "reused_samples.tsv")
        self.assertEqual(rows[0]["source_eggnog_status"], "skipped")

    def test_source_summary_must_agree_with_published_status(self) -> None:
        self.published_cohort(["A"])
        args = self.arguments(["A"])
        path = self.source / "samples/A/prokka/prokka.log"
        path.write_text("exit_code=1\n")
        with self.assertRaisesRegex(
            update.CohortUpdateError, "prokka_status disagrees"
        ):
            update.run_prepare(args)

    def test_source_summary_rejects_duplicate_accessions_and_malformed_values(
        self,
    ) -> None:
        self.published_cohort(["A"])
        args = self.arguments(["A"])
        path = self.source / "samples/A/checkm2/checkm2_summary.tsv"
        header, rows = update.read_table(path)
        table(path, header, [rows[0], rows[0]])
        with self.assertRaisesRegex(update.CohortUpdateError, "exactly one row"):
            update.run_prepare(args)
        rows[0]["Completeness_gcode4"] = "nan"
        table(path, header, rows)
        with self.assertRaisesRegex(
            update.CohortUpdateError, "Invalid published Completeness"
        ):
            update.run_prepare(args)

    def test_missing_requested_busco_lineage_fails(self) -> None:
        self.published_cohort(["A"])
        args = self.arguments(["A"])
        args.busco_lineage = ["another_odb12"]
        with self.assertRaisesRegex(update.CohortUpdateError, "BUSCO lineages absent"):
            update.run_prepare(args)

    def test_removed_samples_do_not_require_their_artefacts(self) -> None:
        self.published_cohort(["A", "B"])
        args = self.arguments(["B"])
        shutil.rmtree(self.source / "samples/A")
        update.run_prepare(args)

    def test_source_and_destination_must_not_overlap(self) -> None:
        self.published_cohort(["A"])
        args = self.arguments(["A"])
        for destination in (self.source, self.source / "child", self.root):
            args.destination = destination
            with self.assertRaisesRegex(update.CohortUpdateError, "non-overlapping"):
                update.run_prepare(args)

    def test_nested_symlink_cannot_reintroduce_work_dependencies(self) -> None:
        self.published_cohort(["A"])
        args = self.arguments(["A"])
        (self.source / "samples/A/dangling-work-file").symlink_to(
            self.root / "missing-work"
        )
        with self.assertRaisesRegex(update.CohortUpdateError, "is a symlink"):
            update.run_prepare(args)

    def test_resume_requires_the_same_inputs_and_settings(self) -> None:
        self.published_cohort(["A"])
        args = self.arguments(["A"])
        update.run_prepare(args)
        args.previous_update = args.outdir / "cohort_update_run.json"
        update.run_prepare(args)
        args.settings.write_text(json.dumps({"ani_threshold": 0.99}))
        with self.assertRaisesRegex(
            update.CohortUpdateError, "different inputs or settings"
        ):
            update.run_prepare(args)

    def test_fingerprint_rejects_empty_duplicate_and_malformed_records(self) -> None:
        path = self.root / "invalid.fasta"
        for content in ("", ">id\n", ">id\nACGT\n>id\nACGT\n", "ACGT\n", ">id\nBAD?\n"):
            path.write_text(content)
            with self.assertRaises(update.CohortUpdateError):
                update.genome_fingerprint(path)


if __name__ == "__main__":
    unittest.main()
