"""Task reuse, failure propagation, and conserved cohort counts for v0.4."""

from __future__ import annotations

import os
import subprocess
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import test_annotation_bundle as bundle_fixture
from aggregate_annotations import aggregate
from annotation_commands import commands, shell_script
from annotation_common import (
    TOOLS,
    AnnotationError,
    bundle_proteins,
    digest,
    identity,
    read_tsv,
    write_json,
    write_tsv,
)
from annotation_result import inventory, validate_result
from annotation_summary import (
    FIELDS,
    MATRICES,
    accepted_count,
    feature_counts,
    genome_summary,
)
from annotation_tasks import (
    POLICY,
    code_identity,
    plan_task,
    runtime_identity,
    task_identity,
)
from prepare_annotation_tasks import plan


class AnnotationTaskTests(unittest.TestCase):
    def setUp(self):
        self.fixture = bundle_fixture.AnnotationBundleTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root
        self.bundle = self.fixture.build("A")
        self.metadata, self.proteins = bundle_proteins(self.bundle)

    def entry(self, tool="kofam"):
        search = dict(
            tool=tool,
            runtime_id="sha256:" + "1" * 64,
            resource_id="resource1",
            command=commands(tool, 2, 32),
            native_code={},
        )
        interpretation = dict(
            schema_version=1, code=code_identity(tool), policy=POLICY[tool]
        )
        return dict(
            search_method=search,
            interpretation=interpretation,
            method_id=identity(dict(search=search, interpretation=interpretation)),
            resource=dict(
                resource_id="resource1", contract=dict(version="fixture", settings={})
            ),
            resource_manifest_sha256="resource_manifest",
            container="sha256:" + "1" * 64,
        )

    def result(self, tool="kofam", features=None, errors=None):
        record = task_identity(self.metadata, tool, self.entry(tool))
        gene = self.proteins[0]["gene_id"]
        record.update(
            status="success",
            reason="validated_native_output",
            action="run",
            exit_code=0,
            reported_hits=2,
            accepted_proteins=1,
            evidence=dict(
                features=features or {"ko": {gene: ["K00001", "K00001", "K00002"]}},
                definitions={},
                mapped_genes=[gene],
                accepted_genes=[gene],
                field_errors=errors or [],
            ),
        )
        return record

    def write_result(self, record, name):
        root = self.root / name
        (root / "raw").mkdir(parents=True)
        (root / "raw" / "native.txt").write_text("native evidence\n")
        (root / "normalized").mkdir()
        write_tsv(
            root / "normalized" / "evidence.tsv",
            ("accession", "gene_id"),
            [dict(accession=record["accession"], gene_id=self.proteins[0]["gene_id"])],
        )
        record["raw_files"] = inventory(root / "raw")
        record["normalized_files"] = inventory(root / "normalized")
        record["result_id"] = identity(record)
        write_json(root / "result.json", record)
        return root

    def test_immutable_runtime_and_allocated_memory(self):
        with self.assertRaisesRegex(AnnotationError, "Mutable"):
            runtime_identity("tool:latest")
        self.assertEqual(
            runtime_identity("repo/tool@sha256:" + "a" * 64), "sha256:" + "a" * 64
        )
        sif = self.root / "tool.sif"
        sif.write_bytes(b"immutable runtime fixture")
        self.assertTrue(runtime_identity(str(sif)).startswith("sif-sha256:"))
        with self.assertRaisesRegex(AnnotationError, "memory"):
            commands("eggnog", 8, 8)
        command = commands("eggnog", 4, 32)["native_stages"]["search"]
        self.assertEqual(command[command.index("--tax_scope") + 1], "auto")
        self.assertEqual(command[command.index("--dmnd_sensmode") + 1], "sensitive")
        self.assertEqual(command[command.index("--dmnd_iterate") + 1], "yes")
        self.assertIn("--dmnd_block_size", command)
        self.assertIn("--dmnd_index_chunks", command)

    def test_eggnog_sensitivity_change_invalidates_native_search(self):
        previous_entry = self.entry("eggnog")
        previous_command = previous_entry["search_method"]["command"]["native_stages"][
            "search"
        ]
        previous_command[previous_command.index("--dmnd_sensmode") + 1] = (
            "ultra-sensitive"
        )
        previous_entry["method_id"] = identity(
            {
                "search": previous_entry["search_method"],
                "interpretation": previous_entry["interpretation"],
            }
        )
        previous = self.result("eggnog")
        previous.update(task_identity(self.metadata, "eggnog", previous_entry))
        current_entry = self.entry("eggnog")
        planned = task_identity(self.metadata, "eggnog", current_entry)
        self.assertNotEqual(
            planned["search_fingerprint"], previous["search_fingerprint"]
        )
        self.assertEqual(
            current_entry["interpretation"], previous_entry["interpretation"]
        )
        with self.assertRaisesRegex(AnnotationError, "shared batch planner"):
            plan_task(
                self.bundle,
                "eggnog",
                current_entry,
                self.root / "sensitive",
                bundle_data=(self.metadata, self.proteins),
            )

    def test_kofam_uses_task_temporary_storage_with_invalid_inherited_directory(self):
        command = commands("kofam", 2, 32)
        command["version_commands"] = [[sys.executable, "--version"]]
        command["steps"] = [
            [
                sys.executable,
                "-c",
                """
import tempfile
from pathlib import Path
with tempfile.NamedTemporaryFile() as temporary:
    assert Path(temporary.name).parent == (Path.cwd() / 'scratch').resolve()
""",
            ]
        ]
        work = self.root / "native-temp"
        work.mkdir()
        script = work / "run.sh"
        script.write_text(shell_script(command))
        subprocess.run(
            ["bash", str(script)],
            cwd=work,
            check=True,
            env={**os.environ, "TMPDIR": str(work / "missing")},
        )
        self.assertEqual((work / "raw/exit_code.txt").read_text().strip(), "0")

    def test_search_reuse_interpretation_change_and_tampering(self):
        record = self.result()
        old = self.write_result(record, "old")
        planned = plan_task(
            self.bundle,
            "kofam",
            self.entry(),
            self.root / "reuse",
            old,
            bundle_data=(self.metadata, self.proteins),
        )
        self.assertEqual(planned["action"], "reuse")
        self.assertFalse((self.root / "reuse/input.faa").exists())
        self.assertFalse((self.root / "reuse/run.sh").exists())
        changed = self.entry()
        changed["interpretation"]["policy"] = {"test_policy": "changed"}
        planned = plan_task(
            self.bundle,
            "kofam",
            changed,
            self.root / "renormalize",
            old,
            bundle_data=(self.metadata, self.proteins),
        )
        self.assertEqual(planned["action"], "renormalize")
        self.assertFalse((self.root / "renormalize/input.faa").exists())
        self.assertFalse((self.root / "renormalize/run.sh").exists())
        changed = self.entry()
        changed["search_method"]["resource_id"] = "new_resource"
        planned = plan_task(
            self.bundle,
            "kofam",
            changed,
            self.root / "rerun",
            old,
            bundle_data=(self.metadata, self.proteins),
        )
        self.assertEqual(planned["action"], "run")
        (old / "raw" / "native.txt").write_text("changed native evidence\n")
        with self.assertRaisesRegex(AnnotationError, "changed"):
            validate_result(old)

    def test_unverified_execution_wrapper_identity_requires_a_new_search(self):
        previous = self.write_result(self.result(), "unverified-wrapper")
        current = self.entry()
        current["search_method"]["native_code"]["annotation_commands.py"] = digest(
            Path(__file__).resolve().parents[1] / "bin/annotation_commands.py"
        )
        current["method_id"] = identity(
            dict(
                search=current["search_method"],
                interpretation=current["interpretation"],
            )
        )
        planned = plan_task(
            self.bundle,
            "kofam",
            current,
            self.root / "verified-wrapper",
            previous,
            bundle_data=(self.metadata, self.proteins),
        )
        self.assertEqual(planned["action"], "run")

    def test_planner_validates_each_bundle_once_across_enabled_tools(self):
        samples = self.root / "samples.tsv"
        write_tsv(samples, ("accession",), [{"accession": "A"}])
        enabled = ["cogclassifier", "pfam", "kofam"]
        receipt = self.root / "preflight.json"
        entries = {}
        for tool in enabled:
            entry = self.entry(tool)
            entry["resource_path"] = str(self.root / "resource")
            entries[tool] = entry
        write_json(
            receipt,
            {"enabled_tools": enabled, "preflight_id": "fixture", "tools": entries},
        )
        with (
            patch("prepare_annotation_tasks.validate_preflight"),
            patch(
                "prepare_annotation_tasks.bundle_proteins", wraps=bundle_proteins
            ) as validate,
        ):
            planned = plan(samples, [self.bundle], receipt, self.root / "planned")
        self.assertEqual(validate.call_count, 1)
        self.assertEqual(sum(row["action"] == "run" for row in planned["tasks"]), 3)
        (self.bundle / "proteins.faa").write_text(">changed\nM\n")
        with (
            patch("prepare_annotation_tasks.validate_preflight"),
            self.assertRaises(AnnotationError),
        ):
            plan(samples, [self.bundle], receipt, self.root / "tampered-plan")

    def test_field_errors_invalidate_every_dependent_count(self):
        gene = self.proteins[0]["gene_id"]
        for field in FIELDS["eggnog"]:
            with self.subTest(field=field):
                features = {name: {gene: ["feature"]} for name in FIELDS["eggnog"]}
                error = dict(
                    gene_id=gene, field=field, raw_value="bad", reason="invalid"
                )
                record = self.result("eggnog", features, [error])
                self.assertIsNone(feature_counts(record, field))
                self.assertIsNone(accepted_count(record))
                other = "go" if field != "go" else "ko"
                self.assertEqual(feature_counts(record, other), {"feature": 1})
                group = {
                    tool: dict(tool=tool, status="skipped_disabled") for tool in TOOLS
                }
                group["eggnog"] = record
                summary = genome_summary(group, 1, True)
                self.assertIsNone(summary["eggnog_accepted_proteins"])
                self.assertIsNone(summary["eggnog_accepted_fraction"])
                if ("eggnog", field) in MATRICES.values():
                    self.assertIsNone(summary[f"eggnog_{field}_proteins"])
                    self.assertIsNone(summary[f"eggnog_{field}_assignments"])
                self.assertEqual(summary["eggnog_mapped_proteins"], 1)

    def test_category_fraction_and_system_totals_are_not_member_counts(self):
        gene = self.proteins[0]["gene_id"]
        group = {tool: dict(tool=tool, status="skipped_disabled") for tool in TOOLS}
        group["cogclassifier"] = self.result(
            "cogclassifier",
            {"cog": {gene: ["COG0001"]}, "categories": {gene: ["R", "S"]}},
        )
        group["padloc"] = self.result(
            "padloc", {"systems": {gene: ["A::padloc::1"], "second": ["A::padloc::1"]}}
        )
        summary = genome_summary(group, 2, True)
        self.assertEqual(summary["cog_category_R_gene_count"], 1)
        self.assertEqual(summary["cog_category_S_fractional_gene_count"], "0.5")
        self.assertEqual(summary["padloc_systems"], 1)
        group["cogclassifier"]["evidence"]["field_errors"] = [dict(field="categories")]
        summary = genome_summary(group, 2, True)
        self.assertIsNone(summary["cog_category_R_gene_count"])
        self.assertIsNone(summary["cog_category_S_fractional_gene_count"])
        self.assertEqual(summary["cogclassifier_cog_proteins"], 1)
        self.assertEqual(summary["cogclassifier_accepted_proteins"], 1)

    def test_aggregation_keeps_failed_genomes_and_zero_hit_axes(self):
        record = self.result()
        result_dir = self.write_result(record, "result")
        master, status, plan_path = [
            self.root / name for name in ("master.tsv", "status.tsv", "plan.json")
        ]
        write_tsv(
            master,
            ("Accession", "Gcode"),
            [dict(Accession="B", Gcode="11"), dict(Accession="A", Gcode="4")],
        )
        write_tsv(
            status,
            ("accession", "gcode_status"),
            [
                dict(accession="A", gcode_status="done"),
                dict(accession="B", gcode_status="failed"),
            ],
        )
        tasks = []
        for accession in ("A", "B"):
            for tool in TOOLS:
                selected = tool == "kofam"
                tasks.append(
                    dict(
                        accession=accession,
                        tool=tool,
                        action="run" if selected and accession == "A" else "skip",
                        status="planned"
                        if selected and accession == "A"
                        else "upstream_failed"
                        if selected
                        else "skipped_disabled",
                        reason="fixture",
                        input_proteins=1 if accession == "A" else None,
                    )
                )
        plan = dict(
            schema_version=1,
            accessions=["A", "B"],
            enabled_tools=["kofam"],
            tasks=tasks,
        )
        for row in tasks:
            if row["action"] == "run":
                row.update(
                    {
                        name: record[name]
                        for name in (
                            "search_fingerprint",
                            "normalization_fingerprint",
                            "method_id",
                        )
                    }
                )
        plan["plan_id"] = identity(plan)
        write_json(plan_path, plan)
        output = self.root / "report"
        manifest = aggregate(
            plan_path, master, status, [self.bundle], [result_dir], output
        )
        self.assertFalse(manifest["complete"])
        matrix = read_tsv(output / "tables/functional_matrices/kofam_ko_counts.tsv")
        self.assertEqual(
            matrix,
            [
                dict(accession="B", K00001="NA", K00002="NA"),
                dict(accession="A", K00001="1", K00002="1"),
            ],
        )
        masters = read_tsv(output / "tables/master_table.tsv")
        self.assertEqual(masters[1]["kofam_ko_assignments"], "2")
        self.assertEqual(masters[1]["kofam_ko_proteins"], "1")
        self.assertEqual(
            len(read_tsv(output / "tables/protein_function_summary.tsv")), 1
        )
        self.assertEqual(len(read_tsv(output / "tables/annotation_status.tsv")), 10)
        self.assertEqual(
            (output / "tables/functional_matrices/pfam_gene_counts.tsv").read_text(),
            "accession\nB\nA\n",
        )
        # A missing result for a planned search is a failure, even with valid input.
        output = self.root / "missing"
        aggregate(plan_path, master, status, [self.bundle], [], output)
        rows = read_tsv(output / "tables/annotation_status.tsv")
        missing = next(
            row for row in rows if row["accession"] == "A" and row["tool"] == "kofam"
        )
        self.assertEqual(missing["status"], "failed")
        self.assertEqual(missing["reason"], "missing_planned_result")
        # Even one successful result must belong to the requested method.
        record["method_id"] = "different_analysis"
        record.pop("result_id")
        record["result_id"] = identity(record)
        write_json(result_dir / "result.json", record)
        with self.assertRaisesRegex(AnnotationError, "planned analysis"):
            aggregate(
                plan_path,
                master,
                status,
                [self.bundle],
                [result_dir],
                self.root / "stale",
            )


if __name__ == "__main__":
    unittest.main()
