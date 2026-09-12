"""Exact ordering, failure semantics and bounded evidence memory during aggregation."""

from __future__ import annotations

import gc
import json
import sys
import tracemalloc
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import aggregate_annotations
import test_annotation_tasks as task_fixture
from annotation_common import (
    COORDINATE_COLUMNS,
    PROTEIN_COLUMNS,
    TOOLS,
    AnnotationError,
    bundle_proteins,
    digest,
    identity,
    read_json,
    read_tsv,
    write_json,
    write_tsv,
)
from annotation_result import inventory
from annotation_summary import FIELDS


class AnnotationAggregationScalingTests(unittest.TestCase):
    def setUp(self):
        self.fixture = task_fixture.AnnotationTaskTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root

    def cohort(self, accessions, tools=("kofam",), error_size=0):
        """Create valid one-protein bundles and explicit synthetic tool evidence."""
        bundles, results, tasks = [], [], []
        for accession in accessions:
            bundle = self.fixture.fixture.build(accession, name=f"bundle-{accession}")
            bundles.append(bundle)
            self.fixture.metadata, self.fixture.proteins = bundle_proteins(bundle)
            gene = self.fixture.proteins[0]["gene_id"]
            for tool in TOOLS:
                intended = dict(
                    accession=accession,
                    tool=tool,
                    action="run" if tool in tools else "skip",
                    status="planned" if tool in tools else "skipped_disabled",
                    reason="fixture",
                    input_proteins=1,
                )
                if tool in tools:
                    field = FIELDS[tool][0]
                    errors = (
                        [
                            dict(
                                gene_id=gene,
                                field=field,
                                raw_value="x" * error_size,
                                reason="invalid",
                            )
                        ]
                        if error_size
                        else []
                    )
                    record = self.fixture.result(
                        tool, {field: {gene: ["F1", "F1", "F2"]}}, errors
                    )
                    directory = self.fixture.write_result(
                        record, f"result-{accession}-{tool}"
                    )
                    # Multiple tools deliberately share this table and its sort key.
                    write_tsv(
                        directory / "normalized/evidence.tsv",
                        ("accession", "gene_id", "marker"),
                        [
                            dict(
                                accession=accession,
                                gene_id=gene,
                                marker=f"{tool}:second",
                            ),
                            dict(
                                accession=accession,
                                gene_id=gene,
                                marker=f"{tool}:first",
                            ),
                        ],
                    )
                    self.reidentify(directory)
                    results.append(directory)
                    intended.update(
                        {
                            name: record[name]
                            for name in (
                                "search_fingerprint",
                                "normalization_fingerprint",
                                "method_id",
                            )
                        }
                    )
                tasks.append(intended)
        master, status, plan_path = [
            self.root / name for name in ("master.tsv", "status.tsv", "plan.json")
        ]
        write_tsv(master, ("Accession",), [dict(Accession=acc) for acc in accessions])
        write_tsv(
            status,
            ("accession",),
            [dict(accession=acc) for acc in reversed(accessions)],
        )
        plan = dict(
            schema_version=1,
            accessions=accessions,
            enabled_tools=list(tools),
            tasks=tasks,
        )
        plan["plan_id"] = identity(plan)
        write_json(plan_path, plan)
        return plan_path, master, status, bundles, results

    def reidentify(self, directory, **changes):
        record = read_json(directory / "result.json")
        record.update(changes)
        record["normalized_files"] = inventory(directory / "normalized")
        record.pop("result_id")
        record["result_id"] = identity(record)
        write_json(directory / "result.json", record)

    def reidentify_bundle(self, bundle):
        metadata = read_json(bundle / "bundle.json")
        rows = read_tsv(bundle / "protein_manifest.tsv")
        metadata["files"] = {name: digest(bundle / name) for name in metadata["files"]}
        metadata["coordinate_id"] = metadata["files"]["gene_coordinates.tsv"]
        metadata["input_id"] = identity(
            dict(
                accession=metadata["accession"],
                genetic_code=metadata["genetic_code"],
                proteins=[
                    {name: row[name] for name in ("gene_id", "tool_id", "sha256")}
                    for row in rows
                ],
            )
        )
        write_json(bundle / "bundle.json", metadata)

    def test_stable_source_order_and_single_bundle_validation(self):
        args = list(self.cohort(["C", "A", "B"], ("kofam", "pfam")))
        args[3].reverse()
        args[4].reverse()
        output = self.root / "ordered"
        with patch.object(
            aggregate_annotations, "bundle_proteins", wraps=bundle_proteins
        ) as validate:
            manifest = aggregate_annotations.aggregate(*args, output)
        self.assertEqual(validate.call_count, 3)
        rows = read_tsv(output / "tables/evidence.tsv")
        self.assertEqual(
            [(row["accession"], row["marker"]) for row in rows],
            [
                (acc, marker)
                for acc in ("C", "A", "B")
                for marker in (
                    "kofam:second",
                    "kofam:first",
                    "pfam:second",
                    "pfam:first",
                )
            ],
        )
        self.assertEqual(
            (output / "tables/functional_matrices/kofam_ko_counts.tsv").read_text(),
            "accession\tF1\tF2\nC\t1\t1\nA\t1\t1\nB\t1\t1\n",
        )
        self.assertEqual(
            [
                row["accession"]
                for row in read_tsv(output / "tables/protein_function_summary.tsv")
            ],
            ["C", "A", "B"],
        )
        self.assertEqual(
            manifest["tables"],
            {
                str(path.relative_to(output)): digest(path)
                for path in (output / "tables").rglob("*.tsv")
            },
        )
        self.assertFalse(list(self.root.glob(".annotation-aggregation-*")))

    def test_coordinate_numeric_order_and_stable_ties(self):
        args = self.cohort(["A"], tools=())
        bundle = args[3][0]
        row = read_tsv(bundle / "gene_coordinates.tsv")[0]
        write_tsv(
            bundle / "gene_coordinates.tsv",
            COORDINATE_COLUMNS,
            [
                dict(row, segment=segment, attributes=label)
                for segment, label in (
                    (10, "last"),
                    (2, "first tie"),
                    (2, "second tie"),
                )
            ],
        )
        self.reidentify_bundle(bundle)
        output = self.root / "coordinates"
        aggregate_annotations.aggregate(*args, output)
        self.assertEqual(
            [
                row["attributes"]
                for row in read_tsv(output / "tables/gene_coordinates.tsv")
            ],
            ["first tie", "second tie", "last"],
        )

    def test_global_gene_collision_rejects_the_cohort_before_publication(self):
        args = self.cohort(["A", "B"], tools=())
        bundle = args[3][1]
        rows = read_tsv(bundle / "protein_manifest.tsv")
        rows[0]["gene_id"] = "A::gene_1"
        write_tsv(bundle / "protein_manifest.tsv", PROTEIN_COLUMNS, rows)
        self.reidentify_bundle(bundle)
        output = self.root / "collision"
        with self.assertRaisesRegex(AnnotationError, "collide across genomes"):
            aggregate_annotations.aggregate(*args, output)
        self.assertFalse(output.exists())
        self.assertFalse(list(self.root.glob(".annotation-aggregation-*")))

    def test_conflicting_definition_and_mixed_schema_fail_closed(self):
        args = self.cohort(["A", "B"])
        for directory, definition in zip(args[4], ("first", "second")):
            record = read_json(directory / "result.json")
            record["evidence"]["definitions"] = {"ko": {"F1": definition}}
            self.reidentify(directory, evidence=record["evidence"])
        with self.assertRaisesRegex(AnnotationError, "Conflicting definition"):
            aggregate_annotations.aggregate(*args, self.root / "definitions")
        directory = args[4][1]
        write_tsv(
            directory / "normalized/evidence.tsv",
            ("accession", "gene_id", "different"),
            [],
        )
        self.reidentify(directory)
        with self.assertRaisesRegex(AnnotationError, "Mixed normalized table schemas"):
            aggregate_annotations.aggregate(*args, self.root / "schemas")
        self.assertFalse((self.root / "definitions").exists())
        self.assertFalse((self.root / "schemas").exists())

    def test_result_integrity_gene_membership_and_shared_method_are_required(self):
        args = self.cohort(["A", "B"])
        directory = args[4][1]
        record = read_json(directory / "result.json")
        original_evidence = record["evidence"]
        self.reidentify(
            directory,
            evidence={**original_evidence, "mapped_genes": ["A::gene_1"]},
        )
        with self.assertRaisesRegex(AnnotationError, "undeclared protein"):
            aggregate_annotations.aggregate(*args, self.root / "unknown")
        self.reidentify(directory, evidence=original_evidence, method_id="other-method")
        plan = read_json(args[0])
        for task in plan["tasks"]:
            if task["accession"] == "B" and task["tool"] == "kofam":
                task["method_id"] = "other-method"
        plan.pop("plan_id")
        plan["plan_id"] = identity(plan)
        write_json(args[0], plan)
        with self.assertRaisesRegex(AnnotationError, "Mixed annotation methods"):
            aggregate_annotations.aggregate(*args, self.root / "methods")
        (directory / "raw/native.txt").write_text("changed evidence\n")
        with self.assertRaisesRegex(
            AnnotationError, "Native annotation evidence changed"
        ):
            aggregate_annotations.aggregate(*args, self.root / "raw")
        self.assertFalse(list(self.root.glob(".annotation-aggregation-*")))

    def test_evidence_memory_does_not_grow_with_cohort_rows(self):
        # The 250 kB malformed field mirrors a real large-field failure case; it
        # must be preserved verbatim while only one sample's evidence is retained.
        args = self.cohort([f"S{index:02}" for index in range(24)], error_size=250_000)
        gc.collect()
        tracemalloc.start()
        try:
            aggregate_annotations.aggregate(*args, self.root / "many")
            _, many_peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        # Reuse exactly the same source sample, changing only the cohort metadata.
        plan = read_json(args[0])
        plan["accessions"] = ["S00"]
        plan["tasks"] = [row for row in plan["tasks"] if row["accession"] == "S00"]
        plan.pop("plan_id")
        plan["plan_id"] = identity(plan)
        write_json(args[0], plan)
        write_tsv(args[1], ("Accession",), [dict(Accession="S00")])
        write_tsv(args[2], ("accession",), [dict(accession="S00")])
        gc.collect()
        tracemalloc.start()
        try:
            aggregate_annotations.aggregate(
                *args[:3], args[3][:1], args[4][:1], self.root / "one"
            )
            _, one_peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        self.assertLess(many_peak, one_peak + 2_000_000, (one_peak, many_peak))
        matrix = self.root / "many/tables/functional_matrices/kofam_ko_counts.tsv"
        self.assertEqual(
            matrix.read_text(),
            "accession\tF1\tF2\n"
            + "".join(f"S{index:02}\tNA\tNA\n" for index in range(24)),
        )
        statuses = read_tsv(self.root / "many/tables/annotation_status.tsv")

        errors = json.loads(
            next(row["field_errors"] for row in statuses if row["tool"] == "kofam")
        )
        self.assertEqual(errors[0]["raw_value"], "x" * 250_000)
        # All per-sample data rows have the same bytes in the one- and many-sample report.
        for name in (
            "protein_manifest.tsv",
            "gene_coordinates.tsv",
            "protein_function_summary.tsv",
            "evidence.tsv",
        ):
            one = (
                (self.root / "one/tables" / name).read_bytes().splitlines(keepends=True)
            )
            many = (
                (self.root / "many/tables" / name)
                .read_bytes()
                .splitlines(keepends=True)
            )
            self.assertEqual(one, many[: len(one)], name)


if __name__ == "__main__":
    unittest.main()
