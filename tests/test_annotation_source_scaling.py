"""Source-table projection preserves strict CSV validation and bounded memory."""

from __future__ import annotations

import csv
import gc
import shutil
import sys
import tracemalloc
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import test_annotation_aggregation_scaling as aggregation_fixture
import test_annotation_source as source_fixture
from aggregate_annotations import aggregate
from annotation_common import (
    AnnotationError,
    identity,
    read_json,
    read_tsv,
    write_json,
    write_tsv,
)
from annotation_source import source_rows, validate_source
from annotation_workflow_fixture import refresh_tables


class AnnotationSourceScalingTests(unittest.TestCase):
    def setUp(self):
        self.fixture = source_fixture.AnnotationSourceTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root, self.source = self.fixture.root, self.fixture.source

    def test_large_quoted_cells_match_shared_reader_and_restore_csv_limit(self):
        previous_limit = csv.field_size_limit()
        path = self.root / "quoted.tsv"
        expected = [
            dict(
                accession="A",
                tool="eggnog",
                status="success",
                field_errors='A "quoted" field\twith tabs\nand newlines\n' * 6000,
            )
        ]
        write_tsv(path, expected[0], expected)
        self.assertGreater(len(expected[0]["field_errors"]), previous_limit)
        with source_rows(path, ("accession", "tool", "status")) as (header, rows):
            self.assertEqual(header, list(expected[0]))
            self.assertEqual(list(rows), read_tsv(path))
        self.assertEqual(csv.field_size_limit(), previous_limit)
        for value in (
            "accession\ttool\tstatus\tnote\nA\teggnog\tsuccess\n",
            "accession\ttool\tstatus\nA\teggnog\tsuccess\textra\n",
            "accession\ttool\tstatus\tstatus\n",
            "accession\ttool\tstatus\t\n",
        ):
            path.write_text(value)
            with self.subTest(value=value), self.assertRaises(AnnotationError):
                with source_rows(path, ("accession", "tool", "status")) as (_, rows):
                    list(rows)
            self.assertEqual(csv.field_size_limit(), previous_limit)

    def test_duplicate_missing_foreign_and_malformed_status_grid_fail(self):
        path = self.source / "tables/annotation_status.tsv"
        original = read_tsv(path)
        columns = list(original[0])
        for name, rows, message in (
            ("duplicate", [*original, original[0]], "Duplicate published annotation"),
            ("missing", original[:-1], "complete declared grid"),
            (
                "foreign",
                [*original[:-1], dict(original[-1], accession="foreign")],
                "complete declared grid",
            ),
        ):
            with self.subTest(name=name):
                write_tsv(path, columns, rows)
                refresh_tables(self.source)
                with self.assertRaisesRegex(AnnotationError, message):
                    validate_source(self.source)
        write_tsv(path, columns, original)
        with path.open("a") as handle:
            handle.write("A\teggnog\n")
        refresh_tables(self.source)
        previous_limit = csv.field_size_limit()
        with self.assertRaisesRegex(AnnotationError, "Malformed table row"):
            validate_source(self.source)
        self.assertEqual(csv.field_size_limit(), previous_limit)

    def test_projection_still_requires_all_master_columns_and_complete_unused_cells(
        self,
    ):
        path = self.source / "tables/master_table.tsv"
        rows = read_tsv(path)
        columns = list(rows[0])
        omitted = "annotation_field_errors"
        write_tsv(
            path,
            [key for key in columns if key != omitted],
            [
                {key: value for key, value in row.items() if key != omitted}
                for row in rows
            ],
        )
        refresh_tables(self.source)
        with self.assertRaisesRegex(AnnotationError, "Invalid table header"):
            validate_source(self.source)
        write_tsv(path, [*columns, "unused_note"], rows)
        lines = path.read_text().splitlines()
        path.write_text(lines[0] + "\n" + lines[1].removesuffix("\t") + "\n")
        refresh_tables(self.source)
        with self.assertRaisesRegex(AnnotationError, "Malformed table row"):
            validate_source(self.source)

    def publish(self, args, name):
        output = args[0].parent / name
        manifest = aggregate(*args, output)
        for bundle in args[3]:
            accession = read_json(bundle / "bundle.json")["accession"]
            shutil.copytree(bundle, output / manifest["bundles"][accession]["path"])
        for directory in args[4]:
            record = read_json(directory / "result.json")
            shutil.copytree(
                directory,
                output
                / manifest["results"][record["accession"]][record["tool"]]["path"],
            )
        # Both optional metadata tables can carry large quoted notes as well.
        note = 'Quoted "native" metadata\twith a\nsecond line\n' * 6000
        for name in ("master_table.tsv", "sample_status.tsv"):
            path = output / "tables" / name
            rows = read_tsv(path)
            for row in rows:
                row["analysis_note"] = note
            write_tsv(path, rows[0], rows)
        refresh_tables(output)
        return output

    def peak_validation(self, source):
        expected = read_json(source / "annotation_results.json")
        gc.collect()
        previous_limit = csv.field_size_limit()
        tracemalloc.start()
        try:
            actual, batches = validate_source(source)
            _, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        self.assertEqual(actual, expected)
        self.assertEqual(batches, {})
        self.assertEqual(csv.field_size_limit(), previous_limit)
        return peak

    def test_source_validation_heap_does_not_retain_cohort_evidence_cells(self):
        fixture = aggregation_fixture.AnnotationAggregationScalingTests()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        args = fixture.cohort(
            [f"S{index:02}" for index in range(24)], error_size=250_000
        )
        master_rows = read_tsv(args[1])
        write_tsv(
            args[1], ("Accession", "Gcode"), [dict(row, Gcode=4) for row in master_rows]
        )
        many = self.publish(args, "many-source")
        plan = read_json(args[0])
        plan["accessions"] = ["S00"]
        plan["tasks"] = [row for row in plan["tasks"] if row["accession"] == "S00"]
        plan.pop("plan_id")
        plan["plan_id"] = identity(plan)
        write_json(args[0], plan)
        write_tsv(args[1], ("Accession", "Gcode"), [dict(Accession="S00", Gcode=4)])
        write_tsv(args[2], ("accession",), [dict(accession="S00")])
        one = self.publish((*args[:3], args[3][:1], args[4][:1]), "one-source")
        one_peak = self.peak_validation(one)
        many_peak = self.peak_validation(many)
        self.assertLess(many_peak, one_peak + 1_000_000, (one_peak, many_peak))


if __name__ == "__main__":
    unittest.main()
