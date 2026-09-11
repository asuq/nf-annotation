#!/usr/bin/env python3
"""Compare the export with pinned native engine and writer code in the image."""

import importlib.util
import io
import multiprocessing
import random
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from eggnogmapper.annotation import output
from eggnogmapper.annotator.e7 import annotate
from eggnogmapper.go_namespace_export import FIELDS, KEY, write_row


def load_original(name, filename):
    spec = importlib.util.spec_from_file_location(
        name, Path("/opt/nf-annotation/provenance") / filename
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


ORIGINAL = load_original("eggnogmapper.annotator.e7._nf_original", "annotate.py")
ORIGINAL_OUTPUT = load_original("eggnogmapper.annotation._nf_original", "output.py")
NAMESPACES = {
    "GO:0003674": "molecular_function",
    "GO:0003700": "molecular_function",
    "GO:0008150": "biological_process",
    "GO:0006355": "biological_process",
    "GO:0005575": "cellular_component",
    "GO:0005634": "cellular_component",
}


class MemoryDB:
    def __init__(self, rows):
        self.rows = rows

    def get_protein_annotations_bulk(self, ids):
        return {key: dict(self.rows[key]) for key in ids if key in self.rows}

    def get_protein_name(self, name):
        return f"1.protein_{name}"


def engine(module, rows):
    value = module.AnnotationEngine.__new__(module.AnnotationEngine)
    value.db = MemoryDB(rows)
    value.donor_pool = "closest"
    value.lazy_cascade = True
    value._go_namespace_map = NAMESPACES
    return value


def run(module, mode, rows, meta):
    value = engine(module, rows)
    seed = 1
    parsed = value._pre_parse_batch(rows)
    if seed in parsed:
        parsed[seed].pop("pname", None)
    if mode == "eager":
        return value._summarize_annotations(
            rows, meta, "all", parsed=parsed, donor_pool="closest"
        )
    if mode == "masked":
        value.field_presence = value.build_field_presence(
            iter(rows.values()), max(rows) + 1
        )
        method = value._masked_cascade_summarize_batch
    else:
        method = value._lazy_cascade_summarize_batch
    seed_rows = value.db.get_protein_annotations_bulk([seed])
    seed_parsed = value._pre_parse_batch(seed_rows)
    if seed in seed_parsed:
        seed_parsed[seed].pop("pname", None)
    return method([seed], {seed: meta}, "all", seed_rows, seed_parsed, {seed})[seed]


def metadata(kind, depth):
    return {
        "event_id": -1,
        "ev_lca": "",
        "og_lca": "",
        "type": kind,
        "type_tier": annotate.AnnotationEngine.TYPE_TIERS[kind],
        "depth": depth,
        "in_seed_lineage": True,
    }


def tier_case():
    rows = {
        1: {"id": 1, "gos": "GO:0003674", "kegg_ko": "K00001"},
        2: {"id": 2, "gos": "GO:0008150", "kegg_ec": "1.1.1.1"},
        3: {"id": 3, "gos": "GO:0005575", "pfam": "PF00001"},
    }
    meta = {
        1: metadata("self", 40),
        2: metadata("one2many", 30),
        3: metadata("many2many", 20),
    }
    return rows, meta


def annotation_tuple(annotations, confidence):
    return (
        "pquery",
        "1",
        "1e-20",
        "100",
        annotations,
        ("OG", "S", "description"),
        "2",
        ["OG@2|Bacteria"],
        {},
        [],
        confidence,
        "Bacteria",
        "2",
        "Bacteria",
    )


def worker(mode):
    annotations, confidence = run(annotate, mode, *tier_case())
    native = io.StringIO()
    output.output_annotations_row(
        native, annotation_tuple(annotations, confidence), False, {}, MemoryDB({})
    )
    return native.getvalue(), annotations.get(KEY, {}), confidence


class ExportTests(unittest.TestCase):
    def check_equivalence(self, rows, meta):
        expected = run(ORIGINAL, "eager", rows, meta)
        namespace_evidence = None
        for mode in ("eager", "lazy", "masked"):
            # The unmodified three paths must also agree on this control.
            self.assertEqual(run(ORIGINAL, mode, rows, meta), expected)
            actual, confidence = run(annotate, mode, rows, meta)
            sidecar = actual.pop(KEY, {})
            self.assertEqual((actual, confidence), expected)
            if namespace_evidence is None:
                namespace_evidence = sidecar
            self.assertEqual(sidecar, namespace_evidence)
            original_row, patched_row = io.StringIO(), io.StringIO()
            ORIGINAL_OUTPUT.output_annotations_row(
                original_row, annotation_tuple(*expected), False, {}, MemoryDB(rows)
            )
            annotated = {**actual, KEY: sidecar}
            output.output_annotations_row(
                patched_row,
                annotation_tuple(annotated, confidence),
                False,
                {},
                MemoryDB(rows),
            )
            self.assertEqual(original_row.getvalue(), patched_row.getvalue())
            write_row(io.StringIO(), annotation_tuple(annotated, confidence))
        # Run each namespace through the unmodified native selector independently.
        for field in FIELDS:
            value = engine(ORIGINAL, rows)
            value.ANNOTATION_FIELDS = [field]
            assignments, confidence = value._summarize_annotations(
                rows, meta, "all", donor_pool="closest"
            )
            if assignments.get("GOs"):
                self.assertEqual(
                    namespace_evidence[field],
                    {
                        "terms": sorted(assignments["GOs"]),
                        "confidence": confidence["GOs"],
                    },
                )
            else:
                self.assertNotIn(field, namespace_evidence)

    def test_three_confidence_tiers(self):
        rows, meta = tier_case()
        self.check_equivalence(rows, meta)
        assignments, confidence = run(annotate, "eager", rows, meta)
        self.assertEqual(confidence["GOs"], "high")
        self.assertEqual(
            [assignments[KEY][field]["confidence"] for field in FIELDS],
            ["high", "medium", "low"],
        )

    def test_missing_namespace_and_no_go(self):
        rows, meta = tier_case()
        del rows[2]["gos"]
        self.check_equivalence(rows, meta)
        for row in rows.values():
            row.pop("gos", None)
        self.check_equivalence(rows, meta)

    def test_seeded_controls(self):
        rng = random.Random(20260911)
        terms = list(NAMESPACES)
        for _ in range(300):
            rows, meta = {}, {}
            for key in range(1, rng.randint(3, 45)):
                rows[key] = {
                    "id": key,
                    "gos": ",".join(rng.sample(terms, rng.randrange(7))),
                }
                if rng.random() < 0.5:
                    rows[key]["kegg_ko"] = f"K{rng.randrange(1, 20):05d}"
                kind = (
                    "self"
                    if key == 1
                    else rng.choice(["one2one", "one2many", "many2one", "many2many"])
                )
                meta[key] = metadata(kind, rng.randrange(1, 80))
            self.check_equivalence(rows, meta)

    def test_worker_serialization(self):
        modes = ["eager", "lazy", "masked"]
        with multiprocessing.get_context("spawn").Pool(2) as pool:
            results = pool.map(worker, modes)
        self.assertEqual(results, [worker(mode) for mode in modes])

    def test_parent_writer_and_native_bytes(self):
        rows, meta = tier_case()
        annotations = annotation_tuple(*run(annotate, "masked", rows, meta))
        stream = [((None, annotations), False), ((None, None), False)]
        with tempfile.TemporaryDirectory() as temporary:
            actual, expected = [
                str(Path(temporary) / name) for name in ("actual", "expected")
            ]
            with patch.object(output, "get_eggnog_db", return_value=MemoryDB(rows)):
                self.assertEqual(
                    list(
                        output.output_annotations(
                            iter(stream), actual, False, True, False, {}
                        )
                    ),
                    stream,
                )
            with patch.object(
                ORIGINAL_OUTPUT, "get_eggnog_db", return_value=MemoryDB(rows)
            ):
                list(
                    ORIGINAL_OUTPUT.output_annotations(
                        iter(stream), expected, False, True, False, {}
                    )
                )
            self.assertEqual(Path(actual).read_bytes(), Path(expected).read_bytes())
            exported = Path(actual + ".go_namespaces.tsv").read_text().splitlines()
            self.assertEqual(len(exported), 2)
            self.assertEqual(
                exported[1],
                "pquery\tGO:0003674\thigh\tGO:0008150\tmedium\tGO:0005575\tlow",
            )

    def test_missing_export_fails(self):
        annotations, confidence = run(annotate, "eager", *tier_case())
        del annotations[KEY]
        with self.assertRaisesRegex(ValueError, "disagrees"):
            write_row(io.StringIO(), annotation_tuple(annotations, confidence))

    def test_zero_hits_keeps_sidecar_header(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = str(Path(temporary) / "empty")
            with patch.object(output, "get_eggnog_db", return_value=MemoryDB({})):
                list(output.output_annotations(iter([]), path, False, True, False, {}))
            self.assertEqual(
                len(Path(path + ".go_namespaces.tsv").read_text().splitlines()), 1
            )

    def test_tool_resume_rejected(self):
        with self.assertRaisesRegex(ValueError, "fresh tool run"):
            list(output.output_annotations(iter([]), "unused", True, False, False, {}))


if __name__ == "__main__":
    unittest.main(verbosity=2)
