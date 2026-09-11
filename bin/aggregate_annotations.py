#!/usr/bin/env python3
"""Publish authoritative annotation status, evidence, master fields and matrices."""

from __future__ import annotations

import argparse
import csv
import importlib
import logging
import sys
from pathlib import Path
from typing import Any

from annotation_common import (
    COORDINATE_COLUMNS,
    PROTEIN_COLUMNS,
    SCHEMA_VERSION,
    STATUS_COLUMNS,
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
from annotation_normalization import json_cell
from annotation_summary import (
    ANNOTATION_COLUMNS,
    FIELDS,
    MATRICES,
    accepted_count,
    feature_counts,
    genome_summary,
)
from annotation_tasks import inventory

PROVENANCE_COLUMNS = (
    "accession",
    "tool",
    "status",
    "action",
    "search_fingerprint",
    "normalization_fingerprint",
    "method_id",
    "runtime_id",
    "resource_id",
    "parameters",
    "interpretation",
    "result_id",
)
PROTEIN_SUMMARY_COLUMNS = (
    "accession",
    "gene_id",
    "protein_id",
    "contig_id",
    "strand",
    "genetic_code",
    *[f"{tool}_{name}" for tool in TOOLS for name in ("status", *FIELDS[tool])],
)


def header(path: Path) -> list[str]:
    """Read an already validated TSV's declared order, including a zero-row table."""
    with path.open(newline="") as handle:
        return next(csv.reader(handle, delimiter="\t"))


def aggregate(
    plan_path: Path,
    master: Path,
    sample_status: Path,
    bundles: list[Path],
    result_dirs: list[Path],
    output: Path,
) -> dict[str, Any]:
    """Reconcile the complete planned grid before deriving any functional counts."""
    plan = read_json(plan_path)
    if plan.get("schema_version") != SCHEMA_VERSION or plan.get("plan_id") != identity(
        {key: value for key, value in plan.items() if key != "plan_id"}
    ):
        raise AnnotationError("Invalid annotation plan")
    from validate_inputs import detect_metadata_key_column

    master_key = detect_metadata_key_column(header(master))
    master_rows = read_tsv(master, [master_key])
    accessions = [row[master_key] for row in master_rows]
    if len(accessions) != len(set(accessions)) or set(accessions) != set(
        plan["accessions"]
    ):
        raise AnnotationError("Master table and annotation accession sets differ")
    status_rows = read_tsv(sample_status, ["accession"])
    status_index = {row["accession"]: row for row in status_rows}
    if len(status_index) != len(status_rows) or set(status_index) != set(accessions):
        raise AnnotationError("Sample status and annotation accession sets differ")
    planned = {(row["accession"], row["tool"]): row for row in plan["tasks"]}
    if len(planned) != len(plan["tasks"]) or set(planned) != {
        (acc, tool) for acc in accessions for tool in TOOLS
    }:
        raise AnnotationError(
            "Annotation plan does not contain the complete accession/tool grid"
        )
    bundle_index, proteins, coordinates, bundle_records = {}, [], [], {}
    for bundle in bundles:
        metadata = read_json(bundle / "bundle.json")
        accession = metadata.get("accession")
        if accession not in accessions or accession in bundle_index:
            raise AnnotationError("Duplicate or undeclared bundle during aggregation")
        bundle_index[accession] = bundle
        bundle_records[accession] = dict(
            path=f"samples/{accession}/annotation/bundle",
            manifest_sha256=digest(bundle / "bundle.json"),
        )
        if metadata["status"] == "success":
            _, rows = bundle_proteins(bundle)
            proteins.extend(rows)
            coordinates.extend(
                read_tsv(bundle / "gene_coordinates.tsv", COORDINATE_COLUMNS)
            )
    proteins_by_gene = {row["gene_id"]: row for row in proteins}
    if len(proteins_by_gene) != len(proteins):
        raise AnnotationError("Canonical protein IDs collide across genomes")
    records, details, result_manifest = {}, {}, {}
    for tool in TOOLS:
        for name, columns in importlib.import_module(
            "normalize_" + tool
        ).TABLE_COLUMNS.items():
            details[name] = (list(columns), [])
    for directory in result_dirs:
        record = read_json(directory / "result.json")
        key = record.get("accession"), record.get("tool")
        if key not in planned or key in records or planned[key]["action"] == "skip":
            raise AnnotationError("Unexpected, disabled or duplicate annotation result")
        if (
            record.get("schema_version") != SCHEMA_VERSION
            or record.get("result_id")
            != identity(
                {name: value for name, value in record.items() if name != "result_id"}
            )
            or record.get("status") not in ("success", "failed")
        ):
            raise AnnotationError("Invalid annotation task result")
        if any(
            record.get(name) != planned[key].get(name)
            for name in (
                "search_fingerprint",
                "normalization_fingerprint",
                "method_id",
            )
        ):
            raise AnnotationError("Annotation result differs from the planned analysis")
        metadata, rows = bundle_proteins(bundle_index[key[0]])
        if record.get("input_id") != metadata["input_id"] or record.get(
            "input_proteins"
        ) != len(rows):
            raise AnnotationError(
                "Result input identity differs from the canonical bundle"
            )
        if inventory(directory / "raw") != record["raw_files"]:
            raise AnnotationError("Raw evidence changed before aggregation")
        if record["status"] == "success":
            if inventory(directory / "normalized") != record["normalized_files"]:
                raise AnnotationError("Normalized evidence changed before aggregation")
            allowed_genes = {row["gene_id"] for row in rows}
            evidence = record["evidence"]
            returned_genes = set(evidence["mapped_genes"]) | set(
                evidence["accepted_genes"]
            )
            for field, genes in evidence["features"].items():
                if field not in FIELDS[key[1]]:
                    raise AnnotationError(f"Unknown normalized feature field: {field}")
                returned_genes.update(genes)
            returned_genes.update(
                error["gene_id"] for error in evidence["field_errors"]
            )
            if not returned_genes <= allowed_genes:
                raise AnnotationError(
                    "Normalized evidence contains an undeclared protein"
                )
            for path in sorted((directory / "normalized").glob("*.tsv")):
                columns = header(path)
                table = details.setdefault(path.name, (columns, []))
                if columns != table[0]:
                    raise AnnotationError(
                        f"Mixed normalized table schemas: {path.name}"
                    )
                table[1].extend(read_tsv(path))
        records[key] = record
        result_manifest.setdefault(key[0], {})[key[1]] = dict(
            path=f"samples/{key[0]}/annotation/{key[1]}", result_id=record["result_id"]
        )
    for key, intended in planned.items():
        if key not in records:
            records[key] = dict(
                accession=key[0],
                tool=key[1],
                status=intended["status"] if intended["action"] == "skip" else "failed",
                reason=intended["reason"]
                if intended["action"] == "skip"
                else "missing_planned_result",
                input_proteins=intended["input_proteins"],
                evidence=None,
            )
    # Comparing genomes requires a shared method; sample input hashes may differ.
    for tool in TOOLS:
        methods = {
            record["method_id"]
            for (acc, name), record in records.items()
            if name == tool and record["status"] == "success"
        }
        if len(methods) > 1:
            raise AnnotationError(
                f"Mixed annotation methods in the current cohort: {tool}"
            )
    output.mkdir()
    tables = output / "tables"
    tables.mkdir()
    matrix_root = tables / "functional_matrices"
    matrix_root.mkdir()
    ordered_records = [records[(acc, tool)] for acc in accessions for tool in TOOLS]
    statuses, provenance, failures = [], [], []
    for record in ordered_records:
        row = {name: record.get(name) for name in STATUS_COLUMNS}
        row.update(
            schema_version=SCHEMA_VERSION,
            accepted_proteins=accepted_count(record),
            field_errors=json_cell(
                (record.get("evidence") or {}).get("field_errors", [])
            ),
        )
        statuses.append(row)
        method = record.get("search", {}).get("method", {})
        provenance.append(
            dict(
                **{
                    name: record.get(name)
                    for name in PROVENANCE_COLUMNS
                    if name
                    not in ("runtime_id", "resource_id", "parameters", "interpretation")
                },
                runtime_id=method.get("runtime_id"),
                resource_id=method.get("resource_id"),
                parameters=json_cell(method.get("command")),
                interpretation=json_cell(record.get("interpretation")),
            )
        )
        if record["status"] in ("failed", "upstream_failed", "incompatible_input"):
            failures.append(
                dict(
                    accession=record["accession"],
                    tool=record["tool"],
                    reason=record["reason"],
                )
            )
    write_tsv(tables / "annotation_status.tsv", STATUS_COLUMNS, statuses)
    write_tsv(tables / "annotation_provenance.tsv", PROVENANCE_COLUMNS, provenance)
    for accession in accessions:
        group = {tool: records[(accession, tool)] for tool in TOOLS}
        count = sum(row["accession"] == accession for row in proteins) or None
        complete = all(
            result["status"] == "success"
            for tool, result in group.items()
            if tool in plan["enabled_tools"]
        )
        master_row = next(row for row in master_rows if row[master_key] == accession)
        master_row.update(genome_summary(group, count, complete))
        for tool, record in group.items():
            status_index[accession][f"{tool}_status"] = record["status"]
        status_index[accession]["annotation_complete"] = str(complete).lower()
    master_columns = header(master)
    master_columns.extend(
        name for name in ANNOTATION_COLUMNS if name not in master_columns
    )
    status_columns = header(sample_status)
    status_columns.extend(
        name
        for name in (*[f"{tool}_status" for tool in TOOLS], "annotation_complete")
        if name not in status_columns
    )
    write_tsv(tables / "master_table.tsv", master_columns, master_rows)
    write_tsv(
        tables / "sample_status.tsv",
        status_columns,
        [status_index[acc] for acc in accessions],
    )
    order = {acc: index for index, acc in enumerate(accessions)}
    proteins.sort(key=lambda row: (order[row["accession"]], row["tool_id"]))
    coordinates.sort(
        key=lambda row: (order[row["accession"]], row["gene_id"], int(row["segment"]))
    )
    write_tsv(tables / "protein_manifest.tsv", PROTEIN_COLUMNS, proteins)
    write_tsv(tables / "gene_coordinates.tsv", COORDINATE_COLUMNS, coordinates)
    protein_summaries = []
    for protein in proteins:
        row = {name: protein[name] for name in PROTEIN_SUMMARY_COLUMNS[:6]}
        for tool in TOOLS:
            record = records[(protein["accession"], tool)]
            row[f"{tool}_status"] = record["status"]
            evidence = record.get("evidence") or {}
            errors = {
                error["field"]
                for error in evidence.get("field_errors", [])
                if error["gene_id"] == protein["gene_id"]
            }
            for field in FIELDS[tool]:
                value = None
                if record["status"] == "success" and field not in errors:
                    value = json_cell(
                        evidence["features"].get(field, {}).get(protein["gene_id"], [])
                    )
                row[f"{tool}_{field}"] = value
        protein_summaries.append(row)
    write_tsv(
        tables / "protein_function_summary.tsv",
        PROTEIN_SUMMARY_COLUMNS,
        protein_summaries,
    )
    for name, (columns, rows) in details.items():
        rows.sort(key=lambda row: (order[row["accession"]], row.get("gene_id", "")))
        write_tsv(tables / name, columns, rows)
    catalogue = []
    for matrix, (tool, field) in MATRICES.items():
        features, definitions = set(), {}
        method_ids, resource_ids = set(), set()
        for accession in accessions:
            record = records[(accession, tool)]
            if record["status"] == "success":
                evidence = record["evidence"]
                features.update(
                    value
                    for values in evidence["features"].get(field, {}).values()
                    for value in values
                )
                for feature, definition in (
                    evidence["definitions"].get(field, {}).items()
                ):
                    if feature in definitions and definitions[feature] != definition:
                        raise AnnotationError(f"Conflicting definition for {feature}")
                    definitions[feature] = definition
                method_ids.add(record["method_id"])
                resource_ids.add(record["search"]["method"]["resource_id"])
        features = sorted(features)
        matrix_rows = []
        for accession in accessions:
            counts = feature_counts(records[(accession, tool)], field)
            matrix_rows.append(
                dict(
                    accession=accession,
                    **{
                        feature: counts[feature] if counts is not None else None
                        for feature in features
                    },
                )
            )
        write_tsv(
            matrix_root / (matrix + ".tsv"), ("accession", *features), matrix_rows
        )
        for feature in features:
            catalogue.append(
                dict(
                    matrix=matrix + ".tsv",
                    feature_id=feature,
                    definition=definitions.get(feature),
                    source=tool,
                    resource_id=next(iter(resource_ids)),
                    method_id=next(iter(method_ids)),
                )
            )
    write_tsv(
        matrix_root / "feature_catalogue.tsv",
        ("matrix", "feature_id", "definition", "source", "resource_id", "method_id"),
        catalogue,
    )
    write_json(
        output / "annotation_acceptance.json",
        dict(complete=not failures, failures=failures, plan_id=plan["plan_id"]),
    )
    manifest = dict(
        schema_version=SCHEMA_VERSION,
        pipeline_contract="nf-annotation-v0.4",
        accessions=accessions,
        enabled_tools=plan["enabled_tools"],
        complete=not failures,
        plan_id=plan["plan_id"],
        bundles=bundle_records,
        results=result_manifest,
        tables={
            str(path.relative_to(output)): digest(path)
            for path in sorted(tables.rglob("*.tsv"))
        },
    )
    write_json(output / "annotation_results.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("plan", "master", "sample-status", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--bundle", type=Path, action="append", default=[])
    parser.add_argument("--result", type=Path, action="append", default=[])
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        aggregate(
            args.plan,
            args.master,
            args.sample_status,
            args.bundle,
            args.result,
            args.output,
        )
    except (AnnotationError, OSError) as error:
        logging.error("%s", error)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
