#!/usr/bin/env python3
"""Publish authoritative annotation status, evidence, master fields and matrices."""

from __future__ import annotations

import argparse
import csv
import importlib
import json
import logging
import sys
import tempfile
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any, TypedDict

from annotation_common import (
    COORDINATE_COLUMNS,
    PROTEIN_COLUMNS,
    SCHEMA_VERSION,
    STATUS_COLUMNS,
    TOOLS,
    AnnotationError,
    bundle_proteins,
    canonical_gene_id,
    digest,
    identity,
    read_json,
    read_tsv,
    write_json,
    write_tsv,
)
from annotation_normalization import json_cell
from annotation_path_lists import read_path_list, read_publication_paths
from annotation_result import (
    NativeBatch,
    inventory,
    native_batch_index,
    validate_raw_evidence,
)
from annotation_summary import (
    ANNOTATION_COLUMNS,
    FIELDS,
    MATRICES,
    accepted_count,
    feature_counts,
    genome_summary,
)

LOGGER = logging.getLogger(__name__)

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


class MatrixState(TypedDict):
    features: set[str]
    definitions: dict[str, str]
    method_ids: set[str]
    resource_ids: set[str]


class TableWriter:
    """Stream sample-ordered tables with bounded sparse feature counts."""

    def __init__(self, path: Path):
        self.count_handles = {}
        self.directory = path
        self.directory.mkdir()
        self.pending = {}
        self.accession = None
        self.schemas: dict[str, Sequence[str]] = {
            "proteins": PROTEIN_COLUMNS,
            "coordinates": COORDINATE_COLUMNS,
            "summaries": PROTEIN_SUMMARY_COLUMNS,
            "status": STATUS_COLUMNS,
            "provenance": PROVENANCE_COLUMNS,
        }
        self.handles = {}
        self.writers = {}

    def add_genes(self, proteins):
        # A canonical key encodes accession and protein ID. Bundle-local ID
        # uniqueness and cohort accession uniqueness therefore imply global uniqueness.
        if any(
            row["gene_id"] != canonical_gene_id(row["accession"], row["protein_id"])
            for row in proteins
        ):
            raise AnnotationError(
                "Noncanonical protein IDs could collide across genomes"
            )

    def add_counts(self, matrix, index, counts):
        if matrix not in self.count_handles:
            self.count_handles[matrix] = (
                self.directory.parent / (matrix + ".jsonl")
            ).open("w+")
        self.count_handles[matrix].write(json.dumps([index, counts]) + "\n")

    def add_rows(
        self, name, rows, order, *, gene_key="", coordinate=False, source_order=0
    ):
        pending = self.pending.setdefault(name, [])
        for index, row in enumerate(rows):
            if row["accession"] != self.accession:
                raise AnnotationError("Report row belongs to a different accession")
            pending.append(
                (
                    (
                        row.get(gene_key, ""),
                        int(row["segment"]) if coordinate else 0,
                        source_order,
                        index,
                    ),
                    row,
                )
            )

    def finish_sample(self):
        for name, rows in self.pending.items():
            writer = self.writer(name)
            for _, row in sorted(rows, key=lambda item: item[0]):
                writer.writerow(
                    {
                        key: "NA" if value is None else value
                        for key, value in row.items()
                    }
                )
        self.pending.clear()

    def writer(self, name):
        if name not in self.writers:
            filename = {
                "proteins": "protein_manifest.tsv",
                "coordinates": "gene_coordinates.tsv",
                "summaries": "protein_function_summary.tsv",
                "status": "annotation_status.tsv",
                "provenance": "annotation_provenance.tsv",
            }.get(name, name.removeprefix("detail:"))
            handle = (self.directory / filename).open("w", newline="")
            writer = csv.DictWriter(
                handle,
                fieldnames=self.schemas[name],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            self.handles[name], self.writers[name] = handle, writer
        return self.writers[name]

    def close(self):
        for handle in self.handles.values():
            handle.close()
        for handle in self.count_handles.values():
            handle.close()

    def matrix_rows(self, matrix, features, accessions):
        handle = self.count_handles[matrix]
        handle.seek(0)
        for line in handle:
            index, counts = json.loads(line)
            yield dict(
                accession=accessions[index],
                **{
                    feature: counts.get(feature, 0) if counts is not None else None
                    for feature in features
                },
            )


def result_record(
    directory: Path,
    intended: dict[str, Any],
    metadata: dict[str, Any],
    proteins: list[dict[str, str]],
    batches: dict[str, NativeBatch],
) -> dict[str, Any]:
    """Validate one task against its plan, bundle and immutable evidence."""
    record = read_json(directory / "result.json")
    if (
        (record.get("accession"), record.get("tool"))
        != (intended["accession"], intended["tool"])
        or record.get("schema_version") != SCHEMA_VERSION
        or record.get("result_id")
        != identity(
            {name: value for name, value in record.items() if name != "result_id"}
        )
        or record.get("status") not in ("success", "failed")
    ):
        raise AnnotationError("Invalid annotation task result")
    if any(
        record.get(name) != intended.get(name)
        for name in ("search_fingerprint", "normalization_fingerprint", "method_id")
    ):
        raise AnnotationError("Annotation result differs from the planned analysis")
    if record.get("input_id") != metadata["input_id"] or record.get(
        "input_proteins"
    ) != len(proteins):
        raise AnnotationError("Result input identity differs from the canonical bundle")
    validate_raw_evidence(directory, record, batches)
    if record["status"] == "success":
        if inventory(directory / "normalized") != record["normalized_files"]:
            raise AnnotationError("Normalized evidence changed before aggregation")
        allowed_genes = {row["gene_id"] for row in proteins}
        evidence = record["evidence"]
        returned_genes = set(evidence["mapped_genes"]) | set(evidence["accepted_genes"])
        for field, genes in evidence["features"].items():
            if field not in FIELDS[record["tool"]]:
                raise AnnotationError(f"Unknown normalized feature field: {field}")
            returned_genes.update(genes)
        returned_genes.update(error["gene_id"] for error in evidence["field_errors"])
        if not returned_genes <= allowed_genes:
            raise AnnotationError("Normalized evidence contains an undeclared protein")
    return record


def protein_summaries(
    proteins: list[dict[str, str]], group: dict[str, dict[str, Any]]
) -> Iterator[dict[str, Any]]:
    """Derive rows with sample-local field-error indexes instead of repeated scans."""
    errors: dict[str, dict[str, set[str]]] = {}
    for tool, record in group.items():
        errors[tool] = {}
        for error in (record.get("evidence") or {}).get("field_errors", []):
            errors[tool].setdefault(error["gene_id"], set()).add(error["field"])
    for protein in proteins:
        row: dict[str, Any] = {
            name: protein[name] for name in PROTEIN_SUMMARY_COLUMNS[:6]
        }
        gene = protein["gene_id"]
        for tool in TOOLS:
            record = group[tool]
            row[f"{tool}_status"] = record["status"]
            evidence = record.get("evidence") or {}
            for field in FIELDS[tool]:
                value = None
                if record["status"] == "success" and field not in errors[tool].get(
                    gene, ()
                ):
                    value = json_cell(evidence["features"].get(field, {}).get(gene, []))
                row[f"{tool}_{field}"] = value
        yield row


def record_rows(
    group: dict[str, dict[str, Any]],
) -> Iterator[dict[str, Any]]:
    """Serialize status without retaining cohort field-error arrays."""
    for record in group.values():
        row = {name: record.get(name) for name in STATUS_COLUMNS}
        row.update(
            schema_version=SCHEMA_VERSION,
            accepted_proteins=accepted_count(record),
            field_errors=json_cell(
                (record.get("evidence") or {}).get("field_errors", [])
            ),
        )
        yield row


def provenance_rows(
    group: dict[str, dict[str, Any]],
) -> Iterator[dict[str, Any]]:
    """Preserve the published native-search and interpretation provenance."""
    for record in group.values():
        method = record.get("search", {}).get("method", {})
        yield dict(
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


def aggregate(
    plan_path: Path,
    master: Path,
    sample_status: Path,
    bundles: list[Path],
    result_dirs: list[dict[str, Any]],
    output: Path,
    *,
    batch_dirs: Sequence[Path] = (),
) -> dict[str, Any]:
    """Reconcile the planned grid using one sample plus compact cohort metadata."""
    with tempfile.TemporaryDirectory(
        prefix=".annotation-aggregation-", dir=output.parent
    ) as temporary:
        work = Path(temporary)
        spool = TableWriter(work / "tables")
        try:
            return aggregate_samples(
                plan_path,
                master,
                sample_status,
                bundles,
                result_dirs,
                output,
                work,
                spool,
                native_batch_index(list(batch_dirs)),
            )
        finally:
            spool.close()


def aggregate_samples(
    plan_path: Path,
    master: Path,
    sample_status: Path,
    bundles: list[Path],
    result_dirs: list[dict[str, Any]],
    output: Path,
    work: Path,
    spool: TableWriter,
    batches: dict[str, NativeBatch],
) -> dict[str, Any]:
    """Validate and spool every sample before publishing the complete report."""
    plan = read_json(plan_path)
    if plan.get("schema_version") != SCHEMA_VERSION or plan.get("plan_id") != identity(
        {key: value for key, value in plan.items() if key != "plan_id"}
    ):
        raise AnnotationError("Invalid annotation plan")
    from validate_inputs import detect_metadata_key_column

    master_key = detect_metadata_key_column(header(master))
    master_rows = read_tsv(master, [master_key])
    accessions = [row[master_key] for row in master_rows]
    order = {acc: index for index, acc in enumerate(accessions)}
    if len(accessions) != len(order) or set(accessions) != set(plan["accessions"]):
        raise AnnotationError("Master table and annotation accession sets differ")
    for batch in batches.values():
        if not batch.members.keys() <= order.keys():
            raise AnnotationError(
                "Native eggNOG batch contains undeclared cohort members"
            )
    batch_members: dict[str, set[str]] = {batch_id: set() for batch_id in batches}
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
    bundle_index, bundle_records = {}, {}
    for source_order, bundle in enumerate(bundles):
        metadata = read_json(bundle / "bundle.json")
        accession = metadata.get("accession")
        if accession not in order or accession in bundle_index:
            raise AnnotationError("Duplicate or undeclared bundle during aggregation")
        bundle_index[accession] = (source_order, bundle)
        bundle_records[accession] = {
            "path": f"samples/{accession}/annotation/bundle",
            "manifest_sha256": digest(bundle / "bundle.json"),
        }
    # Carry producer identities with paths; read full evidence only for its sample.
    if not isinstance(result_dirs, list):
        raise AnnotationError("Annotation results require a keyed path list")
    result_index = {}
    for source_order, entry in enumerate(result_dirs):
        if not isinstance(entry, dict) or set(entry) != {"accession", "tool", "path"}:
            raise AnnotationError("Annotation results require keyed paths")
        directory = Path(entry["path"])
        key = entry["accession"], entry["tool"]
        if (
            key not in planned
            or key in result_index
            or planned[key]["action"] == "skip"
        ):
            raise AnnotationError("Unexpected, disabled or duplicate annotation result")
        result_index[key] = (source_order, directory)
    details = {}
    for tool in TOOLS:
        for name, columns in importlib.import_module(
            "normalize_" + tool
        ).TABLE_COLUMNS.items():
            details[name] = list(columns)
    matrices: dict[str, MatrixState] = {
        name: {
            "features": set(),
            "definitions": {},
            "method_ids": set(),
            "resource_ids": set(),
        }
        for name in MATRICES
    }
    methods = {tool: set() for tool in TOOLS}
    failures, result_manifest = [], {}
    for master_row in master_rows:
        accession = master_row[master_key]
        spool.accession = accession
        proteins, metadata = [], None
        if accession in bundle_index:
            source_order, bundle = bundle_index[accession]
            metadata = read_json(bundle / "bundle.json")
            if metadata["status"] == "success":
                metadata, proteins = bundle_proteins(bundle)
                spool.add_genes(proteins)
                spool.add_rows("proteins", proteins, order, gene_key="tool_id")
                spool.add_rows(
                    "coordinates",
                    read_tsv(bundle / "gene_coordinates.tsv", COORDINATE_COLUMNS),
                    order,
                    gene_key="gene_id",
                    coordinate=True,
                    source_order=source_order,
                )
        group: dict[str, dict[str, Any]] = {}
        record: dict[str, Any]
        for tool in TOOLS:
            key = accession, tool
            intended = planned[key]
            if key in result_index:
                if metadata is None or metadata["status"] != "success":
                    raise AnnotationError(
                        "Annotation result has no successful canonical bundle"
                    )
                source_order, directory = result_index[key]
                record = result_record(directory, intended, metadata, proteins, batches)
                if "native_batch" in record:
                    batch_members[record["native_batch"]["batch_id"]].add(accession)
                result_manifest.setdefault(accession, {})[tool] = {
                    "path": f"samples/{accession}/annotation/{tool}",
                    "result_id": record["result_id"],
                }
                if record["status"] == "success":
                    methods[tool].add(record["method_id"])
                    if len(methods[tool]) > 1:
                        raise AnnotationError(
                            f"Mixed annotation methods in the current cohort: {tool}"
                        )
                    for path in sorted((directory / "normalized").glob("*.tsv")):
                        columns = header(path)
                        if columns != details.setdefault(path.name, columns):
                            raise AnnotationError(
                                f"Mixed normalized table schemas: {path.name}"
                            )
                        spool.schemas["detail:" + path.name] = columns
                        spool.add_rows(
                            "detail:" + path.name,
                            read_tsv(path),
                            order,
                            gene_key="gene_id",
                            source_order=source_order,
                        )
            else:
                record = {
                    "accession": accession,
                    "tool": tool,
                    "status": intended["status"]
                    if intended["action"] == "skip"
                    else "failed",
                    "reason": intended["reason"]
                    if intended["action"] == "skip"
                    else "missing_planned_result",
                    "input_proteins": intended["input_proteins"],
                    "evidence": None,
                }
            group[tool] = record
            if record["status"] in ("failed", "upstream_failed", "incompatible_input"):
                failures.append(
                    {"accession": accession, "tool": tool, "reason": record["reason"]}
                )
            status_index[accession][f"{tool}_status"] = record["status"]
        complete = all(
            record["status"] == "success"
            for tool, record in group.items()
            if tool in plan["enabled_tools"]
        )
        master_row.update(genome_summary(group, len(proteins) or None, complete))
        status_index[accession]["annotation_complete"] = str(complete).lower()
        spool.add_rows("status", record_rows(group), order)
        spool.add_rows("provenance", provenance_rows(group), order)
        # The spool supplies tool_id ordering for proteins and their summary rows.
        proteins.sort(key=lambda row: row["tool_id"])
        spool.add_rows("summaries", protein_summaries(proteins, group), order)
        for matrix, (tool, field) in MATRICES.items():
            record = group[tool]
            state = matrices[matrix]
            if record["status"] == "success":
                evidence = record["evidence"]
                state["features"].update(
                    value
                    for values in evidence["features"].get(field, {}).values()
                    for value in values
                )
                for feature, definition in (
                    evidence["definitions"].get(field, {}).items()
                ):
                    if (
                        feature in state["definitions"]
                        and state["definitions"][feature] != definition
                    ):
                        raise AnnotationError(f"Conflicting definition for {feature}")
                    state["definitions"][feature] = definition
                state["method_ids"].add(record["method_id"])
                state["resource_ids"].add(record["search"]["method"]["resource_id"])
            spool.add_counts(matrix, order[accession], feature_counts(record, field))
        spool.finish_sample()
        # Do not retain the last tool's evidence while loading the next sample.
        del group, proteins, record
        evidence = None
    for batch_id, members in batch_members.items():
        if not members:
            raise AnnotationError(f"Unused native eggNOG batch artefact: {batch_id}")
        if members != batches[batch_id].members.keys():
            raise AnnotationError(
                f"Missing result references for native eggNOG batch members: {batch_id}"
            )
    tables = work / "tables"
    matrix_root = tables / "functional_matrices"
    matrix_root.mkdir()
    for name in ("status", "provenance", "proteins", "coordinates", "summaries"):
        spool.writer(name)
    for name, columns in details.items():
        spool.schemas["detail:" + name] = columns
        spool.writer("detail:" + name)
    for handle in spool.handles.values():
        handle.flush()
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
        (status_index[acc] for acc in accessions),
    )
    catalogue = []
    for matrix, (tool, field) in MATRICES.items():
        state = matrices[matrix]
        features = sorted(state["features"])
        write_tsv(
            matrix_root / (matrix + ".tsv"),
            ("accession", *features),
            spool.matrix_rows(matrix, features, accessions),
        )
        for feature in features:
            catalogue.append(
                {
                    "matrix": matrix + ".tsv",
                    "feature_id": feature,
                    "definition": state["definitions"].get(feature),
                    "source": tool,
                    "resource_id": next(iter(state["resource_ids"])),
                    "method_id": next(iter(state["method_ids"])),
                }
            )
    write_tsv(
        matrix_root / "feature_catalogue.tsv",
        ("matrix", "feature_id", "definition", "source", "resource_id", "method_id"),
        catalogue,
    )
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "pipeline_contract": "nf-annotation-v0.5",
        "accessions": accessions,
        "enabled_tools": plan["enabled_tools"],
        "complete": not failures,
        "plan_id": plan["plan_id"],
        "bundles": bundle_records,
        "results": result_manifest,
        "cohort_artifacts": {},
        "tables": {
            str(path.relative_to(work)): digest(path)
            for path in sorted(tables.rglob("*.tsv"))
        },
    }
    if batches:
        manifest["native_batches"] = {
            batch_id: {
                "path": f"annotation_batches/{batch_id}",
                "batch_result_id": native.record["batch_result_id"],
            }
            for batch_id, native in batches.items()
        }
    # No acceptance or partial report is published before the cohort validates.
    output.mkdir()
    tables.rename(output / "tables")
    write_json(
        output / "annotation_acceptance.json",
        {"complete": not failures, "failures": failures, "plan_id": plan["plan_id"]},
    )
    write_json(output / "annotation_results.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    for name in (
        "plan",
        "master",
        "sample-status",
        "output",
        "bundle-list",
        "result-list",
        "batch-list",
    ):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--publication", type=Path)
    parser.add_argument(
        "--publication-mode",
        choices=("self_contained", "reference"),
        default="self_contained",
    )
    parser.add_argument("--cohort-entries", type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        manifest = aggregate(
            args.plan,
            args.master,
            args.sample_status,
            read_publication_paths(
                args.bundle_list, args.publication, args.publication_mode
            ),
            read_json(args.result_list),
            args.output,
            batch_dirs=read_path_list(args.batch_list),
        )
        from cohort_evidence import inventory as cohort_inventory
        from cohort_evidence import read_entries

        manifest["cohort_artifacts"] = cohort_inventory(
            read_entries(args.cohort_entries)
        )
        write_json(args.output / "annotation_results.json", manifest)
        if args.publication_mode == "reference":
            from publication_references import (
                REFERENCE_CONTRACT,
                REFERENCE_FILE,
                Publication,
            )

            if args.publication is None:
                raise AnnotationError(
                    "Reference aggregation requires the publication root"
                )
            view = Publication(args.publication, pending=True)
            assert view.document is not None
            if (
                set(view.document["samples"]) != set(manifest["accessions"])
                or manifest["enabled_tools"]
            ):
                raise AnnotationError(
                    "Reference aggregation requires matching membership and disabled tools"
                )
            manifest["pipeline_contract"] = REFERENCE_CONTRACT
            manifest["reference_manifest_sha256"] = digest(view.root / REFERENCE_FILE)
            write_json(args.output / "annotation_results.json", manifest)
    except (AnnotationError, OSError) as error:
        LOGGER.error("%s", error)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
