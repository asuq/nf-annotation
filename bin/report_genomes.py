#!/usr/bin/env python3
"""Publish genome-level reports without rebuilding inherited protein tables."""

from __future__ import annotations

import argparse
from pathlib import Path

from aggregate_annotations import PROVENANCE_COLUMNS, provenance_rows, record_rows
from annotation_common import (
    SCHEMA_VERSION,
    STATUS_COLUMNS,
    TOOLS,
    AnnotationError,
    digest,
    identity,
    read_json,
    read_tsv,
    write_json,
    write_tsv,
)
from annotation_path_lists import read_path_list
from annotation_summary import genome_summary
from cohort_evidence import inventory as cohort_inventory
from cohort_evidence import read_entries
from publication_references import PORTABLE_CONTRACT, REFERENCE_CONTRACT
from validate_inputs import detect_metadata_key_column


def report(
    plan_path: Path,
    master: Path,
    status: Path,
    bundles: list[Path],
    output: Path,
    reference: Path | None = None,
    cohort_entries=(),
) -> dict:
    plan = read_json(plan_path)
    if (
        plan.get("schema_version") != SCHEMA_VERSION
        or plan.get("enabled_tools")
        or plan.get("plan_id")
        != identity({k: v for k, v in plan.items() if k != "plan_id"})
    ):
        raise AnnotationError(
            "Genome reporting requires a valid disabled annotation plan"
        )
    masters, statuses = read_tsv(master), read_tsv(status)
    master_key = detect_metadata_key_column(list(masters[0]))
    accessions = [row[master_key] for row in masters]
    status_index = {row["accession"]: row for row in statuses}
    if (
        len(set(accessions)) != len(accessions)
        or set(accessions) != set(plan["accessions"])
        or set(status_index) != set(accessions)
        or len(status_index) != len(statuses)
    ):
        raise AnnotationError("Genome report membership differs from its inputs")
    intended = {(row["accession"], row["tool"]): row for row in plan["tasks"]}
    if (
        len(intended) != len(plan["tasks"])
        or set(intended) != {(a, t) for a in accessions for t in TOOLS}
        or any(
            row["action"] != "skip" or row["status"] != "skipped_disabled"
            for row in intended.values()
        )
    ):
        raise AnnotationError("Genome report requires the complete disabled tool grid")
    records, bundle_manifest = {}, {}
    for path in bundles:
        metadata = read_json(path)
        accession = metadata.get("accession")
        if accession not in status_index or accession in records:
            raise AnnotationError("Duplicate or foreign genome bundle")
        records[accession] = {
            key: metadata.get(key) for key in ("status", "protein_count")
        }
        bundle_manifest[accession] = {
            "path": f"samples/{accession}/annotation/bundle",
            "manifest_sha256": digest(path),
        }
    annotation_status, provenance = [], []
    for row in masters:
        accession = row[master_key]
        group = {
            tool: {
                "accession": accession,
                "tool": tool,
                "status": "skipped_disabled",
                "reason": intended[accession, tool]["reason"],
                "input_proteins": intended[accession, tool]["input_proteins"],
                "evidence": None,
            }
            for tool in TOOLS
        }
        metadata = records.get(accession, {})
        row.update(
            genome_summary(
                group,
                metadata.get("protein_count")
                if metadata.get("status") == "success"
                else None,
                True,
            )
        )
        status_index[accession].update(
            {f"{tool}_status": "skipped_disabled" for tool in TOOLS}
        )
        status_index[accession]["annotation_complete"] = "true"
        annotation_status.extend(record_rows(group))
        provenance.extend(provenance_rows(group))
    tables = output / "tables"
    tables.mkdir(parents=True)
    write_tsv(tables / "master_table.tsv", list(masters[0]), masters)
    write_tsv(
        tables / "sample_status.tsv",
        list(status_index[accessions[0]]),
        (status_index[a] for a in accessions),
    )
    write_tsv(tables / "annotation_status.tsv", STATUS_COLUMNS, annotation_status)
    write_tsv(tables / "annotation_provenance.tsv", PROVENANCE_COLUMNS, provenance)
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "pipeline_contract": REFERENCE_CONTRACT if reference else PORTABLE_CONTRACT,
        "report_scope": "genome",
        "accessions": accessions,
        "enabled_tools": [],
        "complete": True,
        "plan_id": plan["plan_id"],
        "bundles": bundle_manifest,
        "results": {},
        "cohort_artifacts": cohort_inventory(cohort_entries),
        "tables": {f"tables/{p.name}": digest(p) for p in sorted(tables.iterdir())},
    }
    if reference:
        manifest["reference_manifest_sha256"] = digest(reference)
    write_json(output / "annotation_results.json", manifest)
    write_json(
        output / "annotation_acceptance.json",
        {"complete": True, "failures": [], "plan_id": plan["plan_id"]},
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("plan", "master", "status", "bundles", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--cohort-entries", type=Path)
    args = parser.parse_args()
    try:
        report(
            args.plan,
            args.master,
            args.status,
            read_path_list(args.bundles),
            args.output,
            args.reference,
            read_entries(args.cohort_entries),
        )
    except (AnnotationError, OSError, ValueError) as error:
        parser.exit(1, str(error) + "\n")


if __name__ == "__main__":
    main()
