#!/usr/bin/env python3
"""Preflight, plan and validate the shared functional-annotation workflow."""

from __future__ import annotations

import argparse
import logging
import shutil
import sys
from pathlib import Path
from typing import Any

from annotation_batch_tasks import complete_batch, member_task_identity, plan_batch
from annotation_common import (
    SCHEMA_VERSION,
    TOOLS,
    AnnotationError,
    bundle_proteins,
    identity,
    read_json,
    read_tsv,
    validate_accession,
    write_json,
    write_tsv,
)
from annotation_path_lists import read_path_list
from annotation_resources import AnnotationResourceError
from annotation_result import validate_result
from annotation_source import import_source, validate_source
from annotation_tasks import (
    normalize_task,
    plan_task,
    preflight,
    validate_preflight,
)
from eggnog_batches import prepare_batches

PLAN_COLUMNS = (
    "accession",
    "tool",
    "action",
    "reason",
    "status",
    "input_proteins",
    "task_directory",
    "batch_id",
    "bundle",
    "resource",
    "container",
    "cpus",
    "memory_gib",
    "search_fingerprint",
    "normalization_fingerprint",
    "method_id",
)


def plan(
    samples: Path,
    bundles: list[Path],
    receipt: Path,
    output: Path,
    source: Path | None = None,
) -> dict[str, Any]:
    """Plan every declared accession/tool combination, including absent bundles."""
    checked = read_json(receipt)
    validate_preflight(checked)
    rows = read_tsv(samples, ["accession"])
    accessions = [row["accession"] for row in rows]
    accession_set = set(accessions)
    if not accessions or len(accession_set) != len(accessions):
        raise AnnotationError(
            "Annotation sample manifest is empty or contains duplicates"
        )
    for accession in accessions:
        validate_accession(accession)
    by_accession = {}
    for bundle in bundles:
        record = read_json(bundle / "bundle.json")
        accession = record.get("accession")
        if (
            record.get("schema_version") != SCHEMA_VERSION
            or accession not in accession_set
            or accession in by_accession
        ):
            raise AnnotationError("Unexpected or duplicate annotation bundle")
        if record["status"] not in ("success", "upstream_failed", "incompatible_input"):
            raise AnnotationError("Unsupported bundle status")
        by_accession[accession] = (bundle.resolve(), record)
    previous_batches = validate_source(source)[1] if source is not None else {}
    output.mkdir()
    successful = [
        bundle
        for bundle, record in by_accession.values()
        if record["status"] == "success"
    ]
    batches = (
        prepare_batches(
            successful,
            checked["tools"]["eggnog"]["search_method"],
            output / "batch_inputs",
        )
        if "eggnog" in checked["enabled_tools"] and successful
        else []
    )
    batch_for_accession = {
        member["accession"]: (directory, record)
        for directory, record in batches
        for member in record["members"]
    }
    eggnog_members = {}
    eggnog_rows = {}
    tasks = []
    for accession in accessions:
        bundle, bundle_record = by_accession.get(accession, (None, {}))
        proteins = []
        if bundle_record.get("status") == "success":
            bundle_record, proteins = bundle_proteins(bundle)
            if bundle_record["accession"] != accession:
                raise AnnotationError(
                    "Bundle accession changed during annotation planning"
                )
        for tool in TOOLS:
            row: dict[str, Any] = dict.fromkeys(PLAN_COLUMNS)
            row.update(
                accession=accession,
                tool=tool,
                action="skip",
                status="upstream_failed",
                reason="no_valid_protein_bundle",
                input_proteins=bundle_record.get("protein_count"),
            )
            if tool not in checked["enabled_tools"]:
                row.update(status="skipped_disabled", reason="tool_disabled")
            elif bundle_record.get("status") in (
                "incompatible_input",
                "upstream_failed",
            ):
                row.update(
                    status=bundle_record["status"], reason=bundle_record["reason"]
                )
            elif bundle_record.get("status") == "success":
                entry = checked["tools"][tool]
                name = f"task{len(tasks):08d}"
                old = (
                    source / "samples" / accession / "annotation" / tool
                    if source
                    else None
                )
                if tool == "eggnog":
                    batchdir, batch = batch_for_accession[accession]
                    task = member_task_identity(bundle_record, entry, batch)
                    eggnog_members[accession] = task
                    eggnog_rows[accession] = row
                    name = "eggnog_" + batch["batch_id"]
                    row["batch_id"] = batch["batch_id"]
                    input_directory = batchdir.resolve()
                else:
                    task = plan_task(
                        bundle,
                        tool,
                        entry,
                        output / name,
                        old,
                        bundle_data=(bundle_record, proteins),
                    )
                    input_directory = bundle
                row.update(
                    action=task.get("action"),
                    reason=task.get("reason"),
                    status="planned",
                    task_directory=name,
                    bundle=str(input_directory),
                    resource=entry["resource_path"],
                    container=entry["container"],
                    cpus=entry["search_method"]["command"]["cpus"],
                    memory_gib=entry["search_method"]["command"]["memory_gib"],
                    search_fingerprint=task["search_fingerprint"],
                    normalization_fingerprint=task["normalization_fingerprint"],
                    method_id=task["method_id"],
                )
            tasks.append(row)
    for directory, batch in batches:
        member_accessions = [member["accession"] for member in batch["members"]]
        task = plan_batch(
            directory,
            checked["tools"]["eggnog"],
            [eggnog_members[accession] for accession in member_accessions],
            output / ("eggnog_" + batch["batch_id"]),
            previous_batches.get(batch["batch_id"]),
            {
                accession: source / "samples" / accession / "annotation" / "eggnog"
                for accession in member_accessions
            }
            if source is not None
            else {},
        )
        for accession in member_accessions:
            eggnog_rows[accession].update(action=task["action"], reason=task["reason"])
    record = dict(
        schema_version=SCHEMA_VERSION,
        accessions=accessions,
        enabled_tools=checked["enabled_tools"],
        preflight_id=checked["preflight_id"],
        tasks=tasks,
    )
    record["plan_id"] = identity(record)
    write_json(output / "annotation_plan.json", record)
    write_tsv(output / "annotation_plan.tsv", PLAN_COLUMNS, tasks)
    return record


def reuse(taskdir: Path, outdir: Path) -> None:
    """Republish a checksummed normalized result while recording this run's action."""
    task = read_json(taskdir / "task.json")
    record = validate_result(taskdir / "previous")
    if any(
        record[key] != task[key]
        for key in (
            "accession",
            "tool",
            "search_fingerprint",
            "normalization_fingerprint",
        )
    ):
        raise AnnotationError("Reused result differs from the current task")
    shutil.copytree(taskdir / "previous", outdir)
    record["action"] = "reuse"
    record["result_id"] = identity(
        {key: value for key, value in record.items() if key != "result_id"}
    )
    write_json(outdir / "result.json", record)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    check = sub.add_parser("preflight")
    check.add_argument("--config", type=Path, required=True)
    check.add_argument("--output", type=Path, required=True)
    planner = sub.add_parser("plan", allow_abbrev=False)
    planner.add_argument("--samples", type=Path, required=True)
    planner.add_argument("--bundle-list", type=Path, required=True)
    planner.add_argument("--preflight", type=Path, required=True)
    planner.add_argument("--source", type=Path)
    planner.add_argument("--output", type=Path, required=True)
    normalizer = sub.add_parser("normalize")
    for name in ("task", "raw", "bundle", "output"):
        normalizer.add_argument("--" + name, type=Path, required=True)
    normalizer.add_argument("--resource", type=Path)
    reuser = sub.add_parser("reuse")
    reuser.add_argument("--task", type=Path, required=True)
    reuser.add_argument("--output", type=Path, required=True)
    batcher = sub.add_parser("complete-batch", allow_abbrev=False)
    for name in ("task", "raw", "batch", "resource", "output"):
        batcher.add_argument("--" + name, type=Path, required=True)
    importer = sub.add_parser("import")
    importer.add_argument("--source", type=Path, required=True)
    importer.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        if args.command == "preflight":
            write_json(args.output, preflight(read_json(args.config)))
        elif args.command == "plan":
            plan(
                args.samples,
                read_path_list(args.bundle_list),
                args.preflight,
                args.output,
                args.source,
            )
        elif args.command == "normalize":
            normalize_task(args.task, args.raw, args.bundle, args.resource, args.output)
        elif args.command == "reuse":
            reuse(args.task, args.output)
        elif args.command == "complete-batch":
            complete_batch(args.task, args.raw, args.batch, args.resource, args.output)
        else:
            import_source(args.source, args.output)
    except (AnnotationError, AnnotationResourceError, OSError) as error:
        logging.error("%s", error)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
