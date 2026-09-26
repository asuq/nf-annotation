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
    *,
    tools: tuple[str, ...] = TOOLS,
    previous_batch_ids: list[str] | None = None,
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
    if source is not None and "eggnog" in tools:
        if previous_batch_ids is None:
            previous_batches = validate_source(source)[1]
        else:
            from annotation_result import native_batch_index

            previous_batches = native_batch_index(
                [source / "annotation_batches" / key for key in previous_batch_ids]
            )
    else:
        previous_batches = {}
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
        if "eggnog" in tools and "eggnog" in checked["enabled_tools"] and successful
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
        for tool in tools:
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
                if tools != TOOLS:
                    name = "task_" + identity([accession, tool])
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
                    # Planner tasks may execute in node-local scratch. Keep the
                    # generated batch path relative to the staged plan so it
                    # remains valid after Nextflow copies declared outputs back.
                    input_directory = (
                        batchdir.relative_to(output)
                        if tools != TOOLS
                        else batchdir.resolve()
                    )
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


def execution_groups(samples: Path, receipt: Path, source: Path | None, size: int) -> dict:
    """Freeze small eggNOG readiness groups without waiting for new proteins."""
    if size < 1:
        raise AnnotationError("eggNOG group size must be positive")
    checked = read_json(receipt)
    validate_preflight(checked)
    accessions = [row["accession"] for row in read_tsv(samples, ["accession"])]
    if not accessions or len(set(accessions)) != len(accessions):
        raise AnnotationError("Invalid grouping sample manifest")
    for accession in accessions:
        validate_accession(accession)
    groups = []
    remaining = set(accessions)
    if "eggnog" in checked["enabled_tools"]:
        previous = validate_source(source)[1] if source is not None else {}
        for batch_id, native in sorted(previous.items()):
            members = sorted(remaining.intersection(native.members))
            if members:
                groups.append(dict(accessions=members, previous_batch_ids=[batch_id]))
                remaining.difference_update(members)
        new_members = sorted(remaining)
        for offset in range(0, len(new_members), size):
            groups.append(dict(accessions=new_members[offset:offset + size], previous_batch_ids=[]))
    for group in groups:
        group["key"] = identity(group)
    record = dict(accessions=accessions, preflight_id=checked["preflight_id"], groups=groups)
    record["groups_id"] = identity(record)
    return record


def plan_part(members: Path, bundles: list[Path], receipt: Path, output: Path,
              source: Path | None, scope: str) -> dict:
    """Plan one ready proteome or one ready eggNOG group with existing methods."""
    group = read_json(members)
    accessions = group["accessions"]
    if scope == "individual" and len(accessions) != 1:
        raise AnnotationError("Individual planning requires one accession")
    if scope == "eggnog" and group.get("key") != identity({key: value for key, value in group.items() if key != "key"}):
        raise AnnotationError("Changed eggNOG readiness group")
    sample_file = output.parent / "part_samples.tsv"
    write_tsv(sample_file, ("accession",), [dict(accession=value) for value in accessions])
    tools = ("eggnog",) if scope == "eggnog" else tuple(tool for tool in TOOLS if tool != "eggnog")
    return plan(sample_file, bundles, receipt, output, source, tools=tools,
                previous_batch_ids=group.get("previous_batch_ids", []))


def merge_plans(samples: Path, receipt: Path, fragments: list[Path], output: Path,
                source: Path | None) -> dict:
    """Build the complete cohort accounting after independently scheduled work."""
    checked = read_json(receipt)
    validate_preflight(checked)
    if source is not None:
        validate_source(source)
    accessions = [row["accession"] for row in read_tsv(samples, ["accession"])]
    if not accessions or len(set(accessions)) != len(accessions):
        raise AnnotationError("Invalid cohort planning manifest")
    accession_set = set(accessions)
    observed = {}
    for path in fragments:
        fragment = read_json(path)
        if fragment.get("preflight_id") != checked["preflight_id"] or fragment.get("plan_id") != identity({key: value for key, value in fragment.items() if key != "plan_id"}):
            raise AnnotationError("Changed or incompatible annotation plan fragment")
        for row in fragment["tasks"]:
            key = (row["accession"], row["tool"])
            if key[0] not in accession_set or key[1] not in TOOLS or key in observed:
                raise AnnotationError("Duplicate or foreign sample/tool plan entry")
            observed[key] = row
    tasks = []
    for accession in accessions:
        for tool in TOOLS:
            row = observed.get((accession, tool))
            if row is None:
                row = dict.fromkeys(PLAN_COLUMNS)
                enabled = tool in checked["enabled_tools"]
                row.update(accession=accession, tool=tool, action="skip",
                           status="upstream_failed" if enabled else "skipped_disabled",
                           reason="no_valid_protein_bundle" if enabled else "tool_disabled")
            tasks.append(row)
    record = dict(schema_version=SCHEMA_VERSION, accessions=accessions,
                  enabled_tools=checked["enabled_tools"], preflight_id=checked["preflight_id"], tasks=tasks)
    record["plan_id"] = identity(record)
    output.mkdir()
    write_json(output / "annotation_plan.json", record)
    write_tsv(output / "annotation_plan.tsv", PLAN_COLUMNS, tasks)
    return record


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
    groups = sub.add_parser("groups", allow_abbrev=False)
    groups.add_argument("--samples", type=Path, required=True)
    groups.add_argument("--preflight", type=Path, required=True)
    groups.add_argument("--source", type=Path)
    groups.add_argument("--size", type=int, default=8)
    groups.add_argument("--output", type=Path, required=True)
    part = sub.add_parser("plan-part", allow_abbrev=False)
    part.add_argument("--members", type=Path, required=True)
    part.add_argument("--bundle-list", type=Path, required=True)
    part.add_argument("--preflight", type=Path, required=True)
    part.add_argument("--source", type=Path)
    part.add_argument("--scope", choices=("individual", "eggnog"), required=True)
    part.add_argument("--output", type=Path, required=True)
    merger = sub.add_parser("merge-plans", allow_abbrev=False)
    merger.add_argument("--samples", type=Path, required=True)
    merger.add_argument("--preflight", type=Path, required=True)
    merger.add_argument("--fragment-list", type=Path, required=True)
    merger.add_argument("--source", type=Path)
    merger.add_argument("--output", type=Path, required=True)
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
        elif args.command == "groups":
            write_json(args.output, execution_groups(args.samples, args.preflight, args.source, args.size))
        elif args.command == "plan-part":
            plan_part(args.members, read_path_list(args.bundle_list), args.preflight,
                      args.output, args.source, args.scope)
        elif args.command == "merge-plans":
            merge_plans(args.samples, args.preflight, read_path_list(args.fragment_list), args.output, args.source)
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
