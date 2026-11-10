#!/usr/bin/env python3
"""Recover failed eggNOG annotation stages from a sealed native batch packet."""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True

from annotation_batch_tasks import IDENTITY_FIELDS, _validate_task
from annotation_commands import shell_script
from annotation_common import AnnotationError, digest, identity, read_json, write_json
from annotation_normalization import Normalized
from annotation_resources import RESOURCE_FILE, validate_resource
from annotation_result import inventory, validate_native_batch, validate_raw_evidence
from eggnog_batches import validate_batch
from eggnog_database import stage_database
from eggnog_native import (
    DERIVED_HEADER,
    EXECUTION_SCHEMA,
    NATIVE_CODE_FILES,
    member_name,
    partition_receipt,
    partition_rows,
    query_lookup,
    read_native_seeds,
    seed_payload,
    stage_plan,
    validate_execution,
)
from normalize_eggnog import _normalize_annotations, _seed_evidence, ontology_terms

RECOVERY_SCHEMA = "eggnog-annotation-stage-recovery-v1"


def validate_partial(
    raw: Path, batch: dict, proteins: list[dict], resource: Path
) -> list[dict]:
    """Reuse only a completed ordered prefix with intact per-proteome evidence."""
    execution = read_json(raw / "execution.json")
    if not isinstance(execution, dict):
        raise AnnotationError("Invalid partial eggNOG execution receipt")
    plans = stage_plan(batch)
    stages = execution.get("stages")
    if (
        set(execution) != {"schema", "batch_id", "stages"}
        or execution["schema"] != EXECUTION_SCHEMA
        or execution["batch_id"] != batch["batch_id"]
        or not isinstance(stages, list)
        or not 1 <= len(stages) <= len(plans)
    ):
        raise AnnotationError("Invalid partial eggNOG execution receipt")
    for index, observed in enumerate(stages):
        expected = plans[index]
        if (
            not isinstance(observed, dict)
            or set(observed) != set(expected) | {"exit_code", "wall_seconds"}
            or any(observed[key] != value for key, value in expected.items())
            or type(observed["exit_code"]) is not int
            or not 0 <= observed["exit_code"] <= 255
            or (observed["exit_code"] != 0 and index != len(stages) - 1)
            or type(observed["wall_seconds"]) not in (int, float)
            or not math.isfinite(observed["wall_seconds"])
            or observed["wall_seconds"] < 0
            or (raw / expected["directory"] / "exit_code.txt").read_text()
            != f"{observed['exit_code']}\n"
        ):
            raise AnnotationError(
                "Invalid partial native stage identity or exit status"
            )
    completed = [stage for stage in stages if stage["exit_code"] == 0]
    if (
        not completed
        or completed[0]["name"] != "search"
        or len(completed) == len(plans)
    ):
        raise AnnotationError(
            "Recovery requires a completed search and unfinished annotations"
        )
    if not (raw / "search/eggnog.emapper.hits").is_file():
        raise AnnotationError("Recovery requires the retained native search hits")
    observed_names = {Path(stage["directory"]).name for stage in stages[1:]}
    annotations = raw / "annotations"
    actual_names = (
        {path.name for path in annotations.iterdir()} if annotations.exists() else set()
    )
    if actual_names != observed_names:
        raise AnnotationError("Unrecorded or missing native member outputs")
    names = {member_name(index) for index in range(len(batch["members"]))}
    if {path.name for path in (raw / "derived").iterdir()} != names:
        raise AnnotationError("Incomplete or foreign derived seed member grid")
    lookup = query_lookup(proteins)
    native = raw / "search/eggnog.emapper.seed_orthologs"
    rows = read_native_seeds(native, lookup)
    partitions = partition_rows(batch, proteins, rows)
    source_sha = digest(native)
    for index, member in enumerate(batch["members"]):
        directory = raw / "derived" / member_name(index)
        selected = partitions[member["accession"]]
        payload = DERIVED_HEADER + seed_payload(selected)
        receipt = partition_receipt(batch, member, source_sha, payload, len(selected))
        expected_receipt = (
            json.dumps(receipt, indent=2, sort_keys=True, ensure_ascii=True) + "\n"
        ).encode()
        sorted_file = directory / "seeds.tsv.sorted"
        if (
            digest(directory / "input.faa") != member["fasta_sha256"]
            or (directory / "seeds.tsv").read_bytes() != payload
            or (directory / "partition.json").read_bytes() != expected_receipt
            or (index < len(completed) - 1 and not sorted_file.is_file())
            or (
                index < len(completed) - 1
                and sorted_file.exists()
                and sorted_file.read_bytes()
                != seed_payload(
                    sorted(
                        selected, key=lambda row: (int(row["sseqid"]), row["qseqid"])
                    )
                )
            )
        ):
            raise AnnotationError("Retained proteome or seed partition changed")
    seeds = _seed_evidence(rows, lookup, resource, Normalized())
    ontology = ontology_terms(resource / "go-basic.obo")
    for stage in completed[1:]:
        accession = stage["accession"]
        selected = partitions[accession]
        member_lookup = {
            key: value
            for key, value in lookup.items()
            if value["accession"] == accession
        }
        _normalize_annotations(
            raw / stage["directory"],
            member_lookup,
            {row["qseqid"]: seeds[row["qseqid"]] for row in selected},
            Normalized(),
            ontology,
        )
        orthologs = (raw / stage["directory"] / "eggnog.emapper.orthologs").read_text()
        lines = orthologs.splitlines()
        counts = [
            line
            for line in lines
            if line.startswith("## ") and line.endswith(" queries scanned")
        ]
        annotation_counts = [
            line
            for line in (raw / stage["directory"] / "eggnog.emapper.annotations")
            .read_text()
            .splitlines()
            if line.startswith("## ") and line.endswith(" queries scanned")
        ]
        if (
            not orthologs.endswith("\n")
            or [
                line
                for line in lines
                if line.startswith("#") and not line.startswith("##")
            ]
            != ["#query\torth_type\tspecies\torthologs"]
            or counts != annotation_counts
            or any(
                len(line.split("\t")) != 4 or line.split("\t")[0] not in member_lookup
                for line in lines
                if line and not line.startswith("#")
            )
        ):
            raise AnnotationError("Incomplete or foreign native ortholog evidence")
    return completed


def validate_source(
    taskdir: Path, batchdir: Path, previous: Path, members: list[Path], resource: Path
):
    """Match the frozen failed packet and every member to the unchanged current task."""
    batch, proteins = validate_batch(batchdir)
    task = read_json(taskdir / "task.json")
    entry, planned = _validate_task(task, batch)
    old = validate_native_batch(previous)
    if (
        task["action"] != "run"
        or old.record["batch"] != batch
        or old.record["exit_code"] == 0
    ):
        raise AnnotationError("Recovery source is not the failed current batch")
    if (
        (taskdir / "run.sh").read_text()
        != shell_script(entry["search_method"]["command"])
        or digest(taskdir / "input.faa") != batch["files"]["input.faa"]
        or read_json(taskdir / "native_commands.json")
        != batch["search_method"]["command"]["native_stages"]
        or any(
            digest(taskdir / name) != batch["search_method"]["native_code"].get(name)
            or digest(Path(__file__).with_name(name))
            != batch["search_method"]["native_code"].get(name)
            for name in NATIVE_CODE_FILES
        )
    ):
        raise AnnotationError("Recovery task input or execution method changed")
    prepared = validate_resource(resource, "eggnog", verify_checksums=False)
    if (
        prepared["resource_id"] != batch["search_method"]["resource_id"]
        or digest(resource / RESOURCE_FILE) != entry["resource_manifest_sha256"]
    ):
        raise AnnotationError("Recovery resource identity changed")
    records = []
    for path in members:
        record = read_json(path / "result.json")
        if (
            record.get("result_id")
            != identity(
                {key: value for key, value in record.items() if key != "result_id"}
            )
            or record.get("status") != "failed"
            or record.get("exit_code") != old.record["exit_code"]
        ):
            raise AnnotationError("Invalid failed eggNOG member receipt")
        validate_raw_evidence(path, record, {batch["batch_id"]: old})
        records.append(record)
    records.sort(key=lambda record: record["accession"])
    expected = sorted(planned, key=lambda record: record["accession"])
    if len(records) != len(expected) or any(
        any(record.get(key) != member[key] for key in IDENTITY_FIELDS)
        for record, member in zip(records, expected, strict=True)
    ):
        raise AnnotationError("Recovery member identity or coordinate mapping changed")
    return old, batch, proteins, prepared


def recover(
    taskdir: Path, batchdir: Path, previous: Path, members: list[Path], resource: Path
) -> int:
    """Execute only failed/unstarted member commands, retaining real attempt records."""
    old, batch, proteins, prepared = validate_source(
        taskdir, batchdir, previous, members, resource
    )
    raw = Path("raw")
    raw.mkdir()
    Path("scratch").mkdir()
    # Copy into this attempt before interpreting any prior evidence; never mutate the source.
    shutil.copytree(previous / "raw", raw, dirs_exist_ok=True)
    if inventory(raw) != old.record["raw_files"]:
        raise AnnotationError("Recovery copy differs from the sealed source packet")
    archive = raw / "failed_attempts" / old.record["batch_result_id"]
    archive.mkdir(parents=True)
    write_json(archive / "batch_result.json", old.record)
    for name in (
        "execution.json",
        "exit_code.txt",
        "versions.txt",
        "tool.log",
        "recovery.json",
    ):
        if (raw / name).exists():
            shutil.copyfile(raw / name, archive / name)
    local_database = Path("eggnog_data")
    database_files = stage_database(resource.resolve(), prepared, local_database)
    completed = validate_partial(raw, batch, proteins, local_database)
    receipt = {
        "schema": RECOVERY_SCHEMA,
        "source_batch_result_id": old.record["batch_result_id"],
        "code_sha256": digest(Path(__file__)),
        "runtime_id": batch["search_method"]["runtime_id"],
        "reused_stages": [stage["directory"] for stage in completed],
        "executed_stages": [],
        "database_files": database_files,
    }
    write_json(raw / "recovery.json", receipt)
    execution = {
        "schema": EXECUTION_SCHEMA,
        "batch_id": batch["batch_id"],
        "stages": completed,
    }
    write_json(raw / "execution.json", execution)
    # Keep a real nonzero status until every remaining native stage has completed.
    (raw / "exit_code.txt").write_text("1\n")
    with (raw / "tool.log").open("w") as log:
        environment = dict(
            os.environ, **batch["search_method"]["command"]["environment"]
        )
        for stage in stage_plan(batch)[len(completed) :]:
            sorted_seeds = Path(
                stage["command"][stage["command"].index("--annotate_hits_table") + 1]
                + ".sorted"
            )
            if sorted_seeds.exists():
                saved_seeds = archive / sorted_seeds.relative_to(raw)
                saved_seeds.parent.mkdir(parents=True, exist_ok=True)
                sorted_seeds.rename(saved_seeds)
            destination = raw / stage["directory"]
            if destination.exists():
                saved = archive / stage["directory"]
                saved.parent.mkdir(parents=True, exist_ok=True)
                destination.rename(saved)
            destination.mkdir(parents=True, exist_ok=True)
            temporary = Path(stage["command"][stage["command"].index("--temp_dir") + 1])
            temporary.mkdir(parents=True)
            started = time.monotonic()
            with (destination / "tool.log").open("w") as output:
                result = subprocess.run(
                    stage["command"],
                    env=environment,
                    stdout=output,
                    stderr=subprocess.STDOUT,
                    check=False,
                )
            status = (
                result.returncode if result.returncode >= 0 else 128 - result.returncode
            )
            (destination / "exit_code.txt").write_text(f"{status}\n")
            record = dict(
                stage, exit_code=status, wall_seconds=time.monotonic() - started
            )
            execution["stages"].append(record)
            receipt["executed_stages"].append(stage["directory"])
            write_json(raw / "execution.json", execution)
            write_json(raw / "recovery.json", receipt)
            print(json.dumps(record, sort_keys=True), file=log, flush=True)
            if status:
                (raw / "exit_code.txt").write_text(f"{status}\n")
                return status
        validate_execution(raw, batch, proteins)
        (raw / "exit_code.txt").write_text("0\n")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    for name in ("task", "batch", "previous", "resource"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--members", type=Path, nargs="+", required=True)
    args = parser.parse_args()
    if any(Path(name).exists() for name in ("raw", "scratch", "eggnog_data")):
        print("Recovery requires a fresh task directory", file=sys.stderr)
        return 1
    try:
        return recover(
            args.task.resolve(strict=True),
            args.batch.resolve(strict=True),
            args.previous.resolve(strict=True),
            args.members,
            args.resource.resolve(strict=True),
        )
    except (OSError, ValueError, RuntimeError) as error:
        print(str(error), file=sys.stderr)
        raw = Path("raw")
        if raw.is_dir():
            (raw / "exit_code.txt").write_text("1\n")
            with (raw / "tool.log").open("a") as handle:
                print(str(error), file=handle)
        return 1


if __name__ == "__main__":
    sys.exit(main())
