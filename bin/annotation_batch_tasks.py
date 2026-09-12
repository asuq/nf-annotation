"""Plan and complete whole-proteome eggNOG batches with shared native evidence."""

from __future__ import annotations

import re
import shutil
from pathlib import Path
from typing import Any

import normalize_eggnog
from annotation_commands import shell_script
from annotation_common import (
    SCHEMA_VERSION,
    AnnotationError,
    digest,
    identity,
    read_json,
    write_json,
)
from annotation_resources import RESOURCE_FILE, validate_resource
from annotation_result import (
    NativeBatch,
    batch_search_identity,
    inventory,
    native_exit_code,
    validate_native_batch,
    validate_result,
)
from annotation_tasks import code_identity, task_identity
from eggnog_batches import validate_batch

IDENTITY_FIELDS = frozenset(
    {
        "schema_version",
        "accession",
        "tool",
        "input_id",
        "coordinate_id",
        "input_proteins",
        "search",
        "search_fingerprint",
        "normalization_fingerprint",
        "method_id",
        "interpretation",
        "resource",
        "container",
    }
)


def member_task_identity(
    manifest: dict[str, Any], entry: dict[str, Any], batch_record: dict[str, Any]
) -> dict[str, Any]:
    """Preserve the sample identity while naming its exact pooled query context."""
    if entry["search_method"] != batch_record["search_method"]:
        raise AnnotationError("Member method differs from the current eggNOG batch")
    member = task_identity(manifest, "eggnog", entry)
    member["search"]["batch"] = batch_search_identity(batch_record)
    member["search_fingerprint"] = identity(member["search"])
    member["normalization_fingerprint"] = identity(
        dict(
            search_fingerprint=member["search_fingerprint"],
            interpretation=member["interpretation"],
        )
    )
    return member


def _validate_members(
    batch: dict[str, Any], members: list[dict[str, Any]], entry: dict[str, Any]
) -> None:
    """Require the complete ordered member grid and its current interpretation."""
    if (
        not isinstance(members, list)
        or len(members) != len(batch["members"])
        or entry["search_method"] != batch["search_method"]
        or entry["interpretation"].get("schema_version") != SCHEMA_VERSION
        or entry["interpretation"].get("code") != code_identity("eggnog")
        or entry["method_id"]
        != identity(
            dict(search=entry["search_method"], interpretation=entry["interpretation"])
        )
    ):
        raise AnnotationError("Invalid or changed eggNOG batch member method/code")
    for member, canonical in zip(members, batch["members"], strict=True):
        if (
            not isinstance(member, dict)
            or set(member) != IDENTITY_FIELDS
            or not isinstance(member["coordinate_id"], str)
            or not re.fullmatch(r"[a-f0-9]{64}", member["coordinate_id"])
        ):
            raise AnnotationError("Malformed eggNOG batch member task")
        manifest = dict(canonical, coordinate_id=member["coordinate_id"])
        if member != member_task_identity(manifest, entry, batch):
            raise AnnotationError("Task member differs from the current eggNOG batch")


def _entry_from_member(member: dict[str, Any]) -> dict[str, Any]:
    """Recover the shared planned entry without adding external path dependencies."""
    try:
        resource = member["resource"]
        return dict(
            search_method=member["search"]["method"],
            interpretation=member["interpretation"],
            method_id=member["method_id"],
            container=member["container"],
            resource_manifest_sha256=resource["manifest_sha256"],
            resource=dict(
                resource_id=resource["resource_id"],
                contract=dict(
                    version=resource["version"], settings=resource["settings"]
                ),
            ),
        )
    except (KeyError, TypeError) as error:
        raise AnnotationError("Malformed eggNOG batch member provenance") from error


def _copy_native(previous: NativeBatch, destination: Path) -> NativeBatch:
    """Check both sides of copying an archived packet selected by its immutable ID."""
    current = validate_native_batch(previous.root)
    if current.record != previous.record:
        raise AnnotationError("Native eggNOG batch changed after selection")
    destination.mkdir(parents=True)
    for name in ("inputs", "raw"):
        shutil.copytree(previous.root / name, destination / name)
    shutil.copyfile(
        previous.root / "batch_result.json", destination / "batch_result.json"
    )
    copied = validate_native_batch(destination)
    if copied.record != previous.record:
        raise AnnotationError(
            "Copied native eggNOG batch differs from the selected packet"
        )
    return copied


def plan_batch(
    batchdir: Path,
    entry: dict[str, Any],
    members: list[dict[str, Any]],
    outdir: Path,
    previous: NativeBatch | None,
    previous_results: dict[str, Path],
) -> dict[str, Any]:
    """Choose one run/reuse/renormalize action for the entire validated batch."""
    batch, _ = validate_batch(batchdir)
    _validate_members(batch, members, entry)
    action, reason = "run", "no_compatible_native_batch"
    selected: dict[str, tuple[Path, dict[str, Any]]] = {}
    compatible = (
        previous is not None
        and previous.record["batch"] == batch
        and previous.record["exit_code"] == 0
    )
    if compatible:
        action, reason = "renormalize", "member_interpretation_or_result_changed"
        for member in members:
            old_path = previous_results.get(member["accession"])
            if old_path is None or not (old_path / "result.json").is_file():
                continue
            old = read_json(old_path / "result.json")
            if (
                not isinstance(old, dict)
                or old.get("schema_version") != SCHEMA_VERSION
                or old.get("result_id")
                != identity(
                    {key: value for key, value in old.items() if key != "result_id"}
                )
            ):
                raise AnnotationError("Invalid previous eggNOG member result")
            if (
                old.get("status") != "success"
                or old.get("search_fingerprint") != member["search_fingerprint"]
                or "native_batch" not in old
            ):
                continue
            reuse_normalized = (
                old.get("normalization_fingerprint")
                == member["normalization_fingerprint"]
            )
            checked = validate_result(
                old_path,
                normalized=reuse_normalized,
                batches={batch["batch_id"]: previous},
            )
            if (
                checked.get("accession") != member["accession"]
                or checked.get("tool") != "eggnog"
                or checked.get("search") != member["search"]
                or (
                    reuse_normalized
                    and any(
                        checked.get(key) != member[key]
                        for key in ("interpretation", "method_id", "resource")
                    )
                )
            ):
                raise AnnotationError(
                    "Previous eggNOG result has a mismatched member identity"
                )
            if reuse_normalized:
                selected[member["accession"]] = old_path, checked
        if len(selected) == len(members):
            action, reason = "reuse", "matching_batch_search_and_member_interpretation"
    task = dict(
        schema_version=SCHEMA_VERSION,
        kind="eggnog_batch",
        batch_id=batch["batch_id"],
        action=action,
        reason=reason,
        members=[dict(member, action=action, reason=reason) for member in members],
    )
    outdir.mkdir(parents=True)
    if action == "run":
        shutil.copyfile(batchdir / "input.faa", outdir / "input.faa")
        if digest(outdir / "input.faa") != batch["files"]["input.faa"]:
            raise AnnotationError("Copied eggNOG batch input differs from its receipt")
        (outdir / "run.sh").write_text(shell_script(entry["search_method"]["command"]))
    else:
        copied = _copy_native(previous, outdir / "previous_batch")
        task["native_batch"] = dict(
            batch_id=batch["batch_id"], batch_result_id=copied.record["batch_result_id"]
        )
        if action == "reuse":
            task["previous_result_ids"] = {}
            for accession, (old_path, old) in selected.items():
                destination = outdir / "previous_results" / accession
                destination.mkdir(parents=True)
                shutil.copytree(old_path / "normalized", destination / "normalized")
                shutil.copyfile(old_path / "result.json", destination / "result.json")
                checked = validate_result(
                    destination, batches={batch["batch_id"]: copied}
                )
                if checked != old:
                    raise AnnotationError(
                        "Copied normalized member differs from its selected result"
                    )
                task["previous_result_ids"][accession] = old["result_id"]
    task["task_id"] = identity(task)
    write_json(outdir / "task.json", task)
    return task


def _validate_task(
    task: dict[str, Any], batch: dict[str, Any]
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Validate the checksummed task and all action-specific integrity anchors."""
    if not isinstance(task, dict):
        raise AnnotationError("Invalid eggNOG batch task")
    required = {
        "schema_version",
        "kind",
        "batch_id",
        "action",
        "reason",
        "members",
        "task_id",
    }
    if task.get("action") in ("reuse", "renormalize"):
        required.add("native_batch")
    if task.get("action") == "reuse":
        required.add("previous_result_ids")
    if (
        set(task) != required
        or task.get("schema_version") != SCHEMA_VERSION
        or task.get("kind") != "eggnog_batch"
        or task.get("batch_id") != batch["batch_id"]
        or task.get("action") not in ("run", "reuse", "renormalize")
        or task.get("task_id")
        != identity({key: value for key, value in task.items() if key != "task_id"})
        or not isinstance(task.get("members"), list)
        or not task["members"]
    ):
        raise AnnotationError("Invalid or mismatched eggNOG batch task")
    members = []
    for member in task["members"]:
        if (
            not isinstance(member, dict)
            or set(member) != IDENTITY_FIELDS | {"action", "reason"}
            or member["action"] != task["action"]
            or member["reason"] != task["reason"]
        ):
            raise AnnotationError("Mixed or malformed eggNOG batch member tasks")
        members.append(
            {key: value for key, value in member.items() if key in IDENTITY_FIELDS}
        )
    entry = _entry_from_member(members[0])
    _validate_members(batch, members, entry)
    if task["action"] != "run":
        reference = task["native_batch"]
        if (
            not isinstance(reference, dict)
            or set(reference) != {"batch_id", "batch_result_id"}
            or reference["batch_id"] != batch["batch_id"]
            or not isinstance(reference["batch_result_id"], str)
            or not re.fullmatch(r"[a-f0-9]{64}", reference["batch_result_id"])
        ):
            raise AnnotationError("Invalid planned native eggNOG batch reference")
    if task["action"] == "reuse":
        selected = task["previous_result_ids"]
        if (
            not isinstance(selected, dict)
            or set(selected) != {member["accession"] for member in members}
            or any(
                not isinstance(value, str) or not re.fullmatch(r"[a-f0-9]{64}", value)
                for value in selected.values()
            )
        ):
            raise AnnotationError("Invalid planned normalized member result identities")
    return entry, members


def complete_batch(
    taskdir: Path, raw: Path, batchdir: Path, resource: Path, outdir: Path
) -> dict[str, Any]:
    """Publish one immutable native packet and either all successful or failed members."""
    # Nextflow intentionally stages these directory roots as symlinks. Resolve
    # only the declared roots; validators continue to reject links within them.
    taskdir, raw, batchdir, resource = (
        path.resolve(strict=True) for path in (taskdir, raw, batchdir, resource)
    )
    batch, proteins = validate_batch(batchdir)
    task = read_json(taskdir / "task.json")
    entry, members = _validate_task(task, batch)
    prepared = validate_resource(resource, "eggnog", verify_checksums=False)
    if (
        prepared["resource_id"] != batch["search_method"]["resource_id"]
        or prepared["resource_id"] != entry["resource"]["resource_id"]
        or digest(resource / RESOURCE_FILE) != entry["resource_manifest_sha256"]
        or any(
            prepared["contract"][key] != entry["resource"]["contract"][key]
            for key in ("version", "settings")
        )
    ):
        raise AnnotationError(
            "Task/resource identity differs from the planned eggNOG batch"
        )
    action = task["action"]
    previous, old_results = None, {}
    if action == "run":
        if (
            (taskdir / "input.faa").is_symlink()
            or (taskdir / "run.sh").is_symlink()
            or digest(taskdir / "input.faa") != batch["files"]["input.faa"]
            or (taskdir / "run.sh").read_text()
            != shell_script(entry["search_method"]["command"])
        ):
            raise AnnotationError("Planned eggNOG input or native command changed")
    else:
        previous = validate_native_batch(taskdir / "previous_batch")
        if (
            previous.record["batch"] != batch
            or previous.record["batch_result_id"]
            != task["native_batch"]["batch_result_id"]
            or previous.record["exit_code"] != 0
        ):
            raise AnnotationError(
                "Archived native eggNOG batch differs from the planned packet"
            )
        if action == "reuse":
            for member in members:
                accession = member["accession"]
                old = validate_result(
                    taskdir / "previous_results" / accession,
                    batches={batch["batch_id"]: previous},
                )
                if old["result_id"] != task["previous_result_ids"][accession] or any(
                    old.get(key) != member[key]
                    for key in (
                        "accession",
                        "tool",
                        "search_fingerprint",
                        "normalization_fingerprint",
                    )
                ):
                    raise AnnotationError(
                        "Archived normalized member differs from the planned result"
                    )
                old_results[accession] = old
    raw_files = inventory(raw)
    exit_code = native_exit_code(raw)
    if previous is not None and (
        raw_files != previous.record["raw_files"]
        or exit_code != previous.record["exit_code"]
    ):
        raise AnnotationError(
            "Supplied raw evidence differs from the planned native batch"
        )
    outdir.mkdir()
    packet_root = outdir / "annotation_batches" / batch["batch_id"]
    packet_root.mkdir(parents=True)
    shutil.copytree(batchdir, packet_root / "inputs")
    shutil.copytree(raw, packet_root / "raw")
    copied_batch, _ = validate_batch(packet_root / "inputs")
    if copied_batch != batch or inventory(packet_root / "raw") != raw_files:
        raise AnnotationError(
            "Copied native packet differs from the validated batch input/evidence"
        )
    packet = dict(
        schema_version=SCHEMA_VERSION,
        kind="eggnog_native_batch",
        batch=batch,
        exit_code=exit_code,
        raw_files=raw_files,
    )
    packet["batch_result_id"] = identity(packet)
    if previous is not None and packet != previous.record:
        raise AnnotationError("Reused native packet identity changed")
    write_json(packet_root / "batch_result.json", packet)
    error_message, normalized = None, {}
    if exit_code != 0:
        error_message = f"Native eggnog exited with status {exit_code}"
    elif action != "reuse":
        try:
            normalized = normalize_eggnog.normalize_batch(
                packet_root / "raw", proteins, resource
            )
            if list(normalized) != [member["accession"] for member in members]:
                raise AnnotationError(
                    "Batch normalization returned an incomplete or reordered member grid"
                )
        except (AnnotationError, OSError) as error:
            error_message = str(error)
    results = {}
    for member in members:
        accession = member["accession"]
        destination = outdir / "samples" / accession / "annotation" / "eggnog"
        destination.mkdir(parents=True)
        record = dict(
            member,
            action=action,
            status="failed",
            reason=error_message,
            exit_code=exit_code,
            reported_hits=None,
            accepted_proteins=None,
            evidence=None,
            normalized_files={},
            raw_files=raw_files,
            native_batch=dict(
                batch_id=batch["batch_id"], batch_result_id=packet["batch_result_id"]
            ),
        )
        if error_message is not None:
            (destination / "validation_error.txt").write_text(error_message + "\n")
        elif action == "reuse":
            old = old_results[accession]
            shutil.copytree(
                taskdir / "previous_results" / accession / "normalized",
                destination / "normalized",
            )
            if inventory(destination / "normalized") != old["normalized_files"]:
                raise AnnotationError("Copied normalized member evidence changed")
            record.update(
                {
                    key: old[key]
                    for key in (
                        "reported_hits",
                        "accepted_proteins",
                        "evidence",
                        "normalized_files",
                    )
                }
            )
            record.update(status="success", reason="validated_native_output")
        else:
            evidence = normalized[accession]
            summary = evidence.publish(destination / "normalized")
            record.update(
                status="success",
                reason="validated_native_output",
                reported_hits=evidence.reported_hits,
                accepted_proteins=len(evidence.accepted),
                evidence=summary,
                normalized_files=inventory(destination / "normalized"),
            )
        record["result_id"] = identity(record)
        write_json(destination / "result.json", record)
        results[accession] = record
    return dict(
        batch_id=batch["batch_id"],
        batch_result_id=packet["batch_result_id"],
        action=action,
        results=results,
    )
