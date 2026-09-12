"""Validate portable native annotation evidence, including shared eggNOG batches."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from annotation_common import (
    SCHEMA_VERSION,
    AnnotationError,
    digest,
    identity,
    read_json,
)
from eggnog_batches import validate_batch


def inventory(root: Path) -> dict[str, str]:
    """Inventory ordinary evidence files; reject symlinks and empty inventories."""
    if root.is_symlink():
        raise AnnotationError(f"Linked evidence is unsupported: {root}")
    files = {}
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise AnnotationError(f"Linked evidence is unsupported: {path}")
        if path.is_file():
            files[path.relative_to(root).as_posix()] = digest(path)
    if not files:
        raise AnnotationError(f"Missing evidence files: {root}")
    return files


def native_exit_code(raw: Path) -> int:
    """Read the exit status retained by the native command wrapper."""
    value = (raw / "exit_code.txt").read_text().strip()
    if not re.fullmatch(r"[0-9]+", value):
        raise AnnotationError("Missing or invalid native exit code")
    return int(value)


@dataclass(frozen=True)
class NativeBatch:
    """One batch validated during this invocation, with compact member metadata."""

    root: Path
    record: dict[str, Any]
    members: dict[str, dict[str, Any]]


def validate_native_batch(root: Path) -> NativeBatch:
    """Check an archived batch's original input, native evidence and identities."""
    record = read_json(root / "batch_result.json")
    if (
        not isinstance(record, dict)
        or set(record)
        != {
            "schema_version",
            "kind",
            "batch",
            "exit_code",
            "raw_files",
            "batch_result_id",
        }
        or record.get("schema_version") != SCHEMA_VERSION
        or record.get("kind") != "eggnog_native_batch"
        or record.get("batch_result_id")
        != identity(
            {key: value for key, value in record.items() if key != "batch_result_id"}
        )
    ):
        raise AnnotationError(f"Invalid native batch result: {root}")
    if (root / "inputs").is_symlink():
        raise AnnotationError("Linked native batch inputs are unsupported")
    batch, _ = validate_batch(root / "inputs")
    if record["batch"] != batch:
        raise AnnotationError("Native batch input receipt differs from its result")
    if inventory(root / "raw") != record["raw_files"]:
        raise AnnotationError(f"Native batch evidence changed: {root}")
    if type(record["exit_code"]) is not int or record["exit_code"] != native_exit_code(
        root / "raw"
    ):
        raise AnnotationError(
            "Native batch exit status differs from its retained evidence"
        )
    return NativeBatch(
        root, record, {member["accession"]: member for member in batch["members"]}
    )


def native_batch_index(paths: list[Path]) -> dict[str, NativeBatch]:
    """Validate each staged shared artefact once; never silently replace duplicates."""
    batches = {}
    for path in paths:
        batch = validate_native_batch(path)
        batch_id = batch.record["batch"]["batch_id"]
        if batch_id in batches:
            raise AnnotationError("Duplicate native eggNOG batch artefact")
        batches[batch_id] = batch
    return batches


def batch_search_identity(batch: dict[str, Any]) -> dict[str, str]:
    """Name the pooled query context included in every member's search identity."""
    return {key: batch[key] for key in ("batch_id", "input_id", "search_fingerprint")}


def validate_raw_evidence(
    root: Path,
    record: dict[str, Any],
    batches: dict[str, NativeBatch] | None = None,
) -> None:
    """Validate individual raw files or an explicit, already checked shared batch."""
    if record.get("status") == "success" and (
        type(record.get("exit_code")) is not int or record["exit_code"] != 0
    ):
        raise AnnotationError(
            "Successful annotation result has a nonzero or missing native exit status"
        )
    if "native_batch" not in record:
        if "batch" in record.get("search", {}):
            raise AnnotationError("Shared native eggNOG batch reference is missing")
        if inventory(root / "raw") != record["raw_files"]:
            raise AnnotationError(f"Native annotation evidence changed: {root}")
        return
    reference = record["native_batch"]
    if (
        record.get("tool") != "eggnog"
        or not isinstance(reference, dict)
        or set(reference) != {"batch_id", "batch_result_id"}
        or not isinstance(reference["batch_id"], str)
        or (root / "raw").exists()
        or (root / "raw").is_symlink()
    ):
        raise AnnotationError("Malformed shared native batch reference")
    if batches is None or reference["batch_id"] not in batches:
        raise AnnotationError("Shared native eggNOG batch is missing")
    native = batches[reference["batch_id"]]
    batch = native.record["batch"]
    member = native.members.get(record.get("accession"))
    if (
        reference["batch_result_id"] != native.record["batch_result_id"]
        or record.get("raw_files") != native.record["raw_files"]
        or member is None
        or record.get("input_id") != member["input_id"]
        or type(record.get("input_proteins")) is not int
        or record.get("input_proteins") != member["protein_count"]
        or record.get("search", {}).get("batch") != batch_search_identity(batch)
        or record.get("search", {}).get("method") != batch["search_method"]
        or record.get("search", {}).get("input_id") != member["input_id"]
        or type(record.get("search", {}).get("genetic_code")) is not int
        or record.get("search", {}).get("genetic_code") != member["genetic_code"]
        or record.get("search", {}).get("source_genome_sha256")
        != member["source_genome_sha256"]
        or type(record.get("exit_code")) is not int
        or record.get("exit_code") != native.record["exit_code"]
    ):
        raise AnnotationError("Member result differs from its shared native batch")


def validate_result(
    root: Path,
    *,
    normalized: bool = True,
    batches: dict[str, NativeBatch] | None = None,
) -> dict[str, Any]:
    """Validate a successful published result before reusing its evidence."""
    record = read_json(root / "result.json")
    if (
        not isinstance(record, dict)
        or record.get("schema_version") != SCHEMA_VERSION
        or record.get("result_id")
        != identity({key: value for key, value in record.items() if key != "result_id"})
    ):
        raise AnnotationError(f"Invalid annotation result record: {root}")
    if record["status"] != "success":
        raise AnnotationError(f"Cannot reuse unsuccessful annotation result: {root}")
    validate_raw_evidence(root, record, batches)
    if normalized and inventory(root / "normalized") != record["normalized_files"]:
        raise AnnotationError(f"Normalized annotation evidence changed: {root}")
    return record
