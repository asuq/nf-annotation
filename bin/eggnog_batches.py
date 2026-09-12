"""Bounded whole-proteome eggNOG inputs with portable, checked provenance."""

from __future__ import annotations

import hashlib
import math
import re
import shutil
import sqlite3
import tempfile
from collections.abc import Iterator
from contextlib import closing
from pathlib import Path
from typing import Any

from annotation_common import (
    PROTEIN_COLUMNS,
    AnnotationError,
    bundle_proteins,
    digest,
    identity,
    read_json,
    read_tsv,
    validate_accession,
    write_json,
    write_tsv,
)

TARGET_FASTA_BYTES = 4 * 1024 * 1024
BATCH_SCHEMA = "eggnog-whole-proteome-batch-v1"
MAPPING_COLUMNS = ("ordinal", "accession", "gene_id", "tool_id", "length", "sha256")
INPUT_FILES = ("input.faa", "protein_manifest.tsv", "mapping.tsv")


def packing_code_identity() -> dict[str, str]:
    """Identify the code that validates and constructs the native query batch."""
    return {
        name: digest(Path(__file__).with_name(name))
        for name in ("eggnog_batches.py", "annotation_common.py")
    }


def fasta_records(path: Path) -> Iterator[tuple[str, str]]:
    """Read canonical protein records without retaining a cohort FASTA in memory."""
    identifier, sequence = None, []
    with path.open(encoding="ascii", newline="") as handle:
        for line in handle:
            if not line.endswith("\n") or "\r" in line:
                raise AnnotationError(
                    "Batch FASTA requires complete LF-terminated lines"
                )
            value = line[:-1]
            if value.startswith(">"):
                if identifier is not None:
                    if not sequence:
                        raise AnnotationError("Empty batch protein sequence")
                    yield identifier, "".join(sequence)
                identifier, sequence = value[1:], []
                if not identifier or any(char.isspace() for char in identifier):
                    raise AnnotationError("Malformed canonical batch FASTA identifier")
            else:
                if (
                    identifier is None
                    or not value
                    or any(char.isspace() for char in value)
                ):
                    raise AnnotationError("Malformed canonical batch protein sequence")
                sequence.append(value)
    if identifier is None or not sequence:
        raise AnnotationError("Missing batch protein sequence")
    yield identifier, "".join(sequence)


def validate_proteins(path: Path, proteins: list[dict[str, str]]) -> None:
    """Require exact identity, order, sequence lengths and sequence SHA-256 values."""
    records = fasta_records(path)
    for protein in proteins:
        if set(protein) != set(PROTEIN_COLUMNS):
            raise AnnotationError("Unsupported batch protein-manifest schema")
        record = next(records, None)
        if record is None:
            raise AnnotationError(
                "Batch FASTA contains fewer proteins than its mapping"
            )
        identifier, sequence = record
        gene = f"{protein['accession']}::{protein['protein_id']}"
        if (
            protein["gene_id"] != gene
            or protein["tool_id"]
            != "p" + hashlib.sha256(gene.encode()).hexdigest()[:24]
            or identifier != protein["tool_id"]
            or protein["length"] != str(len(sequence))
            or protein["sha256"] != hashlib.sha256(sequence.encode("ascii")).hexdigest()
        ):
            raise AnnotationError(
                "Batch protein ID, order, length or sequence checksum mismatch"
            )
    if next(records, None) is not None:
        raise AnnotationError("Batch FASTA contains more proteins than its mapping")


def mapping_rows(proteins: list[dict[str, str]]) -> list[dict[str, Any]]:
    """Map each unchanged native query ID to one accession and canonical gene."""
    return [
        dict(ordinal=ordinal, **{key: row[key] for key in MAPPING_COLUMNS[1:]})
        for ordinal, row in enumerate(proteins)
    ]


def _packing(target: int) -> dict[str, Any]:
    if isinstance(target, bool) or not isinstance(target, int) or target < 1:
        raise AnnotationError("Batch target FASTA bytes must be a positive integer")
    return {
        "algorithm": "accession_sorted_next_fit_whole_proteomes",
        "target_fasta_bytes": target,
        "oversize_policy": "standalone_unmodified_proteome",
        "sequence_policy": "preserve_bytes_ids_order_no_deduplication",
    }


def validate_search_method(method: dict[str, Any]) -> None:
    """Require the immutable resource/runtime and explicit command from preflight."""
    if (
        not isinstance(method, dict)
        or set(method)
        != {"tool", "runtime_id", "resource_id", "command", "native_code"}
        or method["tool"] != "eggnog"
    ):
        raise AnnotationError("Batch inputs require the validated eggNOG search method")
    if (
        not isinstance(method["runtime_id"], str)
        or not re.fullmatch(r"(?:sif-)?sha256:[a-f0-9]{64}", method["runtime_id"])
        or not isinstance(method["resource_id"], str)
        or not re.fullmatch(r"[a-f0-9]{64}", method["resource_id"])
        or not isinstance(method["native_code"], dict)
    ):
        raise AnnotationError(
            "Batch search requires immutable runtime and resource identities"
        )
    command = method["command"]
    if not isinstance(command, dict) or set(command) != {
        "steps",
        "version_commands",
        "environment",
        "cpus",
        "memory_gib",
    }:
        raise AnnotationError("Batch search requires its explicit native command")
    if (
        type(command["cpus"]) is not int
        or command["cpus"] < 1
        or type(command["memory_gib"]) not in (int, float)
        or not math.isfinite(command["memory_gib"])
        or command["memory_gib"] <= 0
        or not isinstance(command["environment"], dict)
    ):
        raise AnnotationError("Invalid batch native execution settings")
    for name in ("steps", "version_commands"):
        if (
            not isinstance(command[name], list)
            or not command[name]
            or any(
                not isinstance(step, list)
                or not step
                or any(not isinstance(arg, str) or not arg for arg in step)
                for step in command[name]
            )
        ):
            raise AnnotationError("Malformed batch native command arguments")


def _write_batch(
    pending: list[tuple[Path, dict[str, Any], list[dict[str, str]]]],
    method: dict[str, Any],
    target: int,
    path: Path,
) -> dict[str, Any]:
    path.mkdir()
    proteins, members, total_bytes = [], [], 0
    with (path / "input.faa").open("wb") as output:
        for bundle, manifest, rows in pending:
            size = (bundle / "proteins.faa").stat().st_size
            members.append(
                dict(
                    **{
                        key: manifest[key]
                        for key in (
                            "accession",
                            "input_id",
                            "protein_count",
                            "genetic_code",
                            "source_genome_sha256",
                        )
                    },
                    fasta_sha256=manifest["files"]["proteins.faa"],
                    protein_manifest_sha256=manifest["files"]["protein_manifest.tsv"],
                    fasta_bytes=size,
                    byte_offset=total_bytes,
                    protein_offset=len(proteins),
                )
            )
            with (bundle / "proteins.faa").open("rb") as source:
                shutil.copyfileobj(source, output)
            total_bytes += size
            proteins.extend(rows)
    write_tsv(path / "protein_manifest.tsv", PROTEIN_COLUMNS, proteins)
    write_tsv(path / "mapping.tsv", MAPPING_COLUMNS, mapping_rows(proteins))
    inputs = {
        "schema": BATCH_SCHEMA,
        "packing": _packing(target),
        "packing_code": packing_code_identity(),
        "oversized": total_bytes > target,
        "fasta_bytes": total_bytes,
        "protein_count": len(proteins),
        "members": members,
        "files": {name: digest(path / name) for name in INPUT_FILES},
    }
    input_id = identity(inputs)
    record = dict(
        **inputs,
        input_id=input_id,
        search_method=method,
        search_fingerprint=identity({"input_id": input_id, "method": method}),
    )
    record["batch_id"] = identity(record)
    write_json(path / "batch.json", record)
    validate_batch(path)
    return record


def prepare_batches(
    bundles: list[Path],
    search_method: dict[str, Any],
    outdir: Path,
    *,
    target_fasta_bytes: int = TARGET_FASTA_BYTES,
) -> list[tuple[Path, dict[str, Any]]]:
    """Pack sorted complete proteomes, with an explicit oversized singleton rule.

    The caller supplies its validated immutable eggNOG search method. No native
    search or annotation is performed. A disk-backed uniqueness index bounds
    cohort-wide ID checking; only one proteome group is held in memory. Existing
    output directories are never reused or overwritten.
    """
    _packing(target_fasta_bytes)
    validate_search_method(search_method)
    ordered = [
        (read_json(bundle / "bundle.json").get("accession"), bundle)
        for bundle in bundles
    ]
    for accession, _ in ordered:
        validate_accession(accession)
    if not ordered or len({accession for accession, _ in ordered}) != len(ordered):
        raise AnnotationError(
            "Batch input cohort is empty or contains duplicate accessions"
        )
    ordered.sort(key=lambda item: item[0])
    outdir.mkdir(parents=True)
    result, pending, pending_bytes = [], [], 0

    def emit_batch() -> None:
        nonlocal pending, pending_bytes
        path = outdir / f"batch{len(result):08d}"
        result.append(
            (path, _write_batch(pending, search_method, target_fasta_bytes, path))
        )
        pending, pending_bytes = [], 0

    with (
        tempfile.TemporaryDirectory(
            prefix="eggnog-batch-identities-", dir=outdir
        ) as temporary,
        closing(sqlite3.connect(Path(temporary) / "identifiers.sqlite")) as database,
    ):
        database.execute(
            "CREATE TABLE proteins (tool_id TEXT PRIMARY KEY, gene_id TEXT UNIQUE)"
        )
        for accession, bundle in ordered:
            manifest, proteins = bundle_proteins(bundle)
            if manifest["accession"] != accession:
                raise AnnotationError("Bundle accession changed during batch planning")
            validate_proteins(bundle / "proteins.faa", proteins)
            try:
                database.executemany(
                    "INSERT INTO proteins (tool_id, gene_id) VALUES (?, ?)",
                    ((row["tool_id"], row["gene_id"]) for row in proteins),
                )
            except sqlite3.IntegrityError as error:
                raise AnnotationError(
                    "Duplicate cross-sample batch protein or gene ID"
                ) from error
            size = (bundle / "proteins.faa").stat().st_size
            if pending and pending_bytes + size > target_fasta_bytes:
                emit_batch()
            pending.append((bundle, manifest, proteins))
            pending_bytes += size
            if pending_bytes > target_fasta_bytes:
                emit_batch()
        if pending:
            emit_batch()
    return result


def validate_batch(
    path: Path, *, expected_search_fingerprint: str | None = None
) -> tuple[dict[str, Any], list[dict[str, str]]]:
    """Validate portable batch files and every member-to-sequence correspondence."""
    record = read_json(path / "batch.json")
    if (
        not isinstance(record, dict)
        or set(record)
        != {
            "schema",
            "packing",
            "packing_code",
            "oversized",
            "fasta_bytes",
            "protein_count",
            "members",
            "files",
            "input_id",
            "search_method",
            "search_fingerprint",
            "batch_id",
        }
        or record.get("schema") != BATCH_SCHEMA
        or record.get("batch_id")
        != identity({key: value for key, value in record.items() if key != "batch_id"})
    ):
        raise AnnotationError("Invalid eggNOG batch receipt")
    if (
        not isinstance(record["search_method"], dict)
        or not isinstance(record["packing"], dict)
        or type(record["protein_count"]) is not int
        or record["protein_count"] < 1
        or type(record["fasta_bytes"]) is not int
        or record["fasta_bytes"] < 1
    ):
        raise AnnotationError("Malformed eggNOG batch receipt")
    if record.get("packing_code") != packing_code_identity():
        raise AnnotationError("Batch packing or validation code changed")
    inputs = {
        key: value
        for key, value in record.items()
        if key not in ("input_id", "search_method", "search_fingerprint", "batch_id")
    }
    if record.get("input_id") != identity(inputs) or record.get(
        "search_fingerprint"
    ) != identity({"input_id": record["input_id"], "method": record["search_method"]}):
        raise AnnotationError("Invalid eggNOG batch input or search identity")
    if (
        expected_search_fingerprint is not None
        and record["search_fingerprint"] != expected_search_fingerprint
    ):
        raise AnnotationError("Batch search differs from the planned fingerprint")
    validate_search_method(record["search_method"])
    packing = record.get("packing", {})
    target = packing.get("target_fasta_bytes")
    if packing != _packing(target):
        raise AnnotationError("Unsupported eggNOG batch packing policy")
    children = list(path.iterdir())
    if (
        {child.name for child in children} != set(INPUT_FILES) | {"batch.json"}
        or any(child.is_symlink() or not child.is_file() for child in children)
        or record.get("files") != {name: digest(path / name) for name in INPUT_FILES}
    ):
        raise AnnotationError(
            "EggNOG batch input files changed or contain linked evidence"
        )
    proteins = read_tsv(path / "protein_manifest.tsv", PROTEIN_COLUMNS)
    mapping = read_tsv(path / "mapping.tsv", MAPPING_COLUMNS)
    expected_mapping = [
        {key: str(value) for key, value in row.items()}
        for row in mapping_rows(proteins)
    ]
    if (
        mapping != expected_mapping
        or len(proteins) != record["protein_count"]
        or not proteins
    ):
        raise AnnotationError("Batch mapping or protein count mismatch")
    for key in ("gene_id", "tool_id"):
        if len({row[key] for row in proteins}) != len(proteins):
            raise AnnotationError("Duplicate cross-sample batch protein or gene ID")
    validate_proteins(path / "input.faa", proteins)
    members = record.get("members")
    if (
        not isinstance(members, list)
        or not members
        or any(
            not isinstance(member, dict)
            or set(member)
            != {
                "accession",
                "input_id",
                "protein_count",
                "genetic_code",
                "source_genome_sha256",
                "fasta_sha256",
                "protein_manifest_sha256",
                "fasta_bytes",
                "byte_offset",
                "protein_offset",
            }
            for member in members
        )
    ):
        raise AnnotationError("Missing eggNOG batch membership")
    accessions = [member["accession"] for member in members]
    for accession in accessions:
        validate_accession(accession)
    if accessions != sorted(set(accessions)):
        raise AnnotationError("Duplicate or unsorted eggNOG batch membership")
    byte_offset, protein_offset = 0, 0
    with (path / "input.faa").open("rb") as source:
        for member in members:
            validate_accession(member["accession"])
            count, size = member["protein_count"], member["fasta_bytes"]
            if (
                type(count) is not int
                or count < 1
                or type(size) is not int
                or size < 1
                or type(member["protein_offset"]) is not int
                or type(member["byte_offset"]) is not int
                or type(member["genetic_code"]) is not int
                or member["genetic_code"] < 1
                or member["protein_offset"] != protein_offset
                or member["byte_offset"] != byte_offset
            ):
                raise AnnotationError("Invalid batch member boundaries")
            rows = proteins[protein_offset : protein_offset + count]
            if len(rows) != count or any(
                row["accession"] != member["accession"]
                or row["genetic_code"] != str(member["genetic_code"])
                for row in rows
            ):
                raise AnnotationError(
                    "Batch member protein or genetic-code mapping mismatch"
                )
            expected_input = identity(
                {
                    "accession": member["accession"],
                    "genetic_code": member["genetic_code"],
                    "proteins": [
                        {key: row[key] for key in ("gene_id", "tool_id", "sha256")}
                        for row in rows
                    ],
                }
            )
            if member["input_id"] != expected_input or any(
                not isinstance(member[name], str)
                or not re.fullmatch(r"[a-f0-9]{64}", member[name])
                for name in (
                    "source_genome_sha256",
                    "protein_manifest_sha256",
                    "fasta_sha256",
                )
            ):
                raise AnnotationError("Invalid batch member provenance")
            hasher, remaining = hashlib.sha256(), size
            while remaining:
                chunk = source.read(min(remaining, 1024 * 1024))
                if not chunk:
                    raise AnnotationError("Truncated batch member FASTA")
                hasher.update(chunk)
                remaining -= len(chunk)
            if hasher.hexdigest() != member["fasta_sha256"]:
                raise AnnotationError(
                    "Batch member FASTA bytes differ from its validated bundle"
                )
            byte_offset += size
            protein_offset += count
        if source.read(1):
            raise AnnotationError("Batch FASTA contains unassigned bytes")
    if byte_offset != record["fasta_bytes"] or protein_offset != len(proteins):
        raise AnnotationError("Batch members do not conserve input bytes or proteins")
    oversized = byte_offset > target
    if (
        type(record["oversized"]) is not bool
        or record["oversized"] != oversized
        or (oversized and len(members) != 1)
    ):
        raise AnnotationError(
            "Invalid oversized batch; only complete standalone proteomes are permitted"
        )
    return record, proteins
