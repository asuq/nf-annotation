"""Validate shared native searches and isolated per-proteome eggNOG annotation."""

from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

from annotation_common import (
    AnnotationError,
    digest,
    integer,
    number,
    query,
    read_json,
    validate_accession,
    write_json,
)

NATIVE_CODE_FILES = (
    "annotation_commands.py",
    "run_eggnog_batch.py",
    "eggnog_native.py",
    "eggnog_batches.py",
    "annotation_common.py",
)
EXECUTION_SCHEMA = "eggnog-shared-search-per-proteome-annotation-v1"
PARTITION_SCHEMA = "eggnog-derived-native-seed-partition-v1"
DERIVED_HEADER = b"# Derived partition of validated native batch seeds; no native completion footer.\n"
SEED_NATIVE_COLUMNS = (
    "qseqid",
    "sseqid",
    "evalue",
    "bitscore",
    "qstart",
    "qend",
    "sstart",
    "send",
    "pident",
    "qcov",
    "scov",
)


def query_lookup(proteins: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    """Resolve only canonical tool IDs, preserving sample-local original IDs."""
    if not proteins:
        raise AnnotationError("eggNOG input contains no declared proteins")
    lookup, genes = {}, set()
    for protein in proteins:
        identifier, gene = protein.get("tool_id"), protein.get("gene_id")
        validate_accession(protein.get("accession"))
        if (
            not isinstance(identifier, str)
            or not identifier
            or any(char.isspace() for char in identifier)
            or not isinstance(gene, str)
            or not gene
        ):
            raise AnnotationError("Invalid canonical eggNOG query or gene ID")
        if identifier in lookup or gene in genes:
            raise AnnotationError("Duplicate canonical eggNOG query or gene ID")
        lookup[identifier] = protein
        genes.add(gene)
    return lookup


def read_native_seeds(
    path: Path, lookup: dict[str, dict[str, str]]
) -> list[dict[str, str]]:
    """Read the complete native search table without changing its rows or footer."""
    text = path.read_text()
    lines = text.splitlines()
    if not text.endswith("\n") or [
        line for line in lines if line.startswith("#") and not line.startswith("##")
    ] != ["#" + "\t".join(SEED_NATIVE_COLUMNS)]:
        raise AnnotationError(
            "Missing native eggNOG seed header or complete final line"
        )
    data = [line.split("\t") for line in lines if line and not line.startswith("#")]
    if [line for line in lines if re.fullmatch(r"## [0-9]+ queries scanned", line)] != [
        f"## {len(data)} queries scanned"
    ]:
        raise AnnotationError(
            "Missing or inconsistent native eggNOG seed completion count"
        )
    seen, rows = set(), []
    for values in data:
        if len(values) != len(SEED_NATIVE_COLUMNS):
            raise AnnotationError("Truncated eggNOG seed ortholog row")
        row = dict(zip(SEED_NATIVE_COLUMNS, values, strict=True))
        protein = query(lookup, row["qseqid"])
        if row["qseqid"] in seen:
            raise AnnotationError("Duplicate eggNOG seed query")
        seen.add(row["qseqid"])
        integer(row["sseqid"], "eggNOG integer seed ID")
        for key in ("evalue", "bitscore"):
            number(row[key], key, minimum=0)
        for key in ("pident", "qcov", "scov"):
            if number(row[key], key, minimum=0) > 100:
                raise AnnotationError(
                    "eggNOG seed identity/coverage exceeds 100 percent"
                )
        for key in ("qstart", "qend", "sstart", "send"):
            integer(row[key], key, minimum=1)
        if not 1 <= int(row["qstart"]) <= int(row["qend"]) <= int(protein["length"]):
            raise AnnotationError("eggNOG seed coordinates exceed the input protein")
        rows.append(row)
    return rows


def member_name(index: int) -> str:
    """Keep native filenames independent of accession punctuation or spaces."""
    return f"member{index:08d}"


def stage_plan(batch: dict[str, Any]) -> list[dict[str, Any]]:
    """Resolve the recorded annotation template for each complete proteome."""
    commands = batch["search_method"]["command"]["native_stages"]
    stages = [
        {
            "name": "search",
            "accession": None,
            "directory": "search",
            "command": commands["search"],
        }
    ]
    for index, member in enumerate(batch["members"]):
        name = member_name(index)
        stages.append(
            {
                "name": "annotation",
                "accession": member["accession"],
                "directory": f"annotations/{name}",
                "command": [
                    arg.replace("{member}", name) for arg in commands["annotation"]
                ],
            }
        )
    for index, stage in enumerate(stages):
        command = stage["command"]
        expected = {
            "--cpu": str(batch["search_method"]["command"]["cpus"]),
            "--output_dir": "raw/" + stage["directory"],
            "--temp_dir": "scratch/" + stage["directory"],
        }
        if index == 0:
            expected.update({"-m": "diamond", "-i": "task/input.faa"})
            valid_mode = (
                command.count("--no_annot") == 1 and "--report_orthologs" not in command
            )
        else:
            member = member_name(index - 1)
            expected.update(
                {
                    "-m": "no_search",
                    "-i": f"raw/derived/{member}/input.faa",
                    "--annotate_hits_table": f"raw/derived/{member}/seeds.tsv",
                }
            )
            valid_mode = (
                command.count("--report_orthologs") == 1 and "--no_annot" not in command
            )
        if not valid_mode or any(
            command.count(flag) != 1
            or command.index(flag) + 1 == len(command)
            or command[command.index(flag) + 1] != value
            for flag, value in expected.items()
        ):
            raise AnnotationError(
                "Native phase does not use its declared search or proteome inputs"
            )
    return stages


def partition_rows(
    batch: dict[str, Any], proteins: list[dict[str, str]], rows: list[dict[str, str]]
) -> dict[str, list[dict[str, str]]]:
    """Assign every validated native seed to exactly one declared proteome."""
    lookup = query_lookup(proteins)
    partitions = {member["accession"]: [] for member in batch["members"]}
    for row in rows:
        accession = query(lookup, row["qseqid"])["accession"]
        if accession not in partitions:
            raise AnnotationError("Native seed belongs to an undeclared batch member")
        partitions[accession].append(row)
    if sum(map(len, partitions.values())) != len(rows):
        raise AnnotationError("Seed partitions do not conserve the native search rows")
    return partitions


def seed_payload(rows: list[dict[str, str]]) -> bytes:
    """Serialize unchanged native row values, with no invented native footer."""
    return b"".join(
        ("\t".join(row[name] for name in SEED_NATIVE_COLUMNS) + "\n").encode("ascii")
        for row in rows
    )


def partition_receipt(
    batch: dict[str, Any],
    member: dict[str, Any],
    source_sha: str,
    payload: bytes,
    count: int,
) -> dict[str, Any]:
    """Bind a derived seed table to its native search and unchanged input proteome."""
    return {
        "schema": PARTITION_SCHEMA,
        "batch_id": batch["batch_id"],
        "batch_input_id": batch["input_id"],
        "accession": member["accession"],
        "member_input_id": member["input_id"],
        "protein_count": member["protein_count"],
        "fasta_sha256": member["fasta_sha256"],
        "seed_rows": count,
        "native_seed_file": "../../search/eggnog.emapper.seed_orthologs",
        "native_seed_sha256": source_sha,
        "mapping_sha256": batch["files"]["mapping.tsv"],
        "derived_seed_file": "seeds.tsv",
        "derived_seed_sha256": hashlib.sha256(payload).hexdigest(),
    }


def write_partitions(
    inputs: Path, raw: Path, batch: dict[str, Any], proteins: list[dict[str, str]]
) -> None:
    """Retain exact input FASTA slices and explicit partitions of a native search."""
    source = raw / "search/eggnog.emapper.seed_orthologs"
    seeds = read_native_seeds(source, query_lookup(proteins))
    partitions = partition_rows(batch, proteins, seeds)
    source_sha = digest(source)
    (raw / "derived").mkdir()
    with (inputs / "input.faa").open("rb") as fasta:
        for index, member in enumerate(batch["members"]):
            destination = raw / "derived" / member_name(index)
            destination.mkdir()
            fasta.seek(member["byte_offset"])
            sequence = fasta.read(member["fasta_bytes"])
            if hashlib.sha256(sequence).hexdigest() != member["fasta_sha256"]:
                raise AnnotationError(
                    "Derived annotation FASTA differs from its proteome"
                )
            (destination / "input.faa").write_bytes(sequence)
            rows = partitions[member["accession"]]
            payload = DERIVED_HEADER + seed_payload(rows)
            (destination / "seeds.tsv").write_bytes(payload)
            write_json(
                destination / "partition.json",
                partition_receipt(batch, member, source_sha, payload, len(rows)),
            )


def validate_execution(
    raw: Path, batch: dict[str, Any], proteins: list[dict[str, str]]
) -> dict[str, list[dict[str, str]]]:
    """Require every native phase and verify all derived inputs against the search."""
    execution = read_json(raw / "execution.json")
    plans = stage_plan(batch)
    if (
        not isinstance(execution, dict)
        or set(execution) != {"schema", "batch_id", "stages"}
        or execution["schema"] != EXECUTION_SCHEMA
        or execution["batch_id"] != batch["batch_id"]
        or not isinstance(execution["stages"], list)
        or len(execution["stages"]) != len(plans)
    ):
        raise AnnotationError("Incomplete native search/annotation phase receipt")
    for observed, expected in zip(execution["stages"], plans, strict=True):
        if (
            not isinstance(observed, dict)
            or set(observed) != set(expected) | {"exit_code", "wall_seconds"}
            or any(observed[key] != value for key, value in expected.items())
            or type(observed["exit_code"]) is not int
            or observed["exit_code"] != 0
            or type(observed["wall_seconds"]) not in (int, float)
            or not math.isfinite(observed["wall_seconds"])
            or observed["wall_seconds"] < 0
            or (raw / expected["directory"] / "exit_code.txt").read_text() != "0\n"
        ):
            raise AnnotationError(
                "Native phase identity, command or exit status differs"
            )
    names = {member_name(index) for index in range(len(batch["members"]))}
    for directory in ("derived", "annotations"):
        parent = raw / directory
        if (
            not parent.is_dir()
            or parent.is_symlink()
            or {path.name for path in parent.iterdir()} != names
            or any(not path.is_dir() or path.is_symlink() for path in parent.iterdir())
        ):
            raise AnnotationError("Incomplete or foreign native annotation member grid")
    source = raw / "search/eggnog.emapper.seed_orthologs"
    rows = read_native_seeds(source, query_lookup(proteins))
    partitions = partition_rows(batch, proteins, rows)
    source_sha = digest(source)
    for index, member in enumerate(batch["members"]):
        directory = raw / "derived" / member_name(index)
        selected = partitions[member["accession"]]
        payload = DERIVED_HEADER + seed_payload(selected)
        expected = partition_receipt(batch, member, source_sha, payload, len(selected))
        receipt = (
            json.dumps(expected, indent=2, sort_keys=True, ensure_ascii=True) + "\n"
        ).encode()
        sorted_rows = sorted(
            selected, key=lambda row: (int(row["sseqid"]), row["qseqid"])
        )
        if (
            digest(directory / "input.faa") != member["fasta_sha256"]
            or (directory / "seeds.tsv").read_bytes() != payload
            or (directory / "partition.json").read_bytes() != receipt
            or (directory / "seeds.tsv.sorted").read_bytes()
            != seed_payload(sorted_rows)
        ):
            raise AnnotationError(
                "Derived or mapper-sorted seed inputs differ from the native partition"
            )
    return partitions
