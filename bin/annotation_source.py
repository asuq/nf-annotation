"""Validate native published v0.4 results without a Nextflow work directory."""

from __future__ import annotations

import csv
import re
import shutil
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from annotation_common import (
    SCHEMA_VERSION,
    TOOLS,
    AnnotationError,
    bundle_proteins,
    digest,
    identity,
    read_json,
    validate_accession,
)
from annotation_result import NativeBatch, native_batch_index, validate_raw_evidence
from annotation_summary import ANNOTATION_COLUMNS, MATRICES
from publication_references import REFERENCE_CONTRACTS, Publication, strict_path
from validate_inputs import detect_metadata_key_column


def source_roots(source: Path | list[Path] | None) -> list[Path]:
    """Normalise source paths without changing their caller-visible order."""
    roots = [] if source is None else source if isinstance(source, list) else [source]
    resolved = [path.resolve(strict=True) for path in roots]
    if len(set(resolved)) != len(resolved):
        raise AnnotationError("Duplicate published source directories")
    return resolved


def validate_sources(source: Path | list[Path] | None) -> dict[str, NativeBatch]:
    """Validate independent publications and combine their native batch index."""
    accessions: set[str] = set()
    batches: dict[str, NativeBatch] = {}
    for root in source_roots(source):
        manifest, native = validate_source(root)
        overlap = accessions.intersection(manifest["accessions"])
        if overlap:
            raise AnnotationError(f"Duplicate accessions across published sources: {sorted(overlap)}")
        if batches.keys() & native.keys():
            raise AnnotationError("Duplicate native batch IDs across published sources")
        accessions.update(manifest["accessions"])
        batches.update(native)
    return batches


def sample_source(roots: list[Path], accession: str) -> Path | None:
    """Route a sample to its only source; never choose by source precedence."""
    matches = [root for root in roots if (root / "samples" / accession).exists()]
    if len(matches) > 1:
        raise AnnotationError(f"Duplicate accession across published sources: {accession}")
    return matches[0] if matches else None


@contextmanager
def source_rows(
    path: Path, required: Iterable[str]
) -> Iterator[tuple[list[str], Iterator[dict[str, str]]]]:
    """Validate complete TSV rows while callers retain only required metadata.

    Match the shared TSV reader's quoting, row-shape and large-field rules without
    retaining unused evidence cells or changing normalizer code fingerprints.
    """
    previous_limit = csv.field_size_limit()
    csv.field_size_limit(max(previous_limit, path.stat().st_size))
    try:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = reader.fieldnames
            if (
                not header
                or any(not name for name in header)
                or len(header) != len(set(header))
                or not set(required) <= set(header)
            ):
                raise AnnotationError(
                    f"Invalid table header in {path}; required: {list(required)}"
                )

            def rows() -> Iterator[dict[str, str]]:
                for line, row in enumerate(reader, 2):
                    if None in row or None in row.values():
                        raise AnnotationError(f"Malformed table row in {path}:{line}")
                    yield row

            yield header, rows()
    finally:
        csv.field_size_limit(previous_limit)


def source_path(root: Path, name: str) -> Path:
    """Resolve a declared regular artefact without traversing links or parents."""
    return strict_path(root, name)


def source_batches(
    root: Path, manifest: dict[str, Any], accessions: list[str], publication: Publication | None = None
) -> dict[str, NativeBatch]:
    """Validate each declared portable shared archive once and match packet IDs."""
    path_for = publication.path if publication is not None else lambda name: source_path(root, name)
    entries = manifest.get("native_batches", {})
    if not isinstance(entries, dict):
        raise AnnotationError("Invalid published native batch mapping")
    directory = path_for("annotation_batches")
    if directory.exists() and (
        not directory.is_dir()
        or {path.name for path in directory.iterdir()} != set(entries)
    ):
        raise AnnotationError(
            "Published native batch archives differ from the manifest"
        )
    paths = []
    for batch_id, entry in entries.items():
        if (
            not isinstance(batch_id, str)
            or not re.fullmatch(r"[a-f0-9]{64}", batch_id)
            or not isinstance(entry, dict)
            or set(entry) != {"path", "batch_result_id"}
            or entry["path"] != f"annotation_batches/{batch_id}"
            or not isinstance(entry["batch_result_id"], str)
            or not re.fullmatch(r"[a-f0-9]{64}", entry["batch_result_id"])
        ):
            raise AnnotationError("Invalid published native batch reference or path")
        path = path_for(entry["path"])
        if not path.is_dir():
            raise AnnotationError("Missing published native batch archive")
        path_for(entry["path"] + "/batch_result.json")
        paths.append(path)
    batches = native_batch_index(paths)
    if set(batches) != set(entries):
        raise AnnotationError(
            "Published native batch packet IDs differ from the manifest"
        )
    declared_accessions = set(accessions)
    for batch_id, native in batches.items():
        if native.record["batch_result_id"] != entries[batch_id]["batch_result_id"]:
            raise AnnotationError("Published native batch result identity changed")
        if set(native.members) - declared_accessions:
            raise AnnotationError("Published native batch contains a foreign accession")
    return batches


def validate_source(root: Path, *, _content_checked: bool = False, _validated: dict | None = None) -> tuple[dict[str, Any], dict[str, NativeBatch]]:
    """Verify source artefacts and return its once-validated native batch index."""
    cache_key = root.resolve(strict=True)
    validated = {} if _validated is None else _validated
    if cache_key in validated:
        return validated[cache_key]
    publication = Publication(root)
    if publication.reference:
        if not _content_checked:
            from reuse_source_validation import fingerprint
            fingerprint(root)
        for provider in publication.sources:
            validate_source(provider.root, _content_checked=True, _validated=validated)
    path_for = publication.path if publication.reference else lambda name: source_path(root, name)
    manifest = read_json(root / "annotation_results.json")
    if (
        manifest.get("schema_version") != SCHEMA_VERSION
        or manifest.get("pipeline_contract") not in ("nf-annotation-v0.4", *REFERENCE_CONTRACTS)
    ):
        raise AnnotationError("Annotation reuse requires native v0.4 published results")
    accessions = manifest.get("accessions", [])
    if not accessions or len(accessions) != len(set(accessions)):
        raise AnnotationError("Invalid published annotation accession set")
    for accession in accessions:
        validate_accession(accession)
    if any(
        set(manifest.get(key, {})) - set(accessions) for key in ("bundles", "results")
    ):
        raise AnnotationError("Published manifest contains an undeclared accession")
    required = {
        "tables/master_table.tsv",
        "tables/sample_status.tsv",
        "tables/annotation_status.tsv",
        "tables/annotation_provenance.tsv",
        "tables/protein_manifest.tsv",
        "tables/gene_coordinates.tsv",
        "tables/protein_function_summary.tsv",
        "tables/eggnog_seed_hits.tsv",
        "tables/eggnog_annotations.tsv",
        "tables/cog_assignments.tsv",
        "tables/pfam_domains.tsv",
        "tables/kofam_hits.tsv",
        "tables/defence_systems.tsv",
        "tables/defence_genes.tsv",
        "tables/functional_matrices/feature_catalogue.tsv",
        *[f"tables/functional_matrices/{name}.tsv" for name in MATRICES],
    }
    if not required <= set(manifest.get("tables", {})):
        raise AnnotationError("Published annotation manifest omits required tables")
    for name, expected in manifest["tables"].items():
        path = path_for(name)
        if (
            not name.startswith("tables/")
            or not path.is_file()
            or digest(path) != expected
        ):
            raise AnnotationError(f"Published annotation table changed: {name}")
    tool_statuses = [f"{tool}_status" for tool in TOOLS]
    with source_rows(
        root / "tables/master_table.tsv", ["Gcode", *ANNOTATION_COLUMNS]
    ) as (header, rows):
        master_key = detect_metadata_key_column(header)
        masters = [
            {key: row[key] for key in (master_key, "Gcode", *tool_statuses)}
            for row in rows
        ]
    if not masters:
        raise AnnotationError("Published master table has no genomes")
    with source_rows(
        root / "tables/sample_status.tsv",
        ["accession", *tool_statuses],
    ) as (_, rows):
        statuses = [
            {key: row[key] for key in ("accession", *tool_statuses)} for row in rows
        ]
    if [row[master_key] for row in masters] != accessions or [
        row["accession"] for row in statuses
    ] != accessions:
        raise AnnotationError(
            "Published master/status accession order differs from manifest"
        )
    grid = {}
    with source_rows(
        root / "tables/annotation_status.tsv", ["accession", "tool", "status"]
    ) as (_, rows):
        for row in rows:
            key = row["accession"], row["tool"]
            if key in grid:
                raise AnnotationError(
                    "Duplicate published annotation accession/tool status"
                )
            grid[key] = row["status"]
    if set(grid) != {(acc, tool) for acc in accessions for tool in TOOLS}:
        raise AnnotationError(
            "Published annotation status does not contain the complete declared grid"
        )
    batches = source_batches(root, manifest, accessions, publication if publication.reference else None)
    referenced_members = {batch_id: set() for batch_id in batches}
    for master, status in zip(masters, statuses, strict=True):
        accession = master[master_key]
        for tool in TOOLS:
            if (
                grid[(accession, tool)] != status[f"{tool}_status"]
                or status[f"{tool}_status"] != master[f"{tool}_status"]
            ):
                raise AnnotationError(
                    "Published annotation status contradicts a genome summary"
                )
        bundle_record = manifest.get("bundles", {}).get(accession)
        metadata = None
        if bundle_record is not None:
            expected_path = f"samples/{accession}/annotation/bundle"
            if bundle_record.get("path") != expected_path:
                raise AnnotationError("Unexpected published bundle path")
            bundle = path_for(expected_path)
            bundle_json = path_for(expected_path + "/bundle.json")
            if digest(bundle_json) != bundle_record.get("manifest_sha256"):
                raise AnnotationError("Published bundle manifest changed")
            data = read_json(bundle_json)
            if (
                data.get("schema_version") != SCHEMA_VERSION
                or data.get("accession") != accession
                or data.get("status")
                not in ("success", "upstream_failed", "incompatible_input")
            ):
                raise AnnotationError("Invalid published bundle record")
            if data.get("status") == "success":
                metadata, _ = bundle_proteins(bundle)
                if (
                    str(metadata["genetic_code"]) != master["Gcode"]
                    or metadata["accession"] != accession
                ):
                    raise AnnotationError(
                        "Published bundle differs from source genetic code or accession"
                    )
        for tool, result_record in (
            manifest.get("results", {}).get(accession, {}).items()
        ):
            expected_path = f"samples/{accession}/annotation/{tool}"
            if tool not in TOOLS or result_record.get("path") != expected_path:
                raise AnnotationError("Unexpected published annotation result path")
            record = read_json(path_for(expected_path + "/result.json"))
            if (
                record.get("schema_version") != SCHEMA_VERSION
                or record.get("accession") != accession
                or record.get("tool") != tool
                or record.get("result_id") != result_record.get("result_id")
                or record.get("result_id")
                != identity(
                    {key: value for key, value in record.items() if key != "result_id"}
                )
                or record.get("status") != grid[(accession, tool)]
            ):
                raise AnnotationError("Published result identity or status changed")
            validate_raw_evidence(path_for(expected_path), record, batches)
            if "native_batch" in record:
                batch_id = record["native_batch"]["batch_id"]
                member = batches[batch_id].members[accession]
                if metadata is None or any(
                    metadata[key] != member[key]
                    for key in (
                        "input_id",
                        "genetic_code",
                        "source_genome_sha256",
                        "protein_count",
                    )
                ):
                    raise AnnotationError(
                        "Published bundle differs from its native batch member"
                    )
                referenced_members[batch_id].add(accession)
        for tool in TOOLS:
            if grid[(accession, tool)] == "success" and (
                bundle_record is None
                or tool not in manifest.get("results", {}).get(accession, {})
            ):
                raise AnnotationError(
                    "Successful annotation lacks a published bundle or result"
                )
    for batch_id, native in batches.items():
        if not referenced_members[batch_id]:
            raise AnnotationError("Published native batch is unused")
        if referenced_members[batch_id] != set(native.members):
            raise AnnotationError(
                "Published native batch omits a member result reference"
            )
    validated[cache_key] = (manifest, batches)
    return manifest, batches


def import_source(root: Path, output: Path) -> None:
    """Retain published upstream evidence while new functional results are built."""
    publication = Publication(root)
    snapshot = None
    if publication.reference:
        from reuse_source_validation import fingerprint
        snapshot = fingerprint(root)
    manifest, _ = validate_source(root, _content_checked=publication.reference)
    for accession in manifest["accessions"]:
        sample = publication.native_root(accession)
        for path in sample.rglob("*"):
            if path.is_symlink():
                raise AnnotationError(f"Published sample artefact is a symlink: {path}")
    output.mkdir()
    from annotation_common import write_json
    write_json(output / "source_snapshot.json", snapshot or {"root": str(root.resolve()), "mode": "self_contained"})
    inherited = output / "inherited"
    (inherited / "tables").mkdir(parents=True)
    for name, destination in (
        ("master_table.tsv", "upstream_master.tsv"),
        ("sample_status.tsv", "upstream_sample_status.tsv"),
        ("validated_samples.tsv", "source_samples.tsv"),
    ):
        shutil.copyfile(source_path(root, "tables/" + name), output / destination)
    for path in (root / "tables").iterdir():
        if path.is_file() and f"tables/{path.name}" not in manifest["tables"]:
            source_path(root, "tables/" + path.name)
            shutil.copyfile(path, inherited / "tables" / path.name)
    for accession in manifest["accessions"]:
        source_sample = publication.native_root(accession)
        target_sample = inherited / "samples" / accession
        target_sample.mkdir(parents=True)
        for path in source_sample.iterdir():
            if path.name == "annotation":
                continue
            target = target_sample / path.name
            if path.is_dir():
                shutil.copytree(path, target)
            else:
                shutil.copyfile(path, target)
        if accession in manifest["bundles"]:
            shutil.copytree(
                publication.path(f"samples/{accession}/annotation/bundle"), target_sample / "annotation/bundle"
            )
