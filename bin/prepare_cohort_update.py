#!/usr/bin/env python3
"""Validate a published cohort before reusing its per-sample results."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import math
import re
import sys
from pathlib import Path
from typing import Sequence

import build_master_table
import build_sample_status
import collect_versions
import master_table_contract
import summarise_16s
import summarise_busco
import summarise_checkm2
import validate_inputs

LOGGER = logging.getLogger(__name__)
AUDIT_COLUMNS = (
    "accession",
    "action",
    "internal_id",
    "genome_sha256",
    "source_outdir",
    "source_manifest_sha256",
    "source_versions_sha256",
)


class CohortUpdateError(ValueError):
    """Raised when published results cannot be reused without ambiguity."""


def read_table(
    path: Path, required: Sequence[str] = ()
) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV strictly, including header, row-width and required-column checks."""
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if (
            not header
            or any(not name.strip() for name in header)
            or len(set(header)) != len(header)
        ):
            raise CohortUpdateError(f"Missing, empty or duplicate TSV headers: {path}")
        missing = set(required) - set(header)
        if missing:
            raise CohortUpdateError(
                f"Missing columns in {path}: {', '.join(sorted(missing))}"
            )
        rows = []
        for number, values in enumerate(reader, start=2):
            if len(values) != len(header):
                raise CohortUpdateError(f"Malformed TSV row {number} in {path}")
            rows.append(dict(zip(header, values, strict=True)))
    return header, rows


def index_rows(
    rows: Sequence[dict[str, str]], key: str, path: Path
) -> dict[str, dict[str, str]]:
    """Index records without dropping duplicate or empty accessions."""
    result = {}
    for row in rows:
        accession = row[key]
        validate_inputs.ensure_path_safe_accession(accession)
        if accession in result:
            raise CohortUpdateError(f"Duplicate accession {accession!r} in {path}")
        result[accession] = row
    return result


def single_summary(
    path: Path, accession: str, columns: Sequence[str]
) -> dict[str, str]:
    """Validate an accession-keyed per-sample summary."""
    _header, rows = read_table(path, ("accession", *columns))
    if len(rows) != 1 or rows[0]["accession"] != accession:
        raise CohortUpdateError(f"Expected exactly one row for {accession!r} in {path}")
    return rows[0]


def file_sha256(path: Path) -> str:
    """Hash a file without loading it into memory."""
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def genome_fingerprint(path: Path) -> str:
    """Hash ordered FASTA IDs and sequences, ignoring wrapping and descriptions.

    Sequence case is retained because soft masking can affect downstream tools.
    Gzip inputs are recognised by their bytes, including Nextflow-staged files
    whose names no longer have the original extension.
    """
    with path.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else Path.open
    digest = hashlib.sha256()
    seen: set[str] = set()
    current_id = None
    sequence_length = 0
    with opener(path, "rt", encoding="ascii") as handle:
        for number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None and sequence_length == 0:
                    raise CohortUpdateError(
                        f"Empty FASTA record {current_id!r} in {path}"
                    )
                tokens = line[1:].split()
                if not tokens or tokens[0] in seen:
                    raise CohortUpdateError(
                        f"Empty or duplicate FASTA ID at {path}:{number}"
                    )
                current_id = tokens[0]
                seen.add(current_id)
                sequence_length = 0
                digest.update(b"\x00record\x00" + current_id.encode("ascii") + b"\x00")
            else:
                if current_id is None or re.search(
                    r"[^ACGTURYKMSWBDHVNXacgturykmswbdhvnx.\-]", line
                ):
                    raise CohortUpdateError(
                        f"Invalid genome FASTA sequence at {path}:{number}"
                    )
                digest.update(line.encode("ascii"))
                sequence_length += len(line)
    if current_id is None or sequence_length == 0:
        raise CohortUpdateError(f"Empty genome FASTA or final record: {path}")
    return digest.hexdigest()


def allocate_internal_ids(
    requested: Sequence[str],
    source: dict[str, dict[str, str]],
) -> dict[str, str]:
    """Preserve published IDs and allocate deterministic IDs only for additions."""
    reserved = set()
    for accession, row in source.items():
        internal_id = row["internal_id"]
        if not re.fullmatch(r"[A-Za-z0-9_]+", internal_id) or internal_id in reserved:
            raise CohortUpdateError(
                f"Invalid or duplicate published internal_id for {accession!r}"
            )
        reserved.add(internal_id)
    result = {
        accession: source[accession]["internal_id"]
        for accession in requested
        if accession in source
    }
    additions = [accession for accession in requested if accession not in source]
    candidates = validate_inputs.add_collision_suffixes(additions)
    for accession in sorted(additions):
        candidate = candidates[accession]
        if candidate in reserved:
            base = validate_inputs.sanitise_accession(accession)
            suffix = hashlib.sha1(accession.encode("utf-8")).hexdigest()[:8]
            candidate = f"{base}_{suffix}"
        if candidate in reserved:
            raise CohortUpdateError(
                f"Unable to allocate a unique internal ID for {accession!r}"
            )
        result[accession] = candidate
        reserved.add(candidate)
    return result


def require_files(root: Path, relative_paths: Sequence[str]) -> None:
    """Require published artefacts, including zero-byte documented failure outputs."""
    for relative in relative_paths:
        if not (root / relative).is_file():
            raise CohortUpdateError(f"Missing published artefact: {root / relative}")


def read_exit_code(path: Path) -> str:
    """Require the final recorded exit code, never assume success from file presence."""
    codes = [
        line.split("=", 1)[1]
        for line in path.read_text().splitlines()
        if line.startswith("exit_code=")
    ]
    if not codes or not codes[-1].isdigit():
        raise CohortUpdateError(f"Missing or invalid final exit_code in {path}")
    return codes[-1]


def check_status(
    status: dict[str, str], column: str, derived: str, accession: str
) -> None:
    """Reject published summaries that disagree with the completed audit table."""
    if status[column] != derived:
        raise CohortUpdateError(
            f"Published {column} disagrees with artefacts for {accession!r}: "
            f"{status[column]!r} != {derived!r}"
        )


def validate_sixteen_s(
    root: Path, accession: str, status: dict[str, str]
) -> dict[str, str]:
    """Check the stored 16S classification against its selected sequence."""
    sixteen_s = single_summary(
        root / "16s/16S_status.tsv",
        accession,
        (
            "16S",
            "best_16S_header",
            "best_16S_length",
            "warnings",
        ),
    )
    if sixteen_s["16S"] not in {"Yes", "No", "partial", "NA"}:
        raise CohortUpdateError(f"Invalid published 16S value for {accession!r}")
    derived, _warnings = build_sample_status.derive_barrnap_status(
        accession, {accession: sixteen_s}, True
    )
    check_status(status, "barrnap_status", derived, accession)
    selected = summarise_16s.parse_fasta_records(root / "16s/best_16S.fna")
    if sixteen_s["16S"] in {"Yes", "partial"}:
        if (
            len(selected) != 1
            or selected[0].header != sixteen_s["best_16S_header"]
            or str(len(selected[0].sequence)) != sixteen_s["best_16S_length"]
        ):
            raise CohortUpdateError(
                f"Selected 16S sequence disagrees with its summary for {accession!r}"
            )
    elif selected:
        raise CohortUpdateError(f"Unexpected selected 16S sequence for {accession!r}")
    return sixteen_s


def validate_qc(
    root: Path, sample: dict[str, str], status: dict[str, str]
) -> dict[str, str]:
    """Verify the saved QC call without applying a new genetic-code rule."""
    accession = sample["accession"]
    qc = single_summary(
        root / "checkm2/checkm2_summary.tsv",
        accession,
        (*build_master_table.CHECKM2_COLUMNS, "warnings"),
    )
    if qc["Gcode"] not in {"4", "11", "NA"} or qc["Low_quality"] not in {
        "true",
        "false",
        "NA",
    }:
        raise CohortUpdateError(
            f"Invalid published QC classification for {accession!r}"
        )
    for column in master_table_contract.CHECKM2_COLUMNS:
        value = qc[column]
        if value != "NA" and (
            not value or not math.isfinite(float(value)) or float(value) < 0
        ):
            raise CohortUpdateError(f"Invalid published {column} for {accession!r}")
    qc_statuses = build_sample_status.derive_checkm2_statuses(
        accession, {accession: qc}, True
    )
    for column, value in zip(
        ("checkm2_gcode4_status", "checkm2_gcode11_status", "gcode_status"),
        qc_statuses[:3],
        strict=True,
    ):
        check_status(status, column, value, accession)
    check_status(status, "gcode", qc["Gcode"], accession)
    check_status(status, "low_quality", qc["Low_quality"], accession)
    for code in (4, 11):
        if status[f"checkm2_gcode{code}_status"] == "done":
            report_path = root / f"checkm2_gcode{code}/quality_report.tsv"
            _columns, reports = read_table(report_path, ("Name",))
            if len(reports) != 1 or reports[0]["Name"] != sample["internal_id"]:
                raise CohortUpdateError(
                    f"CheckM2 report identity mismatch: {report_path}"
                )
            report = summarise_checkm2.parse_report(report_path)
            for metric, value in report.metrics.items():
                if (
                    summarise_checkm2.format_metric(value)
                    != qc[f"{metric}_gcode{code}"]
                ):
                    raise CohortUpdateError(
                        f"CheckM2 report disagrees with published summary: {report_path}"
                    )
    return qc


def validate_busco(
    root: Path, accession: str, status: dict[str, str], lineages: Sequence[str]
) -> list[dict[str, str]]:
    """Verify each requested lineage against its published summary."""
    summaries = []
    for lineage in lineages:
        require_files(
            root, (f"busco/{lineage}/short_summary.json", f"busco/{lineage}/busco.log")
        )
        summary = single_summary(
            root / f"busco/{lineage}/busco_summary_{lineage}.tsv",
            accession,
            (f"BUSCO_{lineage}", "busco_status", "warnings"),
        )
        if summary["busco_status"] not in {"done", "failed"}:
            raise CohortUpdateError(
                f"Invalid published BUSCO status for {accession!r}: {lineage}"
            )
        check_status(
            status, f"busco_{lineage}_status", summary["busco_status"], accession
        )
        if summary["busco_status"] == "done":
            parsed = summarise_busco.parse_summary(
                root / f"busco/{lineage}/short_summary.json"
            )
            if parsed != summary[f"BUSCO_{lineage}"]:
                raise CohortUpdateError(
                    f"BUSCO report disagrees with published summary for {accession!r}: {lineage}"
                )
        elif summary[f"BUSCO_{lineage}"] != "NA":
            raise CohortUpdateError(
                f"Failed BUSCO report has a non-NA value for {accession!r}: {lineage}"
            )
        summaries.append(summary)
    return summaries


def validate_annotations(
    root: Path, accession: str, status: dict[str, str], gcode: str
) -> list[dict[str, str]]:
    """Keep documented annotation outcomes and reject contradictory evidence."""
    summaries = []
    codetta = single_summary(
        root / "codetta/codetta_summary.tsv",
        accession,
        build_sample_status.CODETTA_SUMMARY_COLUMNS,
    )
    if codetta["codetta_status"] not in {"done", "failed"}:
        raise CohortUpdateError(f"Invalid published Codetta status for {accession!r}")
    check_status(status, "codetta_status", codetta["codetta_status"], accession)
    summaries.append(codetta)
    if gcode == "NA":
        for tool in ("prokka", "ccfinder", "padloc", "eggnog"):
            check_status(status, f"{tool}_status", "skipped", accession)
    else:
        require_files(
            root,
            (
                "prokka/prokka.gff",
                "prokka/prokka.faa",
                "prokka/prokka.gbk",
                "prokka/prokka.log",
                "ccfinder/result.json",
                "ccfinder/ccfinder.log",
                "ccfinder/ccfinder_contigs.tsv",
                "ccfinder/ccfinder_crisprs.tsv",
                "padloc/padloc.log",
            ),
        )
        ccfinder = single_summary(
            root / "ccfinder/ccfinder_strains.tsv",
            accession,
            (*build_master_table.CRISPR_COLUMNS, "ccfinder_status", "warnings"),
        )
        if ccfinder["ccfinder_status"] not in {"done", "failed"}:
            raise CohortUpdateError(
                f"Invalid published CRISPRCasFinder status for {accession!r}"
            )
        check_status(status, "ccfinder_status", ccfinder["ccfinder_status"], accession)
        if ccfinder["ccfinder_status"] == "done":
            for column in build_master_table.CRISPR_COLUMNS:
                number = float(ccfinder[column])
                if not math.isfinite(number) or number < 0:
                    raise CohortUpdateError(
                        f"Invalid published {column} for {accession!r}"
                    )
            with (root / "ccfinder/result.json").open() as handle:
                json.load(handle)
        summaries.append(ccfinder)
        for tool in ("prokka", "padloc", "eggnog"):
            if tool == "eggnog" and status["eggnog_status"] == "skipped":
                continue
            tool_root = root / tool
            exit_code = read_exit_code(tool_root / f"{tool}.log")
            if tool == "prokka":
                present = all(
                    (tool_root / name).stat().st_size > 0
                    for name in ("prokka.gff", "prokka.faa")
                )
            else:
                if not (tool_root / tool).is_dir():
                    raise CohortUpdateError(
                        f"Missing published result directory: {tool_root / tool}"
                    )
                if tool == "eggnog":
                    require_files(root, ("eggnog/eggnog_annotations.tsv",))
                    present = (tool_root / "eggnog_annotations.tsv").stat().st_size > 0
                else:
                    present = any(
                        path.is_file() for path in (tool_root / tool).iterdir()
                    )
            derived = "done" if exit_code == "0" and present else "failed"
            check_status(status, f"{tool}_status", derived, accession)
    return summaries


def validate_sample_outputs(
    root: Path,
    sample: dict[str, str],
    status: dict[str, str],
    master: dict[str, str],
    lineages: Sequence[str],
) -> str:
    """Validate reusable summaries, outcome evidence and the published genome."""
    accession = sample["accession"]
    if status["validation_status"] != "done":
        raise CohortUpdateError(
            f"Source sample was not successfully validated: {accession!r}"
        )
    if (
        status["internal_id"] != sample["internal_id"]
        or status["is_new"] != sample["is_new"]
    ):
        raise CohortUpdateError(
            f"Published manifest/status identity mismatch for {accession!r}"
        )
    # Published copies must be independent of work directories. Do not follow
    # nested links that could escape the sample tree or create copy cycles.
    for path in root.rglob("*"):
        if path.is_symlink():
            raise CohortUpdateError(f"Published sample artefact is a symlink: {path}")
    require_files(
        root,
        (
            f"staged/{sample['internal_id']}.fasta",
            "barrnap/rrna.gff",
            "barrnap/rrna.fa",
            "barrnap/barrnap.log",
            "16s/best_16S.fna",
            "16s/16S_status.tsv",
            "checkm2/checkm2_summary.tsv",
            "checkm2_gcode4/quality_report.tsv",
            "checkm2_gcode4/checkm2.log",
            "checkm2_gcode11/quality_report.tsv",
            "checkm2_gcode11/checkm2.log",
            "codetta/codetta_summary.tsv",
            "codetta/codetta.log",
        ),
    )
    sixteen_s = validate_sixteen_s(root, accession, status)
    qc = validate_qc(root, sample, status)
    master_summaries = [qc, sixteen_s]
    master_summaries.extend(validate_busco(root, accession, status, lineages))
    master_summaries.extend(validate_annotations(root, accession, status, qc["Gcode"]))
    derived_columns = set(master_table_contract.build_append_columns(lineages))
    for summary in master_summaries:
        for column in derived_columns & summary.keys():
            if master[column] != summary[column]:
                raise CohortUpdateError(
                    f"Published master table disagrees with {column} summary for {accession!r}"
                )
    return genome_fingerprint(root / f"staged/{sample['internal_id']}.fasta")


def write_table(
    path: Path, header: Sequence[str], rows: Sequence[dict[str, str]]
) -> None:
    """Write deterministic LF-terminated update tables."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=header,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row[column] for column in header})


def load_source_tables(
    source: Path, lineages: Sequence[str], requested: set[str]
) -> tuple[
    dict[str, dict[str, str]],
    dict[str, dict[str, str]],
    dict[str, dict[str, str]],
    list[dict[str, str]],
]:
    """Load a completed source cohort and validate the reporting contracts."""
    source_manifest = source / "tables/validated_samples.tsv"
    _header, source_rows = read_table(
        source_manifest, (*validate_inputs.REQUIRED_SAMPLE_COLUMNS, "internal_id")
    )
    if not source_rows:
        raise CohortUpdateError("The published source manifest contains no samples.")
    source_index = index_rows(source_rows, "accession", source_manifest)
    for row in source_rows:
        if row["is_new"] not in {"true", "false"}:
            raise CohortUpdateError(
                f"Invalid published is_new for {row['accession']!r}"
            )
    status_path = source / "tables/sample_status.tsv"
    status_header, status_rows = read_table(status_path)
    source_lineages = (
        master_table_contract.extract_busco_lineages_from_sample_status_columns(
            status_header
        )
    )
    if requested & source_index.keys() and not set(lineages) <= set(source_lineages):
        raise CohortUpdateError(
            "The update requests BUSCO lineages absent from the published source."
        )
    statuses = index_rows(status_rows, "accession", status_path)
    master_path = source / "tables/master_table.tsv"
    master_header, master_rows = read_table(
        master_path, master_table_contract.build_append_columns(source_lineages)
    )
    master_key = validate_inputs.detect_metadata_key_column(master_header)
    masters = index_rows(master_rows, master_key, master_path)
    if set(statuses) != set(source_index) or set(masters) != set(source_index):
        raise CohortUpdateError(
            "Published validated, master and status tables must contain the same accessions."
        )
    versions_path = source / "tables/tool_and_db_versions.tsv"
    versions_header, versions = read_table(
        versions_path, collect_versions.OUTPUT_COLUMNS
    )
    if versions_header != list(collect_versions.OUTPUT_COLUMNS) or not versions:
        raise CohortUpdateError(
            f"Invalid or empty published provenance report: {versions_path}"
        )
    return source_index, statuses, masters, versions


def write_update_tables(
    args: argparse.Namespace,
    *,
    sample_header: list[str],
    samples: list[dict[str, str]],
    added: list[dict[str, str]],
    reused: list[dict[str, str]],
    audit: list[dict[str, str]],
    versions: list[dict[str, str]],
    identity: dict[str, object],
) -> None:
    """Publish current validation state, sample routing and labelled provenance."""
    internal_ids = {row["accession"]: row["internal_id"] for row in samples}
    map_header, mapping = read_table(
        args.accession_map, validate_inputs.ACCESSION_MAP_COLUMNS
    )
    initial_header, initial = read_table(
        args.initial_status, ("accession", "internal_id")
    )
    if initial_header != master_table_contract.build_sample_status_columns(
        args.busco_lineage
    ):
        raise CohortUpdateError(
            "Current validation status columns do not match the requested BUSCO lineages."
        )
    warning_header, warnings = read_table(
        args.validation_warnings, validate_inputs.VALIDATION_WARNING_COLUMNS
    )
    for rows, path in ((mapping, args.accession_map), (initial, args.initial_status)):
        if set(index_rows(rows, "accession", path)) != set(internal_ids):
            raise CohortUpdateError(
                f"Current validation table has inconsistent sample membership: {path}"
            )
    for row in [*mapping, *initial]:
        row["internal_id"] = internal_ids[row["accession"]]
    # Validation calculated IDs without the source mapping. Replace only its
    # collision warnings; all other current-manifest warnings remain current.
    warnings = [
        row
        for row in warnings
        if row["warning_code"] != "internal_id_collision_resolved"
    ]
    for sample in samples:
        accession = sample["accession"]
        if sample["internal_id"] != validate_inputs.sanitise_accession(accession):
            warnings.append(
                {
                    "accession": accession,
                    "warning_code": "internal_id_collision_resolved",
                    "message": "Preserved or assigned a deterministic internal ID to avoid an accession collision.",
                }
            )
    mapping_index = index_rows(mapping, "accession", args.accession_map)
    records = [
        validate_inputs.SampleRecord(
            {key: value for key, value in sample.items() if key != "internal_id"},
            sample["internal_id"],
            mapping_index[sample["accession"]]["metadata_present"] == "true",
        )
        for sample in samples
    ]
    initial = validate_inputs.build_initial_sample_status_rows(
        records,
        [validate_inputs.ValidationWarning(**row) for row in warnings],
        initial_header,
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    for name, header, rows in (
        ("validated_samples.tsv", sample_header, samples),
        ("accession_map.tsv", map_header, mapping),
        ("sample_status.tsv", initial_header, initial),
        ("validation_warnings.tsv", warning_header, warnings),
        ("new_samples.tsv", sample_header, added),
        (
            "reused_samples.tsv",
            [
                *validate_inputs.REQUIRED_SAMPLE_COLUMNS,
                "internal_id",
                "source_gcode",
                "source_eggnog_status",
            ],
            reused,
        ),
        ("cohort_update.tsv", AUDIT_COLUMNS, audit),
    ):
        write_table(args.outdir / name, header, rows)
    for row in versions:
        row["notes"] = f"reused source {identity['source_outdir']}: {row['notes']}"
    write_table(
        args.outdir / "inherited_versions.tsv",
        collect_versions.OUTPUT_COLUMNS,
        versions,
    )
    (args.outdir / "cohort_update_run.json").write_text(
        json.dumps(identity, indent=2, sort_keys=True) + "\n"
    )


def run_prepare(args: argparse.Namespace) -> None:
    """Preflight every retained sample before exposing additions to the workflow."""
    source = args.source_results.resolve(strict=True)
    if any(
        not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", lineage)
        for lineage in args.busco_lineage
    ):
        raise CohortUpdateError("BUSCO lineages must be simple directory names.")
    destination = args.destination.resolve()
    if (
        source == destination
        or source in destination.parents
        or destination in source.parents
    ):
        raise CohortUpdateError(
            "--update_from and --outdir must be separate, non-overlapping directories."
        )
    sample_header, samples = read_table(
        args.validated_samples,
        (*validate_inputs.REQUIRED_SAMPLE_COLUMNS, "internal_id"),
    )
    if not samples:
        raise CohortUpdateError(
            "A cohort update requires at least one requested sample."
        )
    requested = index_rows(samples, "accession", args.validated_samples)
    source_index, statuses, masters, versions = load_source_tables(
        source, args.busco_lineage, set(requested)
    )
    source_manifest = source / "tables/validated_samples.tsv"
    versions_path = source / "tables/tool_and_db_versions.tsv"
    internal_ids = allocate_internal_ids(list(requested), source_index)
    candidates = sorted(args.genome_inputs.glob("genome*"))
    if len(candidates) != len(samples):
        raise CohortUpdateError(
            "Staged candidate-genome count does not match the requested manifest."
        )
    manifest_hash = file_sha256(source_manifest)
    versions_hash = file_sha256(versions_path)
    audit = []
    reused = []
    added = []
    for sample, candidate in zip(samples, candidates, strict=True):
        accession = sample["accession"]
        fingerprint = genome_fingerprint(candidate)
        if accession in source_index:
            source_fingerprint = validate_sample_outputs(
                source / "samples" / accession,
                source_index[accession],
                statuses[accession],
                masters[accession],
                args.busco_lineage,
            )
            if source_fingerprint != fingerprint:
                raise CohortUpdateError(
                    f"Genome content changed for retained accession {accession!r}; run a full analysis or assign a new accession."
                )
            reused.append(
                {
                    **sample,
                    "source_gcode": statuses[accession]["gcode"],
                    "source_eggnog_status": statuses[accession]["eggnog_status"],
                }
            )
            action = "reused"
        else:
            added.append(sample)
            action = "added"
        sample["internal_id"] = internal_ids[accession]
        audit.append(
            dict(
                zip(
                    AUDIT_COLUMNS,
                    (
                        accession,
                        action,
                        internal_ids[accession],
                        fingerprint,
                        str(source),
                        manifest_hash,
                        versions_hash,
                    ),
                    strict=True,
                )
            )
        )
    for accession in sorted(set(source_index) - set(requested)):
        audit.append(
            dict(
                zip(
                    AUDIT_COLUMNS,
                    (
                        accession,
                        "removed",
                        source_index[accession]["internal_id"],
                        "NA",
                        str(source),
                        manifest_hash,
                        versions_hash,
                    ),
                    strict=True,
                )
            )
        )
    # Reused rows were copied before remapping. Retained IDs always come from
    # the source, including collisions newly introduced into this cohort.
    for sample in reused:
        sample["internal_id"] = internal_ids[sample["accession"]]
    identity = {
        "source_outdir": str(source),
        "source_manifest_sha256": manifest_hash,
        "source_versions_sha256": versions_hash,
        "metadata_sha256": file_sha256(args.metadata),
        "samples": samples,
        "genomes": {
            row["accession"]: row["genome_sha256"]
            for row in audit
            if row["action"] != "removed"
        },
        "settings": json.loads(args.settings.read_text()),
    }
    if (
        args.previous_update
        and json.loads(args.previous_update.read_text()) != identity
    ):
        raise CohortUpdateError(
            "The existing update output belongs to different inputs or settings; choose a fresh --outdir."
        )
    write_update_tables(
        args,
        sample_header=sample_header,
        samples=samples,
        added=added,
        reused=reused,
        audit=audit,
        versions=versions,
        identity=identity,
    )
    LOGGER.info(
        "Cohort update validated: %d reused, %d added, %d removed.",
        len(reused),
        len(added),
        len(source_index.keys() - requested.keys()),
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Run the update preflight CLI."""
    parser = argparse.ArgumentParser(description=__doc__)
    for option in (
        "source-results",
        "destination",
        "validated-samples",
        "accession-map",
        "initial-status",
        "validation-warnings",
        "genome-inputs",
        "metadata",
        "settings",
        "outdir",
    ):
        parser.add_argument(f"--{option}", type=Path, required=True)
    parser.add_argument("--busco-lineage", action="append", required=True)
    parser.add_argument("--previous-update", type=Path)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    try:
        run_prepare(args)
    except (
        ValueError,
        OSError,
        UnicodeError,
        validate_inputs.ValidationError,
    ) as error:
        LOGGER.error("Cohort update preflight failed: %s", error)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
