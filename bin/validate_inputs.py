#!/usr/bin/env python3
"""Validate pipeline inputs and create canonical sample-mapping tables."""

from __future__ import annotations

import argparse
import csv
import hashlib
import logging
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import master_table_contract


LOGGER = logging.getLogger(__name__)

REQUIRED_SAMPLE_COLUMNS = (
    "accession",
    "is_new",
    "assembly_level",
    "genome_fasta",
)
RESERVED_SAMPLE_COLUMNS = {"internal_id"}
ACCESSION_MAP_COLUMNS = (
    "accession",
    "internal_id",
    "is_new",
    "assembly_level",
    "genome_fasta",
    "metadata_present",
)
VALIDATION_WARNING_COLUMNS = ("accession", "warning_code", "message")
MISSING_VALUE_TOKENS = {"", "na", "n/a", "null", "none"}
TRUE_TOKENS = {"true", "t", "yes", "y", "1"}
FALSE_TOKENS = {"false", "f", "no", "n", "0"}
ASSEMBLY_LEVEL_MAP = {
    "complete genome": "Complete Genome",
    "complete genomes": "Complete Genome",
    "chromosome": "Chromosome",
    "chromosomes": "Chromosome",
    "scaffold": "Scaffold",
    "scaffolds": "Scaffold",
    "contig": "Contig",
    "contigs": "Contig",
}


class ValidationError(RuntimeError):
    """Raised when the sample manifest or metadata table is invalid."""


@dataclass(frozen=True)
class ValidationWarning:
    """Represent one non-fatal validation warning."""

    accession: str
    warning_code: str
    message: str


@dataclass(frozen=True)
class SampleRecord:
    """Store one validated sample-manifest row."""

    values: dict[str, str]
    internal_id: str
    metadata_present: bool


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate the sample manifest and metadata table."
    )
    parser.add_argument(
        "--sample-csv",
        required=True,
        type=Path,
        help="Path to the sample manifest CSV.",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="Path to the metadata table (CSV or TSV).",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Directory for validated output tables.",
    )
    parser.add_argument(
        "--sample-status-columns",
        type=Path,
        help="Optional path to the ordered sample-status column asset override.",
    )
    parser.add_argument(
        "--busco-lineage",
        action="append",
        help="Configured BUSCO lineage name. May be supplied multiple times.",
    )
    parser.add_argument(
        "--defer-genome-fasta-check",
        action="store_true",
        help=(
            "Skip immediate genome_fasta existence checks and defer them to "
            "the downstream workflow staging step."
        ),
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
    )


def is_missing(value: str | None) -> bool:
    """Return True when a field should be treated as missing."""
    if value is None:
        return True
    return value.strip().lower() in MISSING_VALUE_TOKENS


def normalise_optional_value(value: str | None) -> str:
    """Normalise optional fields to a stripped string or `NA`."""
    if is_missing(value):
        return "NA"
    return value.strip()


def normalise_assembly_level(value: str | None, accession: str) -> str:
    """Normalise a sample-manifest assembly level to the canonical ANI token."""
    if is_missing(value):
        return "NA"
    token = value.strip().casefold()
    if token in ASSEMBLY_LEVEL_MAP:
        return ASSEMBLY_LEVEL_MAP[token]
    allowed_values = ", ".join(sorted(set(ASSEMBLY_LEVEL_MAP.values())))
    raise ValidationError(
        f"Sample {accession!r} has invalid assembly_level {value!r}. "
        f"Expected one of: {allowed_values}."
    )


def normalise_boolean(value: str) -> str:
    """Normalise boolean-like values to `true` or `false`."""
    token = value.strip().lower()
    if token in TRUE_TOKENS:
        return "true"
    if token in FALSE_TOKENS:
        return "false"
    raise ValidationError(f"Invalid is_new value: {value!r}")


def ensure_path_safe_accession(accession: str) -> None:
    """Raise when an accession cannot be used safely as a published folder name."""
    if not accession:
        raise ValidationError("Encountered an empty accession.")
    if accession in {".", ".."}:
        raise ValidationError(f"Accession {accession!r} is not filesystem-safe.")
    if "/" in accession or "\\" in accession:
        raise ValidationError(
            f"Accession {accession!r} contains a path separator and is not allowed."
        )
    if any(ord(char) < 32 for char in accession):
        raise ValidationError(
            f"Accession {accession!r} contains control characters and is not allowed."
        )


def sanitise_accession(accession: str) -> str:
    """Convert an accession into a tool-safe internal identifier base."""
    characters: list[str] = []
    previous_was_underscore = False
    for char in accession:
        if char.isascii() and (char.isalnum() or char == "_"):
            if char == "_":
                if previous_was_underscore:
                    continue
                characters.append(char)
                previous_was_underscore = True
                continue
            characters.append(char)
            previous_was_underscore = False
            continue
        if not previous_was_underscore:
            characters.append("_")
            previous_was_underscore = True
    return "".join(characters).strip("_")


def add_collision_suffixes(accessions: Sequence[str]) -> dict[str, str]:
    """Build a deterministic original-accession to internal-ID mapping."""
    grouped_accessions: dict[str, list[str]] = defaultdict(list)
    for accession in accessions:
        internal_base = sanitise_accession(accession)
        if not internal_base:
            raise ValidationError(
                f"Accession {accession!r} sanitises to an empty internal ID."
            )
        grouped_accessions[internal_base].append(accession)

    accession_map: dict[str, str] = {}
    for internal_base, grouped in grouped_accessions.items():
        if len(grouped) == 1:
            accession_map[grouped[0]] = internal_base
            continue
        for accession in sorted(grouped):
            suffix = hashlib.sha1(accession.encode("utf-8")).hexdigest()[:8]
            accession_map[accession] = f"{internal_base}_{suffix}"

    if len(set(accession_map.values())) != len(accession_map):
        raise ValidationError("Unable to resolve internal-ID collisions deterministically.")
    return accession_map


def sniff_delimiter(path: Path) -> str:
    """Guess whether a delimited table is comma- or tab-separated."""
    sample = path.read_text(encoding="utf-8").splitlines()
    sample_text = "\n".join(sample[:10])
    if not sample_text.strip():
        raise ValidationError(f"Input table is empty: {path}")
    try:
        dialect = csv.Sniffer().sniff(sample_text, delimiters=",\t")
    except csv.Error:
        return "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    return dialect.delimiter


def read_delimited_table(path: Path, delimiter: str | None = None) -> tuple[list[str], list[dict[str, str]]]:
    """Read a CSV/TSV file into a header list and row dictionaries."""
    if not path.is_file():
        raise ValidationError(f"Missing input file: {path}")

    table_delimiter = delimiter or sniff_delimiter(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter=table_delimiter)
        rows = [row for row in reader if any(cell.strip() for cell in row)]

    if not rows:
        raise ValidationError(f"Input table is empty: {path}")

    header = [column.strip() for column in rows[0]]
    if any(not column for column in header):
        raise ValidationError(f"Input table has an empty header name: {path}")
    if len(set(header)) != len(header):
        raise ValidationError(f"Input table contains duplicate header names: {path}")

    records: list[dict[str, str]] = []
    for line_number, values in enumerate(rows[1:], start=2):
        if len(values) != len(header):
            raise ValidationError(
                f"Input table {path} has {len(values)} columns on line {line_number}, "
                f"expected {len(header)}."
            )
        records.append(
            {
                column: value.strip()
                for column, value in zip(header, values, strict=True)
            }
        )
    return header, records


def detect_metadata_key_column(header: Sequence[str]) -> str:
    """Return the metadata accession key column name."""
    has_lower = "accession" in header
    has_upper = "Accession" in header
    if has_lower and has_upper:
        raise ValidationError(
            "Metadata table cannot contain both 'accession' and 'Accession' columns."
        )
    if has_upper:
        return "Accession"
    if has_lower:
        return "accession"
    raise ValidationError("Metadata table must contain 'Accession' or 'accession'.")


def validate_sample_header(header: Sequence[str]) -> None:
    """Validate the sample-manifest schema."""
    missing = [column for column in REQUIRED_SAMPLE_COLUMNS if column not in header]
    if missing:
        raise ValidationError(
            f"Sample manifest is missing required columns: {', '.join(missing)}"
        )
    for column in header:
        if column != column.lower():
            raise ValidationError(
                f"Sample manifest column {column!r} must be lower-case."
            )
        if column in RESERVED_SAMPLE_COLUMNS:
            raise ValidationError(
                f"Sample manifest column {column!r} is reserved for derived outputs."
            )


def build_metadata_index(metadata_rows: Sequence[dict[str, str]], key_column: str) -> dict[str, dict[str, str]]:
    """Build a unique accession-keyed metadata index."""
    metadata_index: dict[str, dict[str, str]] = {}
    for row in metadata_rows:
        accession = row[key_column].strip()
        if is_missing(accession):
            raise ValidationError("Metadata table contains a row with a missing accession.")
        if accession in metadata_index:
            raise ValidationError(
                f"Metadata table contains duplicate accession rows for {accession!r}."
            )
        metadata_index[accession] = row
    return metadata_index


def validate_samples(
    sample_header: Sequence[str],
    sample_rows: Sequence[dict[str, str]],
    metadata_index: dict[str, dict[str, str]],
    *,
    defer_genome_fasta_check: bool = False,
) -> tuple[list[SampleRecord], list[ValidationWarning]]:
    """Validate sample rows against manifest and metadata rules."""
    seen_accessions: set[str] = set()
    accessions: list[str] = []

    for row in sample_rows:
        accession = row["accession"].strip()
        if is_missing(accession):
            raise ValidationError("Sample manifest contains a row with a missing accession.")
        ensure_path_safe_accession(accession)
        if accession in seen_accessions:
            raise ValidationError(f"Sample manifest contains duplicate accession {accession!r}.")
        seen_accessions.add(accession)
        accessions.append(accession)

    if not accessions:
        raise ValidationError("Sample manifest must contain at least one sample.")

    accession_map = add_collision_suffixes(accessions)
    warnings: list[ValidationWarning] = []
    collision_groups: dict[str, list[str]] = defaultdict(list)
    for accession in accessions:
        collision_groups[sanitise_accession(accession)].append(accession)
    for grouped in collision_groups.values():
        if len(grouped) <= 1:
            continue
        for accession in grouped:
            warnings.append(
                ValidationWarning(
                    accession=accession,
                    warning_code="internal_id_collision_resolved",
                    message=(
                        "Sanitised accession collided with another sample and received "
                        "a deterministic suffix."
                    ),
                )
            )

    validated_records: list[SampleRecord] = []
    for row in sample_rows:
        accession = row["accession"].strip()
        is_new = normalise_boolean(row["is_new"])
        assembly_level = normalise_assembly_level(row["assembly_level"], accession)
        if is_new == "true" and assembly_level == "NA":
            raise ValidationError(
                f"Sample {accession!r} has is_new=true but no assembly_level."
            )

        genome_value = row["genome_fasta"].strip()
        if is_missing(genome_value):
            raise ValidationError(
                f"Sample {accession!r} is missing a genome_fasta path."
            )
        genome_path = Path(genome_value).expanduser()
        if not defer_genome_fasta_check and not genome_path.exists():
            raise ValidationError(
                f"Genome FASTA for sample {accession!r} does not exist: {genome_value}"
            )

        metadata_present = accession in metadata_index
        if is_new == "false" and not metadata_present:
            raise ValidationError(
                f"Sample {accession!r} has is_new=false but no metadata row exists."
            )
        if is_new == "true" and not metadata_present:
            warnings.append(
                ValidationWarning(
                    accession=accession,
                    warning_code="missing_metadata_for_new_sample",
                    message=(
                        "New sample has no metadata row; downstream metadata fields "
                        "must be filled with NA."
                    ),
                )
            )

        validated_values = {
            column: normalise_optional_value(value)
            for column, value in row.items()
        }
        validated_values["accession"] = accession
        validated_values["is_new"] = is_new
        validated_values["assembly_level"] = assembly_level
        validated_values["genome_fasta"] = str(genome_path.resolve())

        validated_records.append(
            SampleRecord(
                values=validated_values,
                internal_id=accession_map[accession],
                metadata_present=metadata_present,
            )
        )

    return validated_records, warnings


def write_tsv(path: Path, header: Sequence[str], rows: Sequence[dict[str, str]]) -> None:
    """Write a sequence of dictionaries to a TSV file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(header), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_validated_sample_rows(
    sample_header: Sequence[str], records: Sequence[SampleRecord]
) -> list[dict[str, str]]:
    """Convert validated records into validated-manifest rows."""
    rows: list[dict[str, str]] = []
    for record in records:
        row = {column: record.values[column] for column in sample_header}
        row["internal_id"] = record.internal_id
        rows.append(row)
    return rows


def build_accession_map_rows(records: Sequence[SampleRecord]) -> list[dict[str, str]]:
    """Convert validated records into accession-map rows."""
    rows: list[dict[str, str]] = []
    for record in records:
        rows.append(
            {
                "accession": record.values["accession"],
                "internal_id": record.internal_id,
                "is_new": record.values["is_new"],
                "assembly_level": record.values["assembly_level"],
                "genome_fasta": record.values["genome_fasta"],
                "metadata_present": "true" if record.metadata_present else "false",
            }
        )
    return rows


def build_warning_rows(warnings: Sequence[ValidationWarning]) -> list[dict[str, str]]:
    """Convert warning objects into TSV rows."""
    return [
        {
            "accession": warning.accession,
            "warning_code": warning.warning_code,
            "message": warning.message,
        }
        for warning in warnings
    ]


def resolve_sample_status_columns(
    path: Path | None = None,
    busco_lineages: Sequence[str] | None = None,
) -> list[str]:
    """Resolve the sample-status contract from an override, lineages, or defaults."""
    try:
        sample_status_columns, _contract_lineages = (
            master_table_contract.resolve_sample_status_columns(
                path=path,
                busco_lineages=busco_lineages,
            )
        )
    except ValueError as error:
        raise ValidationError(str(error)) from error
    return sample_status_columns


def build_initial_sample_status_row(
    record: SampleRecord,
    warnings: Sequence[ValidationWarning],
    sample_status_columns: Sequence[str],
) -> dict[str, str]:
    """Build the initial sample-status row emitted immediately after validation."""
    row: dict[str, str] = {}
    for column in sample_status_columns:
        if column.endswith("_status") or column == "ani_included":
            row[column] = "na"
        elif column in {"gcode", "low_quality"}:
            row[column] = "NA"
        else:
            row[column] = ""

    row["accession"] = record.values["accession"]
    row["internal_id"] = record.internal_id
    row["is_new"] = record.values["is_new"]
    row["validation_status"] = "done"
    row["warnings"] = ";".join(
        dict.fromkeys(warning.warning_code for warning in warnings)
    )
    row["notes"] = "; ".join(
        dict.fromkeys(warning.message for warning in warnings)
    )
    return row


def build_initial_sample_status_rows(
    records: Sequence[SampleRecord],
    warnings: Sequence[ValidationWarning],
    sample_status_columns: Sequence[str],
) -> list[dict[str, str]]:
    """Build one initial sample-status row per validated sample."""
    warnings_by_accession: dict[str, list[ValidationWarning]] = defaultdict(list)
    for warning in warnings:
        warnings_by_accession[warning.accession].append(warning)

    return [
        build_initial_sample_status_row(
            record=record,
            warnings=warnings_by_accession.get(record.values["accession"], []),
            sample_status_columns=sample_status_columns,
        )
        for record in records
    ]


def run_validation(
    sample_csv: Path,
    metadata: Path,
    outdir: Path,
    sample_status_columns_path: Path | None,
    *,
    busco_lineages: Sequence[str] | None = None,
    defer_genome_fasta_check: bool = False,
) -> None:
    """Validate inputs and write the expected downstream TSV outputs."""
    sample_header, sample_rows = read_delimited_table(sample_csv, delimiter=",")
    validate_sample_header(sample_header)

    metadata_header, metadata_rows = read_delimited_table(metadata)
    metadata_key_column = detect_metadata_key_column(metadata_header)
    metadata_index = build_metadata_index(metadata_rows, metadata_key_column)

    validated_records, warnings = validate_samples(
        sample_header=sample_header,
        sample_rows=sample_rows,
        metadata_index=metadata_index,
        defer_genome_fasta_check=defer_genome_fasta_check,
    )
    sample_status_columns = resolve_sample_status_columns(
        sample_status_columns_path,
        busco_lineages=busco_lineages,
    )

    write_tsv(
        outdir / "validated_samples.tsv",
        [*sample_header, "internal_id"],
        build_validated_sample_rows(sample_header, validated_records),
    )
    write_tsv(
        outdir / "accession_map.tsv",
        ACCESSION_MAP_COLUMNS,
        build_accession_map_rows(validated_records),
    )
    write_tsv(
        outdir / "validation_warnings.tsv",
        VALIDATION_WARNING_COLUMNS,
        build_warning_rows(warnings),
    )
    write_tsv(
        outdir / "sample_status.tsv",
        sample_status_columns,
        build_initial_sample_status_rows(
            records=validated_records,
            warnings=warnings,
            sample_status_columns=sample_status_columns,
        ),
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Run the validation CLI."""
    args = parse_args(argv)
    configure_logging()
    try:
        run_validation(
            args.sample_csv,
            args.metadata,
            args.outdir,
            args.sample_status_columns,
            busco_lineages=args.busco_lineage,
            defer_genome_fasta_check=args.defer_genome_fasta_check,
        )
    except ValidationError as error:
        LOGGER.error(str(error))
        return 1

    LOGGER.info("Validated input manifest and metadata successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
