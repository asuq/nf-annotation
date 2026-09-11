#!/usr/bin/env python3
"""Assemble the final master table from metadata and derived sample summaries."""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import master_table_contract


LOGGER = logging.getLogger(__name__)

VALIDATED_SAMPLE_REQUIRED_COLUMNS = (
    "accession",
    "is_new",
    "assembly_level",
    "genome_fasta",
    "internal_id",
)
VALIDATED_SAMPLE_SUPPLEMENTAL_EXCLUSIONS = (
    "accession",
    "is_new",
    "genome_fasta",
    "internal_id",
)
TAXONOMY_COLUMNS = tuple(master_table_contract.TAXONOMY_COLUMNS)
CHECKM2_COLUMNS = (
    *master_table_contract.CHECKM2_COLUMNS,
    "Gcode",
    *master_table_contract.GCODE_PROVENANCE_COLUMNS,
    "Low_quality",
)
ASSEMBLY_DERIVED_COLUMN_MAP = {"GC_Content": "gc_content"}
CODETTA_COLUMNS = tuple(master_table_contract.CODETTA_COLUMNS)
SIXTEEN_S_COLUMNS = ("16S",)
CRISPR_COLUMNS = tuple(master_table_contract.CRISPR_COLUMNS)
ANI_COLUMNS = tuple(master_table_contract.ANI_COLUMNS)
MISSING_VALUE_TOKENS = {"", "na", "n/a", "null", "none"}
ASSEMBLY_STATS_COLUMNS = ("n50", "scaffolds", "genome_size", "gc_content")
ASSEMBLY_METADATA_COLUMN_MAP = {
    "N50": "n50",
    "Scaffolds": "scaffolds",
    "Genome_Size": "genome_size",
}


@dataclass(frozen=True)
class BuscoLineageSummary:
    """Store one BUSCO lineage summary row for one accession."""

    value: str
    status: str
    warnings: str


class MasterTableError(RuntimeError):
    """Raised when the final master table cannot be built unambiguously."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build the final one-row-per-sample master table."
    )
    parser.add_argument(
        "--validated-samples",
        required=True,
        type=Path,
        help="Path to the validated sample manifest TSV.",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="Path to the metadata table (CSV or TSV).",
    )
    parser.add_argument(
        "--append-columns",
        type=Path,
        help="Optional path to the ordered derived-column contract override.",
    )
    parser.add_argument(
        "--busco-lineage",
        action="append",
        help="Configured BUSCO lineage name. May be supplied multiple times.",
    )
    parser.add_argument(
        "--taxonomy",
        type=Path,
        help="Optional taxonomy-expansion TSV keyed by Tax_ID.",
    )
    parser.add_argument(
        "--checkm2",
        type=Path,
        help="Optional per-sample CheckM2 summary TSV.",
    )
    parser.add_argument(
        "--16s-status",
        dest="sixteen_s_status",
        type=Path,
        help="Optional per-sample 16S status TSV.",
    )
    parser.add_argument(
        "--busco",
        action="append",
        default=[],
        type=Path,
        help="Optional per-lineage BUSCO summary TSV. May be supplied multiple times.",
    )
    parser.add_argument(
        "--codetta-summary",
        dest="codetta_summary",
        type=Path,
        help="Optional per-sample Codetta summary TSV.",
    )
    parser.add_argument(
        "--ccfinder-strains",
        type=Path,
        help="Optional CRISPRCasFinder strain summary TSV.",
    )
    parser.add_argument(
        "--ani",
        type=Path,
        help="Optional ANI summary TSV from select_ani_representatives.py.",
    )
    parser.add_argument(
        "--assembly-stats",
        type=Path,
        help="Optional in-house assembly stats TSV keyed by accession.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Path to the final master_table.tsv output.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def normalise_key(key: str) -> str:
    """Normalise a column name for case-insensitive comparisons."""
    return "".join(character for character in key.lower() if character.isalnum())


def is_missing(value: str | None) -> bool:
    """Return True when a scalar should be treated as missing."""
    if value is None:
        return True
    return value.strip().lower() in MISSING_VALUE_TOKENS


def split_tokens(value: str | None) -> list[str]:
    """Split a semicolon-delimited field into stripped non-empty tokens."""
    if value is None:
        return []
    return [token.strip() for token in value.split(";") if token.strip()]


def join_tokens(tokens: Sequence[str]) -> str:
    """Join tokens into a stable semicolon-delimited string without duplicates."""
    return ";".join(dict.fromkeys(token for token in tokens if token))


def join_notes(notes: Sequence[str]) -> str:
    """Join note strings into a stable semicolon-delimited sentence list."""
    return "; ".join(dict.fromkeys(note for note in notes if note))


def sniff_delimiter(path: Path) -> str:
    """Guess whether a table is comma- or tab-separated."""
    sample = path.read_text(encoding="utf-8").splitlines()
    sample_text = "\n".join(sample[:10])
    if not sample_text.strip():
        raise MasterTableError(f"Input table is empty: {path}")
    try:
        dialect = csv.Sniffer().sniff(sample_text, delimiters=",\t")
    except csv.Error:
        return "\t" if path.suffix.lower() in {".tsv", ".tab"} else ","
    return dialect.delimiter


def read_table(path: Path, delimiter: str | None = None) -> tuple[list[str], list[dict[str, str]]]:
    """Read a delimited table into a header and row dictionaries."""
    if not path.is_file():
        raise MasterTableError(f"Missing input table: {path}")
    table_delimiter = delimiter or sniff_delimiter(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter=table_delimiter)
        rows = [row for row in reader if any(cell.strip() for cell in row)]

    if not rows:
        raise MasterTableError(f"Input table is empty: {path}")

    header = [column.strip() for column in rows[0]]
    if any(not column for column in header):
        raise MasterTableError(f"Input table has an empty header name: {path}")
    if len(set(header)) != len(header):
        raise MasterTableError(f"Input table contains duplicate headers: {path}")

    records: list[dict[str, str]] = []
    for line_number, values in enumerate(rows[1:], start=2):
        if len(values) != len(header):
            raise MasterTableError(
                f"Input table {path} has {len(values)} columns on line {line_number}, "
                f"expected {len(header)}."
            )
        records.append(
            {column: value.strip() for column, value in zip(header, values, strict=True)}
        )
    return header, records


def detect_key_column(header: Sequence[str], aliases: Sequence[str]) -> str:
    """Detect a join key column from a list of aliases."""
    matching_columns = [
        column
        for column in header
        if normalise_key(column) in {normalise_key(alias) for alias in aliases}
    ]
    if len(matching_columns) > 1:
        raise MasterTableError(
            f"Input table contains multiple possible key columns: {matching_columns}"
        )
    if not matching_columns:
        raise MasterTableError(
            f"Input table is missing a required key column. Expected one of: {aliases}"
        )
    return matching_columns[0]


def build_index(
    rows: Sequence[dict[str, str]],
    key_column: str,
    table_name: str,
) -> dict[str, dict[str, str]]:
    """Build a unique keyed row index and fail on duplicates."""
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        key = row.get(key_column, "").strip()
        if not key:
            raise MasterTableError(f"{table_name} contains an empty key value.")
        if key in indexed:
            raise MasterTableError(f"{table_name} contains duplicate key {key!r}.")
        indexed[key] = row
    return indexed


def ensure_index_keys_are_valid(
    index: dict[str, dict[str, str]],
    validated_accessions: set[str],
    table_name: str,
) -> None:
    """Raise when an accession-keyed table contains rows outside the validated set."""
    unexpected = sorted(set(index) - validated_accessions)
    if unexpected:
        raise MasterTableError(
            f"{table_name} contains accession values not present in validated_samples.tsv: "
            + ", ".join(unexpected)
        )


def require_columns(
    header: Sequence[str],
    columns: Sequence[str],
    table_name: str,
) -> None:
    """Ensure a table contains a required set of columns."""
    missing = [column for column in columns if column not in header]
    if missing:
        raise MasterTableError(
            f"{table_name} is missing required columns: {', '.join(missing)}"
        )


def load_validated_samples(path: Path) -> list[dict[str, str]]:
    """Load and validate the canonical validated sample table."""
    header, rows = read_table(path, delimiter="\t")
    require_columns(header, VALIDATED_SAMPLE_REQUIRED_COLUMNS, "validated_samples.tsv")
    samples = build_index(rows, "accession", "validated_samples.tsv")
    internal_ids = [row["internal_id"] for row in samples.values()]
    if len(set(internal_ids)) != len(internal_ids):
        raise MasterTableError("validated_samples.tsv contains duplicate internal_id values.")
    return list(samples.values())


def load_sample_status_columns(
    path: Path = master_table_contract.DEFAULT_SAMPLE_STATUS_COLUMNS_ASSET,
) -> list[str]:
    """Load the maintained sample-status column order asset."""
    if not path.is_file():
        raise MasterTableError(f"Missing sample-status column asset: {path}")
    columns = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not columns:
        raise MasterTableError(f"Sample-status column asset is empty: {path}")
    if len(set(columns)) != len(columns):
        raise MasterTableError(f"Sample-status column asset contains duplicate columns: {path}")
    return columns


def load_metadata(path: Path) -> tuple[list[str], str, dict[str, dict[str, str]]]:
    """Load metadata and preserve the original header order."""
    header, rows = read_table(path)
    key_column = detect_key_column(header, ("Accession", "accession"))
    index = build_index(rows, key_column, "metadata table")
    return header, key_column, index


def build_supplemental_metadata_map(sample_row: dict[str, str]) -> dict[str, str]:
    """Map validated sample extras onto normalised metadata-style names."""
    supplemental: dict[str, str] = {}
    for column, value in sample_row.items():
        if column in VALIDATED_SAMPLE_SUPPLEMENTAL_EXCLUSIONS:
            continue
        if not value or value == "NA":
            continue
        supplemental[normalise_key(column)] = value
    return supplemental


def build_metadata_row(
    sample_row: dict[str, str],
    metadata_header: Sequence[str],
    metadata_key_column: str,
    metadata_index: dict[str, dict[str, str]],
    assembly_stats_row: dict[str, str] | None = None,
) -> dict[str, str]:
    """Build the preserved metadata block for one sample."""
    accession = sample_row["accession"]
    existing_row = metadata_index.get(accession)
    if existing_row is not None:
        metadata_row = {
            column: existing_row.get(column, "NA") or "NA" for column in metadata_header
        }
    else:
        metadata_row = {column: "NA" for column in metadata_header}
        metadata_row[metadata_key_column] = accession
        if sample_row.get("is_new") == "true":
            supplemental = build_supplemental_metadata_map(sample_row)
            for column in metadata_header:
                normalised = normalise_key(column)
                if normalised in supplemental:
                    metadata_row[column] = supplemental[normalised]

    if assembly_stats_row:
        for metadata_column, stats_column in ASSEMBLY_METADATA_COLUMN_MAP.items():
            header_column = find_column_by_normalised_name(metadata_header, metadata_column)
            if header_column is None:
                continue
            if not is_missing(metadata_row.get(header_column)):
                continue
            stats_value = assembly_stats_row.get(stats_column, "NA")
            if is_missing(stats_value):
                continue
            metadata_row[header_column] = stats_value
    return metadata_row


def resolve_assembly_metric_value(
    metadata_row: dict[str, str],
    assembly_stats_row: dict[str, str] | None,
    metadata_column_name: str,
) -> str:
    """Resolve one assembly metric from metadata first, then in-house stats."""
    metadata_value = detect_metadata_value(metadata_row, metadata_column_name)
    if not is_missing(metadata_value):
        return metadata_value

    stats_column = ASSEMBLY_METADATA_COLUMN_MAP[metadata_column_name]
    if assembly_stats_row is None:
        return "NA"

    stats_value = assembly_stats_row.get(stats_column, "NA")
    return "NA" if is_missing(stats_value) else stats_value


def load_taxonomy_index(path: Path | None) -> dict[str, dict[str, str]]:
    """Load taxonomy rows keyed by Tax_ID."""
    if path is None:
        return {}
    header, rows = read_table(path, delimiter="\t")
    key_column = detect_key_column(header, ("Tax_ID", "tax_id"))
    require_columns(header, TAXONOMY_COLUMNS, "taxonomy table")
    return build_index(rows, key_column, "taxonomy table")


def load_optional_accession_index(
    path: Path | None,
    required_columns: Sequence[str],
    table_name: str,
    validated_accessions: set[str],
) -> dict[str, dict[str, str]]:
    """Load an optional accession-keyed TSV with required columns."""
    if path is None:
        return {}
    header, rows = read_table(path, delimiter="\t")
    key_column = detect_key_column(header, ("accession", "Accession"))
    require_columns(header, required_columns, table_name)
    index = build_index(rows, key_column, table_name)
    ensure_index_keys_are_valid(index, validated_accessions, table_name)
    return index


def load_assembly_stats_index(
    path: Path | None,
    validated_accessions: set[str],
) -> dict[str, dict[str, str]]:
    """Load optional in-house assembly stats keyed by accession."""
    return load_optional_accession_index(
        path,
        ASSEMBLY_STATS_COLUMNS,
        "assembly stats table",
        validated_accessions,
    )


def load_busco_index(
    paths: Sequence[Path],
    validated_accessions: set[str],
) -> tuple[dict[str, dict[str, BuscoLineageSummary]], set[str]]:
    """Load per-lineage BUSCO summaries into one accession-keyed mapping."""
    busco_index: dict[str, dict[str, BuscoLineageSummary]] = {}
    seen_busco_columns: set[str] = set()

    for path in paths:
        header, rows = read_table(path, delimiter="\t")
        key_column = detect_key_column(header, ("accession", "Accession"))
        busco_columns = [column for column in header if column.startswith("BUSCO_")]
        if len(busco_columns) != 1:
            raise MasterTableError(
                f"BUSCO table must contain exactly one BUSCO_<lineage> column: {path}"
            )
        busco_column = busco_columns[0]
        if busco_column in seen_busco_columns:
            raise MasterTableError(f"Duplicate BUSCO lineage column supplied: {busco_column}")
        seen_busco_columns.add(busco_column)

        table_index = build_index(rows, key_column, f"BUSCO table {path.name}")
        ensure_index_keys_are_valid(
            table_index,
            validated_accessions,
            f"BUSCO table {path.name}",
        )
        for accession, row in table_index.items():
            busco_index.setdefault(accession, {})
            busco_value = row.get(busco_column, "") or "NA"
            busco_status = row.get("busco_status", "")
            busco_warnings = row.get("warnings", "")
            if not busco_status:
                busco_status = "failed" if busco_value == "NA" else "done"
            if not busco_warnings and busco_status == "failed":
                busco_warnings = "busco_summary_failed"
            busco_index[accession][busco_column] = BuscoLineageSummary(
                value=busco_value,
                status=busco_status,
                warnings=busco_warnings,
            )

    return busco_index, seen_busco_columns


def resolve_append_columns(
    path: Path | None = None,
    busco_lineages: Sequence[str] | None = None,
) -> tuple[list[str], tuple[str, ...]]:
    """Resolve one append-column contract from an override, lineages, or defaults."""
    try:
        return master_table_contract.resolve_append_columns(
            path=path,
            busco_lineages=busco_lineages,
        )
    except ValueError as error:
        raise MasterTableError(str(error)) from error


def derive_taxonomy_values(
    metadata_row: dict[str, str],
    taxonomy_index: dict[str, dict[str, str]],
) -> dict[str, str]:
    """Resolve taxonomy columns by Tax_ID or return NA-filled values."""
    tax_id_column = next(
        (column for column in metadata_row if normalise_key(column) == "taxid"),
        None,
    )
    blank_taxonomy = {column: "NA" for column in TAXONOMY_COLUMNS}
    if tax_id_column is None:
        return blank_taxonomy
    tax_id = metadata_row.get(tax_id_column, "NA")
    if not tax_id or tax_id == "NA":
        return blank_taxonomy
    taxonomy_row = taxonomy_index.get(tax_id)
    if taxonomy_row is None:
        return blank_taxonomy
    return {column: taxonomy_row.get(column, "NA") or "NA" for column in TAXONOMY_COLUMNS}


def find_column_by_normalised_name(
    columns: Sequence[str],
    target: str,
) -> str | None:
    """Find a column by case- and punctuation-insensitive name."""
    normalised_target = normalise_key(target)
    for column in columns:
        if normalise_key(column) == normalised_target:
            return column
    return None


def detect_metadata_value(metadata_row: dict[str, str], column_name: str) -> str:
    """Return one metadata value by normalised column name or `NA`."""
    column = find_column_by_normalised_name(tuple(metadata_row), column_name)
    if column is None:
        return "NA"
    value = metadata_row.get(column, "").strip()
    return "NA" if is_missing(value) else value


def choose_assembly_level(sample_row: dict[str, str], metadata_row: dict[str, str]) -> str:
    """Choose the ANI assembly level, using the sample manifest for new genomes."""
    if sample_row.get("is_new") == "true":
        value = sample_row.get("assembly_level", "NA")
        return "NA" if is_missing(value) else value
    return detect_metadata_value(metadata_row, "Assembly_Level")


def build_master_row(
    sample_row: dict[str, str],
    metadata_header: Sequence[str],
    metadata_key_column: str,
    metadata_index: dict[str, dict[str, str]],
    assembly_stats_index: dict[str, dict[str, str]],
    taxonomy_index: dict[str, dict[str, str]],
    checkm2_index: dict[str, dict[str, str]],
    sixteen_s_index: dict[str, dict[str, str]],
    busco_index: dict[str, dict[str, BuscoLineageSummary]],
    codetta_index: dict[str, dict[str, str]],
    ccfinder_index: dict[str, dict[str, str]],
    ani_index: dict[str, dict[str, str]],
    append_columns: Sequence[str],
) -> dict[str, str]:
    """Build one final master-table row for a validated sample."""
    accession = sample_row["accession"]
    assembly_stats_row = assembly_stats_index.get(accession)
    metadata_row = build_metadata_row(
        sample_row,
        metadata_header,
        metadata_key_column,
        metadata_index,
        assembly_stats_row=assembly_stats_row,
    )
    row = {column: metadata_row[column] for column in metadata_header}

    derived_values = {column: "NA" for column in append_columns}
    derived_values["is_new"] = sample_row["is_new"]
    derived_values.update(derive_taxonomy_values(metadata_row, taxonomy_index))

    for column in CHECKM2_COLUMNS:
        if accession in checkm2_index and column in checkm2_index[accession]:
            derived_values[column] = checkm2_index[accession][column] or "NA"

    if assembly_stats_row is not None:
        for output_column, stats_column in ASSEMBLY_DERIVED_COLUMN_MAP.items():
            stats_value = assembly_stats_row.get(stats_column, "NA")
            derived_values[output_column] = "NA" if is_missing(stats_value) else stats_value

    if accession in sixteen_s_index and "16S" in sixteen_s_index[accession]:
        derived_values["16S"] = sixteen_s_index[accession]["16S"] or "NA"

    if accession in busco_index:
        for column, summary in busco_index[accession].items():
            if column in derived_values:
                derived_values[column] = summary.value or "NA"

    for column in CODETTA_COLUMNS:
        if accession in codetta_index and column in codetta_index[accession]:
            derived_values[column] = codetta_index[accession][column] or "NA"

    for column in CRISPR_COLUMNS:
        if accession in ccfinder_index and column in ccfinder_index[accession]:
            derived_values[column] = ccfinder_index[accession][column] or "NA"

    for column in ANI_COLUMNS:
        if accession in ani_index and column in ani_index[accession]:
            derived_values[column] = ani_index[accession][column] or "NA"

    row.update(derived_values)
    return row


def write_tsv(path: Path, header: Sequence[str], rows: Sequence[dict[str, str]]) -> None:
    """Write a TSV with a fixed header."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(header),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_build(args: argparse.Namespace) -> None:
    """Build the final master table from validated and derived inputs."""
    validated_samples = load_validated_samples(args.validated_samples)
    validated_accessions = {row["accession"] for row in validated_samples}
    metadata_header, metadata_key_column, metadata_index = load_metadata(args.metadata)
    append_columns, busco_lineages = resolve_append_columns(
        args.append_columns,
        busco_lineages=args.busco_lineage,
    )

    try:
        expected_header = master_table_contract.build_master_table_columns(
            metadata_header,
            busco_lineages,
        )
    except ValueError as error:
        raise MasterTableError(str(error)) from error
    if append_columns != expected_header[len(metadata_header) :]:
        raise MasterTableError(
            "Append-column contract does not match the BUSCO-aware master-table contract."
        )

    taxonomy_index = load_taxonomy_index(args.taxonomy)
    assembly_stats_index = load_assembly_stats_index(args.assembly_stats, validated_accessions)
    checkm2_index = load_optional_accession_index(
        args.checkm2,
        CHECKM2_COLUMNS,
        "CheckM2 summary table",
        validated_accessions,
    )
    sixteen_s_index = load_optional_accession_index(
        args.sixteen_s_status,
        SIXTEEN_S_COLUMNS,
        "16S status table",
        validated_accessions,
    )
    codetta_index = load_optional_accession_index(
        args.codetta_summary,
        CODETTA_COLUMNS,
        "Codetta summary table",
        validated_accessions,
    )
    ccfinder_index = load_optional_accession_index(
        args.ccfinder_strains,
        CRISPR_COLUMNS,
        "CRISPR strain summary table",
        validated_accessions,
    )
    ani_index = load_optional_accession_index(
        args.ani,
        ANI_COLUMNS,
        "ANI summary table",
        validated_accessions,
    )
    busco_index, provided_busco_columns = load_busco_index(args.busco, validated_accessions)
    expected_busco_columns = {f"BUSCO_{lineage}" for lineage in busco_lineages}
    unexpected_busco_columns = sorted(provided_busco_columns - expected_busco_columns)
    if unexpected_busco_columns:
        raise MasterTableError(
            "BUSCO summaries contain lineage columns not present in the append-column "
            "contract: " + ", ".join(unexpected_busco_columns)
        )

    rows: list[dict[str, str]] = []
    for sample_row in validated_samples:
        master_row = build_master_row(
            sample_row=sample_row,
            metadata_header=metadata_header,
            metadata_key_column=metadata_key_column,
            metadata_index=metadata_index,
            assembly_stats_index=assembly_stats_index,
            taxonomy_index=taxonomy_index,
            checkm2_index=checkm2_index,
            sixteen_s_index=sixteen_s_index,
            busco_index=busco_index,
            codetta_index=codetta_index,
            ccfinder_index=ccfinder_index,
            ani_index=ani_index,
            append_columns=append_columns,
        )
        rows.append(master_row)

    output_header = [*metadata_header, *append_columns]
    if len(rows) != len({row[metadata_key_column] for row in rows}):
        raise MasterTableError("Final master table contains duplicate accession rows.")

    write_tsv(args.output, output_header, rows)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the master-table builder CLI."""
    args = parse_args(argv)
    configure_logging()
    try:
        run_build(args)
    except MasterTableError as error:
        LOGGER.error(str(error))
        return 1

    LOGGER.info("Wrote master table to %s.", args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
