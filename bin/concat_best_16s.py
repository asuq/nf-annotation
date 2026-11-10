#!/usr/bin/env python3
"""Concatenate per-sample best 16S FASTA files into cohort outputs."""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Sequence

from atypical_warnings import classify_atypical_warnings

LOGGER = logging.getLogger(__name__)

REQUIRED_INPUT_COLUMNS = (
    "accession",
    "internal_id",
    "status_tsv",
    "best_16s_fasta",
    "qc_tsv",
)
REQUIRED_STATUS_COLUMNS = (
    "accession",
    "16S",
    "best_16S_header",
    "best_16S_length",
    "warnings",
)
COHORT_MANIFEST_COLUMNS = (
    "accession",
    "16S",
    "best_16S_header",
    "best_16S_length",
    "include_in_all_best_16S",
    "low_quality",
    "warnings",
)
REQUIRED_QC_COLUMNS = ("accession", "Low_quality")
COHORT_KINDS = (
    "intact",
    "partial",
    "low_quality_intact",
    "low_quality_partial",
)
MISSING_VALUE_TOKENS = {"", "na", "n/a", "null", "none"}


class Cohort16SError(RuntimeError):
    """Raised when the cohort 16S inputs are malformed."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Concatenate cohort-eligible best 16S FASTA files."
    )
    parser.add_argument(
        "--inputs",
        required=True,
        type=Path,
        help=(
            "TSV listing accession, internal_id, status_tsv, best_16s_fasta, "
            "and qc_tsv paths."
        ),
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        help="Optional metadata table used to apply the atypical intact-cohort rule.",
    )
    parser.add_argument(
        "--output-fasta",
        required=True,
        type=Path,
        help="Path to the concatenated cohort 16S FASTA file.",
    )
    parser.add_argument(
        "--output-manifest",
        required=True,
        type=Path,
        help="Path to the cohort manifest TSV.",
    )
    parser.add_argument(
        "--cohort-kind",
        choices=COHORT_KINDS,
        default="intact",
        help="Which cohort subset to build.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def is_missing(value: str | None) -> bool:
    """Return True when one scalar should be treated as missing."""
    if value is None:
        return True
    return value.strip().lower() in MISSING_VALUE_TOKENS


def normalise_key(value: str) -> str:
    """Normalise one header-like token for case-insensitive matching."""
    return "".join(character for character in value.casefold() if character.isalnum())


def sniff_delimiter(path: Path) -> str:
    """Detect whether a metadata table is CSV or TSV."""
    sample = path.read_text(encoding="utf-8").splitlines()
    if not sample:
        raise Cohort16SError(f"Metadata table is empty: {path}")
    sample_text = "\n".join(sample[:5])
    try:
        return csv.Sniffer().sniff(sample_text, delimiters=",\t").delimiter
    except csv.Error:
        if "\t" in sample_text:
            return "\t"
        if "," in sample_text:
            return ","
        raise Cohort16SError(f"Could not detect delimiter for metadata table: {path}")


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV file into a header list and row dictionaries."""
    if not path.is_file():
        raise Cohort16SError(f"Missing cohort 16S input file: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        raise Cohort16SError(f"Empty TSV file: {path}")

    header = [value.strip() for value in rows[0]]
    data_rows: list[dict[str, str]] = []
    for line_number, values in enumerate(rows[1:], start=2):
        if len(values) != len(header):
            raise Cohort16SError(
                f"TSV file {path} has {len(values)} columns on line {line_number}, "
                f"expected {len(header)}."
            )
        data_rows.append(
            {
                column: value.strip()
                for column, value in zip(header, values, strict=True)
            }
        )
    return header, data_rows


def read_metadata_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a CSV or TSV metadata file into header and rows."""
    if not path.is_file():
        raise Cohort16SError(f"Missing metadata table: {path}")

    delimiter = sniff_delimiter(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        rows = [row for row in reader if any(cell.strip() for cell in row)]

    if not rows:
        raise Cohort16SError(f"Metadata table is empty: {path}")

    header = [column.strip() for column in rows[0]]
    if any(not column for column in header):
        raise Cohort16SError(f"Metadata table has an empty header name: {path}")
    if len(set(header)) != len(header):
        raise Cohort16SError(f"Metadata table contains duplicate headers: {path}")

    records: list[dict[str, str]] = []
    for line_number, values in enumerate(rows[1:], start=2):
        if len(values) != len(header):
            raise Cohort16SError(
                f"Metadata table {path} has {len(values)} columns on line {line_number}, "
                f"expected {len(header)}."
            )
        records.append(
            {
                column: value.strip()
                for column, value in zip(header, values, strict=True)
            }
        )
    return header, records


def require_columns(
    header: Sequence[str],
    required_columns: Sequence[str],
    path: Path,
) -> None:
    """Raise when a TSV header is missing required columns."""
    missing = [column for column in required_columns if column not in header]
    if missing:
        raise Cohort16SError(
            f"TSV file {path} is missing required columns: {', '.join(missing)}"
        )


def load_cohort_rows(path: Path) -> list[dict[str, str]]:
    """Load and validate the cohort-input manifest rows."""
    header, rows = read_tsv(path)
    require_columns(header, REQUIRED_INPUT_COLUMNS, path)
    return rows


def load_status_row(path: Path) -> dict[str, str]:
    """Load and validate a single per-sample 16S status row."""
    header, rows = read_tsv(path)
    require_columns(header, REQUIRED_STATUS_COLUMNS, path)
    if len(rows) != 1:
        raise Cohort16SError(f"Status TSV must contain exactly one data row: {path}")
    return rows[0]


def load_qc_row(path: Path) -> dict[str, str]:
    """Load and validate a single per-sample CheckM2 QC row."""
    header, rows = read_tsv(path)
    require_columns(header, REQUIRED_QC_COLUMNS, path)
    if len(rows) != 1:
        raise Cohort16SError(f"QC TSV must contain exactly one data row: {path}")
    return rows[0]


def find_column_by_normalised_name(
    header: Sequence[str],
    column_name: str,
) -> str | None:
    """Find one column in a header by case-insensitive normalised name."""
    expected = normalise_key(column_name)
    matches = [column for column in header if normalise_key(column) == expected]
    if len(matches) > 1:
        raise Cohort16SError(
            f"Metadata table contains duplicate {column_name!r} columns."
        )
    if not matches:
        return None
    return matches[0]


def detect_accession_column(header: Sequence[str]) -> str:
    """Detect the metadata accession column."""
    aliases = {"accession", "Accession"}
    matches = [
        column
        for column in header
        if normalise_key(column) in {normalise_key(alias) for alias in aliases}
    ]
    if len(matches) > 1:
        raise Cohort16SError(
            f"Metadata table contains multiple accession columns: {matches}"
        )
    if not matches:
        raise Cohort16SError(
            "Metadata table is missing an accession column. Expected one of: "
            "Accession, accession"
        )
    return matches[0]


def load_metadata_atypical_index(metadata_path: Path | None) -> dict[str, str]:
    """Load accession to raw atypical warning mappings from metadata."""
    if metadata_path is None:
        return {}

    header, rows = read_metadata_table(metadata_path)
    accession_column = detect_accession_column(header)
    atypical_column = find_column_by_normalised_name(header, "Atypical_Warnings")
    if atypical_column is None:
        return {}

    index: dict[str, str] = {}
    for row in rows:
        accession = row.get(accession_column, "").strip()
        if not accession:
            raise Cohort16SError(
                f"Metadata table {metadata_path} contains an empty accession."
            )
        if accession in index:
            raise Cohort16SError(
                f"Metadata table contains duplicate accession {accession!r}."
            )
        index[accession] = row.get(atypical_column, "")
    return index


def append_fasta_content(
    output_handle,
    fasta_path: Path,
    *,
    internal_id: str,
) -> None:
    """Append one non-empty FASTA file with its record header rewritten."""
    if not fasta_path.is_file():
        raise Cohort16SError(f"Missing best 16S FASTA file: {fasta_path}")
    content = fasta_path.read_text(encoding="utf-8")
    if not content.strip():
        raise Cohort16SError(
            f"Included cohort FASTA file is empty but marked for inclusion: {fasta_path}"
        )
    lines = content.splitlines()
    if not lines or not lines[0].startswith(">"):
        raise Cohort16SError(f"Best 16S FASTA file is malformed: {fasta_path}")
    output_handle.write(f">{internal_id}\n")
    if len(lines) > 1:
        output_handle.write("\n".join(lines[1:]))
        output_handle.write("\n")


def should_include_in_cohort(
    status_row: dict[str, str],
    qc_row: dict[str, str],
    cohort_kind: str,
    metadata_atypical_index: dict[str, str],
) -> bool:
    """Return True when a sample belongs in the requested cohort subset."""
    status_value = status_row["16S"]
    low_quality = qc_row["Low_quality"] == "true"
    if cohort_kind == "low_quality_intact":
        return low_quality and status_value == "Yes"
    if cohort_kind == "low_quality_partial":
        return low_quality and status_value == "partial"
    if low_quality:
        return False
    if cohort_kind == "partial":
        return status_value == "partial"
    if status_value != "Yes":
        return False

    accession = status_row["accession"]
    if metadata_atypical_index:
        is_atypical, is_exception = classify_atypical_warnings(
            metadata_atypical_index.get(accession)
        )
        return not is_atypical or is_exception

    legacy_value = status_row.get("include_in_all_best_16S")
    return legacy_value == "true"


def build_manifest_row(
    status_row: dict[str, str],
    qc_row: dict[str, str],
    *,
    cohort_kind: str,
) -> dict[str, str]:
    """Build one stable cohort-manifest row."""
    return {
        "accession": status_row["accession"],
        "16S": status_row["16S"],
        "best_16S_header": status_row["best_16S_header"],
        "best_16S_length": status_row["best_16S_length"],
        "include_in_all_best_16S": "true" if cohort_kind == "intact" else "false",
        "low_quality": qc_row["Low_quality"] or "NA",
        "warnings": status_row["warnings"],
    }


def write_manifest(path: Path, rows: Sequence[dict[str, str]]) -> None:
    """Write the cohort manifest TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(COHORT_MANIFEST_COLUMNS),
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row[column] for column in COHORT_MANIFEST_COLUMNS})


def build_cohort_fasta(
    inputs_path: Path,
    metadata_path: Path | None,
    output_fasta: Path,
    output_manifest: Path,
    cohort_kind: str,
) -> None:
    """Build the cohort FASTA and manifest from per-sample 16S outputs."""
    manifest_rows: list[dict[str, str]] = []
    cohort_rows = load_cohort_rows(inputs_path)
    metadata_atypical_index = load_metadata_atypical_index(metadata_path)

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with output_fasta.open("w", encoding="utf-8") as fasta_handle:
        for cohort_row in cohort_rows:
            status_path = Path(cohort_row["status_tsv"]).expanduser()
            best_fasta_path = Path(cohort_row["best_16s_fasta"]).expanduser()
            qc_path = Path(cohort_row["qc_tsv"]).expanduser()
            internal_id = cohort_row["internal_id"].strip()
            status_row = load_status_row(status_path)
            qc_row = load_qc_row(qc_path)
            if status_row["accession"] != cohort_row["accession"]:
                raise Cohort16SError(
                    "Cohort input accession does not match the status TSV accession "
                    f"for {status_path}."
                )
            if qc_row["accession"] != cohort_row["accession"]:
                raise Cohort16SError(
                    "Cohort input accession does not match the QC TSV accession "
                    f"for {qc_path}."
                )
            if not internal_id:
                raise Cohort16SError(
                    f"Cohort input row is missing internal_id for {status_path}."
                )
            if not should_include_in_cohort(
                status_row=status_row,
                qc_row=qc_row,
                cohort_kind=cohort_kind,
                metadata_atypical_index=metadata_atypical_index,
            ):
                continue
            append_fasta_content(
                fasta_handle,
                best_fasta_path,
                internal_id=internal_id,
            )
            manifest_rows.append(
                build_manifest_row(status_row, qc_row, cohort_kind=cohort_kind)
            )

    write_manifest(output_manifest, manifest_rows)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the cohort 16S concatenation CLI."""
    args = parse_args(argv)
    configure_logging()
    try:
        build_cohort_fasta(
            inputs_path=args.inputs,
            metadata_path=args.metadata,
            output_fasta=args.output_fasta,
            output_manifest=args.output_manifest,
            cohort_kind=args.cohort_kind,
        )
    except Cohort16SError as error:
        LOGGER.error(str(error))
        return 1

    LOGGER.info("Wrote cohort 16S FASTA to %s.", args.output_fasta)
    return 0


if __name__ == "__main__":
    sys.exit(main())
