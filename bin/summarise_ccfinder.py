#!/usr/bin/env python3
"""Summarise CRISPRCasFinder outputs into strain, contig, and CRISPR tables."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence


LOGGER = logging.getLogger(__name__)

STRAIN_COLUMNS = (
    "accession",
    "CRISPRS",
    "SPACERS_SUM",
    "CRISPR_FRAC",
    "ccfinder_status",
    "warnings",
)
CONTIG_COLUMNS = (
    "accession",
    "contig_id",
    "contig_length",
    "CRISPRS",
    "SPACERS_SUM",
    "CRISPR_FRAC",
)
CRISPR_COLUMNS = (
    "accession",
    "contig_id",
    "crispr_id",
    "evidence_level",
    "spacer_count",
    "start",
    "end",
    "crispr_length",
)


@dataclass(frozen=True)
class CrisprRecord:
    """Represent one retained CRISPR array."""

    accession: str
    contig_id: str
    crispr_id: str
    evidence_level: int
    spacer_count: int
    start: int
    end: int

    @property
    def crispr_length(self) -> int:
        """Return the inclusive genomic span of the CRISPR array."""
        return self.end - self.start + 1


@dataclass(frozen=True)
class ContigSummary:
    """Represent contig-level CRISPR summary statistics."""

    accession: str
    contig_id: str
    contig_length: int | None
    crisprs: int
    spacers_sum: int
    crispr_frac: float | None


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Summarise CRISPRCasFinder result.json for one sample."
    )
    parser.add_argument(
        "--accession",
        required=True,
        help="Original sample accession for the output tables.",
    )
    parser.add_argument(
        "--result-json",
        required=True,
        type=Path,
        help="Path to the CRISPRCasFinder result.json file.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Directory for strain, contig, and CRISPR summary tables.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def format_number(value: float | int | None) -> str:
    """Format optional numeric values for TSV output."""
    if value is None:
        return "NA"
    if isinstance(value, int):
        return str(value)
    if math.isclose(value, round(value), rel_tol=0.0, abs_tol=1e-9):
        return str(int(round(value)))
    return f"{value:.6f}".rstrip("0").rstrip(".")


def load_json(path: Path) -> Any:
    """Read a JSON document from disk."""
    if not path.is_file():
        raise ValueError(f"Missing CRISPRCasFinder result JSON: {path}")
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"Could not parse CRISPRCasFinder JSON: {path}") from error


def find_sequence_records(payload: Any) -> list[dict[str, Any]]:
    """Find the top-level CRISPRCasFinder sequence list in a JSON document."""
    if isinstance(payload, dict):
        for key in ("Sequences", "sequences", "Sequence", "sequence"):
            value = payload.get(key)
            if isinstance(value, list):
                return [item for item in value if isinstance(item, dict)]
        for value in payload.values():
            found = find_sequence_records(value)
            if found:
                return found
    if isinstance(payload, list) and all(isinstance(item, dict) for item in payload):
        return list(payload)
    return []


def find_first_value(record: dict[str, Any], keys: Sequence[str]) -> Any:
    """Return the first present dictionary value from a list of candidate keys."""
    for key in keys:
        if key in record:
            return record[key]
    return None


def parse_int(value: Any) -> int | None:
    """Convert a scalar to int when possible."""
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return None
        try:
            return int(stripped)
        except ValueError:
            return None
    return None


def extract_contig_id(sequence_record: dict[str, Any], index: int) -> str:
    """Extract a stable contig identifier from a sequence record."""
    value = find_first_value(
        sequence_record,
        ("Id", "ID", "SequenceID", "Sequence_Id", "Name", "name"),
    )
    if value is None:
        return f"contig_{index}"
    return str(value)


def extract_contig_length(sequence_record: dict[str, Any]) -> int | None:
    """Extract the contig length when CRISPRCasFinder reports it."""
    return parse_int(
        find_first_value(
            sequence_record,
            ("Length", "length", "SequenceLength", "Sequence_Length"),
        )
    )


def extract_crispr_entries(sequence_record: dict[str, Any]) -> list[dict[str, Any]]:
    """Extract CRISPR array entries from a sequence record."""
    value = find_first_value(
        sequence_record,
        ("Crisprs", "CRISPRs", "crisprs", "crisprs"),
    )
    if isinstance(value, list):
        return [entry for entry in value if isinstance(entry, dict)]
    return []


def extract_spacer_count(crispr_entry: dict[str, Any]) -> int | None:
    """Extract the number of spacers reported for a CRISPR array."""
    spacer_value = find_first_value(
        crispr_entry,
        (
            "Spacers",
            "spacers",
            "Spacer_Count",
            "Spacers_Count",
            "Number_Spacers",
        ),
    )
    if isinstance(spacer_value, list):
        return len(spacer_value)
    return parse_int(spacer_value)


def parse_crispr_records(
    accession: str, payload: Any
) -> tuple[list[CrisprRecord], list[ContigSummary], list[str]]:
    """Parse retained CRISPR arrays plus contig summaries from the JSON payload."""
    sequence_records = find_sequence_records(payload)
    if not sequence_records:
        raise ValueError("CRISPRCasFinder JSON does not contain any sequence records.")

    crispr_records: list[CrisprRecord] = []
    contig_summaries: list[ContigSummary] = []
    warnings: list[str] = []

    for index, sequence_record in enumerate(sequence_records, start=1):
        contig_id = extract_contig_id(sequence_record, index)
        contig_length = extract_contig_length(sequence_record)
        contig_crisprs: list[CrisprRecord] = []
        for crispr_index, crispr_entry in enumerate(
            extract_crispr_entries(sequence_record),
            start=1,
        ):
            evidence_level = parse_int(
                find_first_value(
                    crispr_entry,
                    ("Evidence_Level", "EvidenceLevel", "evidence_level"),
                )
            )
            if evidence_level is None:
                warnings.append("invalid_crispr_entry")
                continue
            if evidence_level == 1:
                continue

            start = parse_int(find_first_value(crispr_entry, ("Start", "start")))
            end = parse_int(find_first_value(crispr_entry, ("End", "end")))
            spacer_count = extract_spacer_count(crispr_entry)
            if start is None or end is None or spacer_count is None or end < start:
                warnings.append("invalid_crispr_entry")
                continue

            crispr_id = find_first_value(
                crispr_entry,
                ("Name", "name", "Id", "ID", "Crispr_Id", "CRISPR_Id"),
            )
            if crispr_id is None:
                crispr_id = f"{contig_id}_crispr_{crispr_index}"

            record = CrisprRecord(
                accession=accession,
                contig_id=contig_id,
                crispr_id=str(crispr_id),
                evidence_level=evidence_level,
                spacer_count=spacer_count,
                start=start,
                end=end,
            )
            crispr_records.append(record)
            contig_crisprs.append(record)

        total_crispr_length = sum(record.crispr_length for record in contig_crisprs)
        contig_summaries.append(
            ContigSummary(
                accession=accession,
                contig_id=contig_id,
                contig_length=contig_length,
                crisprs=len(contig_crisprs),
                spacers_sum=sum(record.spacer_count for record in contig_crisprs),
                crispr_frac=(
                    total_crispr_length / contig_length
                    if contig_length not in (None, 0)
                    else None
                ),
            )
        )

    return crispr_records, contig_summaries, warnings


def build_strain_row(
    accession: str,
    crispr_records: Sequence[CrisprRecord],
    contig_summaries: Sequence[ContigSummary],
    status: str,
    warnings: Sequence[str],
) -> dict[str, str]:
    """Build the strain-level summary row."""
    total_contig_length = sum(
        summary.contig_length
        for summary in contig_summaries
        if summary.contig_length is not None
    )
    total_crispr_length = sum(record.crispr_length for record in crispr_records)

    crispr_frac = None
    if status == "done":
        crispr_frac = (
            total_crispr_length / total_contig_length
            if total_contig_length > 0
            else 0.0
        )

    return {
        "accession": accession,
        "CRISPRS": format_number(len(crispr_records) if status == "done" else None),
        "SPACERS_SUM": format_number(
            sum(record.spacer_count for record in crispr_records)
            if status == "done"
            else None
        ),
        "CRISPR_FRAC": format_number(crispr_frac),
        "ccfinder_status": status,
        "warnings": ";".join(dict.fromkeys(warnings)),
    }


def build_contig_rows(
    contig_summaries: Sequence[ContigSummary],
) -> list[dict[str, str]]:
    """Build contig-level summary rows."""
    return [
        {
            "accession": summary.accession,
            "contig_id": summary.contig_id,
            "contig_length": format_number(summary.contig_length),
            "CRISPRS": format_number(summary.crisprs),
            "SPACERS_SUM": format_number(summary.spacers_sum),
            "CRISPR_FRAC": format_number(
                summary.crispr_frac if summary.crisprs else 0.0
            ),
        }
        for summary in contig_summaries
    ]


def build_crispr_rows(crispr_records: Sequence[CrisprRecord]) -> list[dict[str, str]]:
    """Build one row per retained CRISPR array."""
    return [
        {
            "accession": record.accession,
            "contig_id": record.contig_id,
            "crispr_id": record.crispr_id,
            "evidence_level": format_number(record.evidence_level),
            "spacer_count": format_number(record.spacer_count),
            "start": format_number(record.start),
            "end": format_number(record.end),
            "crispr_length": format_number(record.crispr_length),
        }
        for record in crispr_records
    ]


def write_tsv(
    path: Path, header: Sequence[str], rows: Sequence[dict[str, str]]
) -> None:
    """Write a TSV file with a fixed header."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(header), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CRISPRCasFinder summarisation CLI."""
    args = parse_args(argv)
    configure_logging()
    args.outdir.mkdir(parents=True, exist_ok=True)

    try:
        payload = load_json(args.result_json)
        crispr_records, contig_summaries, warnings = parse_crispr_records(
            args.accession,
            payload,
        )
    except ValueError as error:
        LOGGER.warning(str(error))
        strain_row = build_strain_row(
            args.accession,
            [],
            [],
            status="failed",
            warnings=["ccfinder_summary_failed"],
        )
        contig_rows = []
        crispr_rows = []
    else:
        strain_row = build_strain_row(
            args.accession,
            crispr_records,
            contig_summaries,
            status="done",
            warnings=warnings,
        )
        contig_rows = build_contig_rows(contig_summaries)
        crispr_rows = build_crispr_rows(crispr_records)

    write_tsv(args.outdir / "ccfinder_strains.tsv", STRAIN_COLUMNS, [strain_row])
    write_tsv(args.outdir / "ccfinder_contigs.tsv", CONTIG_COLUMNS, contig_rows)
    write_tsv(args.outdir / "ccfinder_crisprs.tsv", CRISPR_COLUMNS, crispr_rows)
    LOGGER.info("Wrote CRISPRCasFinder summaries for %s.", args.accession)
    return 0


if __name__ == "__main__":
    sys.exit(main())
