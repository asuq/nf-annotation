#!/usr/bin/env python3
"""Summarise Barrnap 16S calls for one sample."""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence, TypeVar


LOGGER = logging.getLogger(__name__)

STATUS_COLUMNS = (
    "accession",
    "16S",
    "best_16S_header",
    "best_16S_length",
    "warnings",
)
SIXTEEN_S_NAME = "16S_rRNA"
FASTA_HEADER_RE = re.compile(
    r"^(?P<rrna_name>[^:]+)::(?P<contig>.+):(?P<start0>\d+)-(?P<end0>\d+)\((?P<strand>[+-])\)$"
)
MatchKey = tuple[str, str, int, int, str]
RecordT = TypeVar("RecordT")


@dataclass(frozen=True)
class FastaRecord:
    """Represent one Barrnap FASTA record."""

    header: str
    sequence: str
    rrna_name: str
    contig: str
    start0: int
    end0: int
    strand: str


@dataclass(frozen=True)
class GffHit:
    """Represent one Barrnap GFF hit."""

    order_index: int
    score: float
    rrna_name: str
    contig: str
    start0: int
    end0: int
    strand: str
    is_partial: bool


@dataclass(frozen=True)
class RrnaHit:
    """Represent one matched Barrnap 16S hit."""

    order_index: int
    score: float
    is_partial: bool
    record: FastaRecord


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Summarise Barrnap GFF/FASTA outputs for one sample."
    )
    parser.add_argument(
        "--accession",
        required=True,
        help="Original sample accession for the output row.",
    )
    parser.add_argument(
        "--rrna-gff",
        required=True,
        type=Path,
        help="Path to the Barrnap GFF file.",
    )
    parser.add_argument(
        "--rrna-fasta",
        required=True,
        type=Path,
        help="Path to the Barrnap FASTA file.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Directory for 16S summary outputs.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def build_match_key(
    rrna_name: str,
    contig: str,
    start0: int,
    end0: int,
    strand: str,
) -> MatchKey:
    """Build a canonical Barrnap record key."""
    return rrna_name, contig, start0, end0, strand


def parse_fasta_header(header: str) -> tuple[str, str, int, int, str]:
    """Parse one Barrnap FASTA header into its canonical fields."""
    match = FASTA_HEADER_RE.fullmatch(header)
    if match is None:
        raise ValueError(f"Malformed Barrnap FASTA header: {header!r}")
    return (
        match.group("rrna_name"),
        match.group("contig"),
        int(match.group("start0")),
        int(match.group("end0")),
        match.group("strand"),
    )


def build_fasta_record(header: str, sequence_lines: Sequence[str]) -> FastaRecord:
    """Build one parsed Barrnap FASTA record from a header and sequence lines."""
    rrna_name, contig, start0, end0, strand = parse_fasta_header(header)
    return FastaRecord(
        header=header,
        sequence="".join(sequence_lines),
        rrna_name=rrna_name,
        contig=contig,
        start0=start0,
        end0=end0,
        strand=strand,
    )


def parse_fasta_records(path: Path) -> list[FastaRecord]:
    """Parse a Barrnap FASTA file into ordered records."""
    if not path.is_file():
        raise ValueError(f"Missing Barrnap FASTA file: {path}")

    records: list[FastaRecord] = []
    header: str | None = None
    sequence_lines: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(build_fasta_record(header, sequence_lines))
                header = line[1:].strip()
                sequence_lines = []
                continue
            if header is None:
                raise ValueError(f"Malformed FASTA file without a header: {path}")
            sequence_lines.append(line.strip())

    if header is not None:
        records.append(build_fasta_record(header, sequence_lines))
    return records


def parse_gff_attributes(attributes_field: str) -> dict[str, str]:
    """Parse one GFF attributes field into a key-value mapping."""
    attributes: dict[str, str] = {}
    for raw_part in attributes_field.split(";"):
        part = raw_part.strip()
        if not part:
            continue
        if "=" not in part:
            raise ValueError(f"Malformed GFF attribute: {part!r}")
        key, value = part.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not key:
            raise ValueError(f"Malformed GFF attribute key: {part!r}")
        attributes[key] = value
    return attributes


def parse_gff_hits(path: Path) -> list[GffHit]:
    """Parse Barrnap GFF rows into ordered hit metadata."""
    if not path.is_file():
        raise ValueError(f"Missing Barrnap GFF file: {path}")

    parsed_hits: list[GffHit] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed GFF line in {path}: {line}")

            attributes = parse_gff_attributes(fields[8])
            rrna_name = attributes.get("Name", "").strip()
            if not rrna_name:
                raise ValueError(f"Missing Name attribute in {path}: {line}")

            score_value = fields[5]
            if score_value == ".":
                score = float("inf")
            else:
                try:
                    score = float(score_value)
                except ValueError as error:
                    raise ValueError(
                        f"Could not parse Barrnap score {score_value!r} in {path}."
                    ) from error

            try:
                start1 = int(fields[3])
                end1 = int(fields[4])
            except ValueError as error:
                raise ValueError(
                    f"Could not parse GFF coordinates in {path}: {line}"
                ) from error

            parsed_hits.append(
                GffHit(
                    order_index=len(parsed_hits) + 1,
                    score=score,
                    rrna_name=rrna_name,
                    contig=fields[0].strip(),
                    start0=start1 - 1,
                    end0=end1,
                    strand=fields[6].strip(),
                    is_partial="partial" in fields[8].lower(),
                )
            )
    return parsed_hits


def build_keyed_records(
    records: Sequence[RecordT],
    path: Path,
    source_name: str,
    key_builder: Callable[[RecordT], MatchKey],
) -> dict[MatchKey, RecordT]:
    """Build a unique keyed mapping for parsed Barrnap records."""
    keyed_records: dict[MatchKey, RecordT] = {}
    for record in records:
        key = key_builder(record)
        if key in keyed_records:
            raise ValueError(
                f"Duplicate {source_name} record in {path} for key {key!r}."
            )
        keyed_records[key] = record
    return keyed_records


def pair_16s_hits(
    gff_hits: Sequence[GffHit],
    fasta_records: Sequence[FastaRecord],
    rrna_gff: Path,
    rrna_fasta: Path,
) -> list[RrnaHit]:
    """Match Barrnap 16S FASTA records to 16S GFF hits by canonical key."""
    build_keyed_records(
        records=gff_hits,
        path=rrna_gff,
        source_name="GFF",
        key_builder=lambda hit: build_match_key(
            hit.rrna_name, hit.contig, hit.start0, hit.end0, hit.strand
        ),
    )
    fasta_by_key = build_keyed_records(
        records=fasta_records,
        path=rrna_fasta,
        source_name="FASTA",
        key_builder=lambda record: build_match_key(
            record.rrna_name,
            record.contig,
            record.start0,
            record.end0,
            record.strand,
        ),
    )

    sixteen_s_gff_hits = [hit for hit in gff_hits if hit.rrna_name == SIXTEEN_S_NAME]
    sixteen_s_fasta_records = [
        record for record in fasta_records if record.rrna_name == SIXTEEN_S_NAME
    ]

    if not sixteen_s_gff_hits and not sixteen_s_fasta_records:
        return []

    gff_keys = {
        build_match_key(hit.rrna_name, hit.contig, hit.start0, hit.end0, hit.strand)
        for hit in sixteen_s_gff_hits
    }
    fasta_keys = {
        build_match_key(
            record.rrna_name,
            record.contig,
            record.start0,
            record.end0,
            record.strand,
        )
        for record in sixteen_s_fasta_records
    }
    if gff_keys != fasta_keys:
        raise ValueError("Barrnap 16S GFF and FASTA outputs do not match.")

    matched_hits: list[RrnaHit] = []
    for hit in sixteen_s_gff_hits:
        key = build_match_key(
            hit.rrna_name, hit.contig, hit.start0, hit.end0, hit.strand
        )
        record = fasta_by_key.get(key)
        if record is None:
            raise ValueError("Barrnap 16S GFF and FASTA outputs do not match.")
        matched_hits.append(
            RrnaHit(
                order_index=hit.order_index,
                score=hit.score,
                is_partial=hit.is_partial,
                record=record,
            )
        )
    return matched_hits


def choose_best_hit(hits: Sequence[RrnaHit]) -> RrnaHit | None:
    """Choose the best 16S hit using the locked Barrnap tie-break rules."""
    if not hits:
        return None

    intact_hits = [hit for hit in hits if not hit.is_partial]
    ranking_pool = intact_hits if intact_hits else list(hits)
    return min(
        ranking_pool,
        key=lambda hit: (
            hit.score,
            -len(hit.record.sequence),
            hit.order_index,
        ),
    )


def write_best_fasta(path: Path, hit: RrnaHit | None) -> None:
    """Write the selected best-hit FASTA record or an empty file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if hit is None:
        path.write_text("", encoding="utf-8")
        return

    wrapped_sequence = [
        hit.record.sequence[index : index + 80]
        for index in range(0, len(hit.record.sequence), 80)
    ]
    path.write_text(
        f">{hit.record.header}\n" + "\n".join(wrapped_sequence) + "\n",
        encoding="utf-8",
    )


def write_status(path: Path, row: dict[str, str]) -> None:
    """Write the one-row 16S status TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(STATUS_COLUMNS),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(row)


def build_invalid_output_row(accession: str, warnings: Sequence[str]) -> dict[str, str]:
    """Build the soft-fail output row for invalid Barrnap outputs."""
    return {
        "accession": accession,
        "16S": "NA",
        "best_16S_header": "NA",
        "best_16S_length": "NA",
        "warnings": ";".join(warnings),
    }


def summarise_hits(
    accession: str,
    rrna_gff: Path,
    rrna_fasta: Path,
    outdir: Path,
) -> None:
    """Summarise one sample's Barrnap outputs into 16S status files."""
    warnings: list[str] = []
    outdir.mkdir(parents=True, exist_ok=True)

    try:
        matched_hits = pair_16s_hits(
            parse_gff_hits(rrna_gff),
            parse_fasta_records(rrna_fasta),
            rrna_gff,
            rrna_fasta,
        )
    except ValueError as error:
        LOGGER.warning(str(error))
        warnings.append("invalid_barrnap_output")
        write_best_fasta(outdir / "best_16S.fna", None)
        write_status(
            outdir / "16S_status.tsv",
            build_invalid_output_row(accession, warnings),
        )
        return

    best_hit = choose_best_hit(matched_hits)
    if not matched_hits:
        status = "No"
    elif any(not hit.is_partial for hit in matched_hits):
        status = "Yes"
    else:
        status = "partial"

    row = {
        "accession": accession,
        "16S": status,
        "best_16S_header": best_hit.record.header if best_hit is not None else "NA",
        "best_16S_length": (
            str(len(best_hit.record.sequence)) if best_hit is not None else "NA"
        ),
        "warnings": "",
    }

    write_best_fasta(outdir / "best_16S.fna", best_hit)
    write_status(outdir / "16S_status.tsv", row)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the 16S summarisation CLI."""
    args = parse_args(argv)
    configure_logging()
    summarise_hits(
        accession=args.accession,
        rrna_gff=args.rrna_gff,
        rrna_fasta=args.rrna_fasta,
        outdir=args.outdir,
    )
    LOGGER.info("Wrote 16S summary for %s.", args.accession)
    return 0


if __name__ == "__main__":
    sys.exit(main())
