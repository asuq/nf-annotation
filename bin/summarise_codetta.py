#!/usr/bin/env python3
"""Summarise Codetta outputs into one accession-keyed reporting row."""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from collections.abc import Sequence
from pathlib import Path

from Bio.Data import CodonTable

LOGGER = logging.getLogger(__name__)
PROFILE_NAME = "Pfam-A_enone.hmm"
OUTPUT_COLUMNS = (
    "accession",
    "Codetta_Genetic_Code",
    "Codetta_NCBI_Table_Candidates",
    "codetta_status",
    "warnings",
)
CODON_ORDER = tuple(
    f"{first}{second}{third}"
    for first in "TCAG"
    for second in "TCAG"
    for third in "TCAG"
)
CODON_LINE_PATTERN = re.compile(r"^(?P<codon>[TCAG]{3})\s+(?P<inference>\S+)\s+")
GENETIC_CODE_PATTERN = re.compile(r"^Genetic code:\s*(?P<code>\S+)\s*$")


class CodettaSummaryError(RuntimeError):
    """Raise when Codetta outputs cannot be summarised safely."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Summarise one Codetta output directory into one TSV row."
    )
    parser.add_argument(
        "--accession",
        required=True,
        help="Sample accession for the summary row.",
    )
    parser.add_argument(
        "--codetta-dir",
        required=True,
        type=Path,
        help="Path to the per-sample Codetta output directory.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output TSV path.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def join_tokens(tokens: Sequence[str]) -> str:
    """Join non-empty warning tokens without duplicates."""
    return ";".join(dict.fromkeys(token for token in tokens if token))


def write_tsv(path: Path, row: dict[str, str]) -> None:
    """Write one-row TSV output."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(OUTPUT_COLUMNS), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def parse_table_genetic_code(inference_path: Path) -> str:
    """Parse the Codetta codon-inference table into one 64-character code string."""
    codon_calls: dict[str, str] = {}
    in_table = False
    for raw_line in inference_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.rstrip()
        if line.startswith("# codon") and "inference" in line:
            in_table = True
            continue
        if not in_table:
            continue
        if line.startswith("#"):
            break
        if not line.strip():
            continue
        match = CODON_LINE_PATTERN.match(line)
        if match is None:
            raise CodettaSummaryError(
                f"Codetta inference table contains an invalid codon row: {line!r}"
            )
        codon = match.group("codon")
        if codon in codon_calls:
            raise CodettaSummaryError(
                f"Codetta inference table contains duplicate rows for {codon!r}."
            )
        inference = match.group("inference").upper()
        if len(inference) != 1:
            raise CodettaSummaryError(
                f"Codetta inference for {codon!r} is not one character: {inference!r}"
            )
        codon_calls[codon] = inference

    missing = [codon for codon in CODON_ORDER if codon not in codon_calls]
    if missing:
        raise CodettaSummaryError(
            "Codetta inference table is missing codons: " + ", ".join(missing[:5])
        )
    return "".join(codon_calls[codon] for codon in CODON_ORDER)


def parse_log_genetic_code(log_path: Path) -> str:
    """Parse the final genetic-code string printed to the Codetta log."""
    for raw_line in log_path.read_text(encoding="utf-8").splitlines():
        match = GENETIC_CODE_PATTERN.match(raw_line.strip())
        if match is None:
            continue
        code = match.group("code").strip().upper()
        if len(code) != 64:
            raise CodettaSummaryError(
                f"Codetta log genetic code has length {len(code)}, expected 64."
            )
        return code
    raise CodettaSummaryError(
        f"Codetta log is missing a 'Genetic code:' line: {log_path}"
    )


def resolve_genetic_code(codetta_dir: Path) -> str:
    """Resolve the final genetic code from the inference table, then the log."""
    inference_path = codetta_dir / "codetta_inference.txt"
    log_path = codetta_dir / "codetta.log"
    if inference_path.is_file() and inference_path.stat().st_size > 0:
        try:
            return parse_table_genetic_code(inference_path)
        except CodettaSummaryError:
            if not log_path.is_file() or log_path.stat().st_size == 0:
                raise
    if log_path.is_file() and log_path.stat().st_size > 0:
        return parse_log_genetic_code(log_path)
    raise CodettaSummaryError(
        f"Codetta output directory is missing usable inference files: {codetta_dir}"
    )


def build_ncbi_table_string(table: CodonTable.CodonTable) -> str:
    """Build one NCBI genetic-code string in Codetta's 64-codon order."""
    stop_codons = set(table.stop_codons)
    return "".join(
        table.forward_table.get(codon, "*" if codon in stop_codons else "?")
        for codon in CODON_ORDER
    )


def match_ncbi_tables(genetic_code: str) -> list[int]:
    """Return all compatible NCBI translation-table ids for one inferred code."""
    candidates: list[int] = []
    for table_id in sorted(CodonTable.generic_by_id):
        candidate_code = build_ncbi_table_string(CodonTable.generic_by_id[table_id])
        compatible = True
        for observed, expected in zip(genetic_code, candidate_code, strict=True):
            if observed == "?":
                continue
            if observed != expected:
                compatible = False
                break
        if compatible:
            candidates.append(table_id)
    return candidates


def build_success_row(accession: str, genetic_code: str) -> dict[str, str]:
    """Build the success output row from one inferred genetic code."""
    warnings: list[str] = []
    candidate_ids = match_ncbi_tables(genetic_code)
    if candidate_ids:
        candidate_value = ";".join(str(table_id) for table_id in candidate_ids)
        if len(candidate_ids) > 1:
            warnings.append("codetta_multiple_ncbi_tables")
    else:
        candidate_value = "unassigned"
        warnings.append("codetta_unassigned")
    return {
        "accession": accession,
        "Codetta_Genetic_Code": genetic_code,
        "Codetta_NCBI_Table_Candidates": candidate_value,
        "codetta_status": "done",
        "warnings": join_tokens(warnings),
    }


def build_failure_row(accession: str) -> dict[str, str]:
    """Build the failure output row."""
    return {
        "accession": accession,
        "Codetta_Genetic_Code": "NA",
        "Codetta_NCBI_Table_Candidates": "NA",
        "codetta_status": "failed",
        "warnings": "codetta_failed",
    }


def run_summary(args: argparse.Namespace) -> dict[str, str]:
    """Summarise one Codetta output directory into one row."""
    codetta_dir = args.codetta_dir.resolve()
    try:
        genetic_code = resolve_genetic_code(codetta_dir)
    except CodettaSummaryError as error:
        LOGGER.warning("%s", error)
        return build_failure_row(args.accession)
    return build_success_row(args.accession, genetic_code)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the Codetta summariser CLI."""
    args = parse_args(argv)
    configure_logging()
    row = run_summary(args)
    write_tsv(args.output, row)
    return 0


if __name__ == "__main__":
    sys.exit(main())
