#!/usr/bin/env python3
"""Merge paired CheckM2 reports and assign gcode and QC status."""

from __future__ import annotations

import argparse
import csv
import logging
import math
import sys
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation, localcontext
from pathlib import Path
from typing import Sequence

from master_table_contract import GCODE_PROVENANCE_COLUMNS

LOGGER = logging.getLogger(__name__)
GCODE_RULE = "mean_gene_length_ratio"
LENGTH_RATIO_THRESHOLD = Decimal("1.5")

OUTPUT_COLUMNS = (
    "accession",
    "Completeness_gcode4",
    "Completeness_gcode11",
    "Contamination_gcode4",
    "Contamination_gcode11",
    "Coding_Density_gcode4",
    "Coding_Density_gcode11",
    "Average_Gene_Length_gcode4",
    "Average_Gene_Length_gcode11",
    "Total_Coding_Sequences_gcode4",
    "Total_Coding_Sequences_gcode11",
    "Gcode",
    *GCODE_PROVENANCE_COLUMNS,
    "Low_quality",
    "checkm2_status",
    "warnings",
)
METRIC_ALIASES = {
    "Completeness": ("Completeness",),
    "Contamination": ("Contamination",),
    "Coding_Density": ("Coding_Density",),
    "Average_Gene_Length": ("Average_Gene_Length", "Avg_Gene_Length"),
    "Total_Coding_Sequences": ("Total_Coding_Sequences",),
}
SHARED_STAT_ALIASES = {
    "Name": ("Name",),
    "Genome_Size": ("Genome_Size",),
    "GC_Content": ("GC_Content",),
    "Contig_N50": ("Contig_N50", "N50"),
}


@dataclass(frozen=True)
class ParsedCheckM2Report:
    """Represent the parsed numeric metrics from one CheckM2 report."""

    metrics: dict[str, Decimal]
    shared_stats: dict[str, Decimal | str]
    translation_table: int


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Summarise paired CheckM2 reports for one sample."
    )
    parser.add_argument(
        "--accession",
        required=True,
        help="Original sample accession for the output row.",
    )
    parser.add_argument(
        "--gcode4-report",
        required=True,
        type=Path,
        help="Path to the CheckM2 quality report for translation table 4.",
    )
    parser.add_argument(
        "--gcode11-report",
        required=True,
        type=Path,
        help="Path to the CheckM2 quality report for translation table 11.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Path to the combined per-sample QC TSV.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def read_report_rows(path: Path) -> list[dict[str, str]]:
    """Read a TSV report into row dictionaries."""
    if not path.is_file():
        raise ValueError(f"Missing CheckM2 report: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames
        if not header or len(header) != len(set(header)):
            raise ValueError(f"Missing or duplicate CheckM2 report columns: {path}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"CheckM2 report is empty: {path}")
    if len(rows) != 1:
        raise ValueError(
            f"CheckM2 report must contain exactly one data row: {path}"
        )
    if None in rows[0] or None in rows[0].values():
        raise ValueError(f"Malformed CheckM2 report row: {path}")
    return rows


def first_present(row: dict[str, str], candidates: Sequence[str]) -> str:
    """Return the first present field from an alias list."""
    for candidate in candidates:
        if candidate in row:
            return row[candidate]
    raise ValueError(f"Missing required CheckM2 field aliases: {', '.join(candidates)}")


def parse_number(value: str, field_name: str, path: Path) -> Decimal:
    """Preserve a native decimal metric and reject non-finite/underflow values."""
    try:
        result = Decimal(value)
    except (InvalidOperation, TypeError, ValueError) as error:
        raise ValueError(
            f"Could not parse {field_name} value {value!r} in {path}."
        ) from error
    if not result.is_finite() or not math.isfinite(float(result)):
        raise ValueError(f"Non-finite {field_name} value {value!r} in {path}.")
    if result != 0 and float(result) == 0:
        raise ValueError(f"Underflow in native {field_name} value {value!r} in {path}.")
    return Decimal(0) if result.is_zero() else result


def parse_report(path: Path) -> ParsedCheckM2Report:
    """Parse the required metrics from a single CheckM2 quality report."""
    row = read_report_rows(path)[0]

    metrics = {
        metric_name: parse_number(
            first_present(row, aliases),
            metric_name,
            path,
        )
        for metric_name, aliases in METRIC_ALIASES.items()
    }
    if any(value < 0 for value in metrics.values()):
        raise ValueError(f"Negative CheckM2 metric in {path}.")
    if metrics["Completeness"] > 100 or metrics["Coding_Density"] > 1:
        raise ValueError(f"CheckM2 completeness or coding density is out of range: {path}")
    if metrics["Average_Gene_Length"] <= 0:
        raise ValueError(f"Average_Gene_Length must be strictly positive: {path}")
    coding_sequences = metrics["Total_Coding_Sequences"]
    if coding_sequences <= 0 or coding_sequences != coding_sequences.to_integral_value():
        raise ValueError(f"Total_Coding_Sequences must be a positive integer: {path}")

    shared_stats: dict[str, Decimal | str] = {}
    for stat_name, aliases in SHARED_STAT_ALIASES.items():
        value = None
        for alias in aliases:
            if alias in row:
                value = row[alias]
                break
        if value is None or not value.strip():
            continue
        if stat_name == "Name":
            shared_stats[stat_name] = value.strip()
        else:
            shared_stats[stat_name] = parse_number(value, stat_name, path)

    if "Name" not in shared_stats:
        raise ValueError(f"Missing CheckM2 genome Name: {path}")
    table = row.get("Translation_Table_Used")
    if table not in {"4", "11"}:
        raise ValueError(f"Missing or unsupported Translation_Table_Used in {path}.")
    return ParsedCheckM2Report(
        metrics=metrics, shared_stats=shared_stats, translation_table=int(table)
    )


def format_metric(value: Decimal | None) -> str:
    """Preserve the parsed decimal's precision, or emit `NA`."""
    if value is None:
        return "NA"
    if not value.is_finite():
        raise ValueError("Cannot publish a non-finite CheckM2 metric.")
    return str(value)


def reports_have_consistent_shared_stats(
    report_four: ParsedCheckM2Report,
    report_eleven: ParsedCheckM2Report,
) -> bool:
    """Check whether shared assembly statistics match between two reports."""
    shared_keys_four = set(report_four.shared_stats)
    shared_keys_eleven = set(report_eleven.shared_stats)
    if shared_keys_four != shared_keys_eleven:
        return False

    shared_keys = shared_keys_four
    for key in shared_keys:
        left = report_four.shared_stats[key]
        right = report_eleven.shared_stats[key]
        if isinstance(left, str) or isinstance(right, str):
            if left != right:
                return False
            continue
        if not math.isclose(left, right, rel_tol=0.0, abs_tol=1e-9):
            return False
    return True


def empty_output_row(accession: str) -> dict[str, str]:
    """Create a default NA-filled output row."""
    row = {column: "NA" for column in OUTPUT_COLUMNS}
    row["accession"] = accession
    row["Gcode_Rule"] = GCODE_RULE
    row["Gcode_Length_Ratio_Threshold"] = str(LENGTH_RATIO_THRESHOLD)
    row["warnings"] = ""
    return row


def assign_low_quality(
    row: dict[str, str],
    report: ParsedCheckM2Report,
) -> None:
    """Populate Low_quality from the chosen report metrics."""
    with localcontext() as context:
        context.prec = max(
            len(value.as_tuple().digits) + abs(value.as_tuple().exponent)
            for value in (report.metrics["Completeness"], report.metrics["Contamination"])
        ) + 10
        low_quality_score = (
            report.metrics["Completeness"] - 5 * report.metrics["Contamination"]
        )
    row["Low_quality"] = "true" if low_quality_score <= 50 else "false"


def assign_gcode_from_valid_pair(
    row: dict[str, str],
    report_four: ParsedCheckM2Report,
    report_eleven: ParsedCheckM2Report,
) -> None:
    """Assign gcode from a valid paired CheckM2 comparison."""
    length_four = report_four.metrics["Average_Gene_Length"]
    length_eleven = report_eleven.metrics["Average_Gene_Length"]
    # Compare decimal input values directly, so a binary-float division cannot
    # place an exact 1.5 boundary on the wrong side of the rule.
    with localcontext() as context:
        context.prec = max(
            len(length_four.as_tuple().digits), len(length_eleven.as_tuple().digits)
        ) + 20
        select_four = length_four > LENGTH_RATIO_THRESHOLD * length_eleven
        row["Gcode_Length_Ratio"] = str(length_four / length_eleven)
    row["Gcode"] = "4" if select_four else "11"
    row["Gcode_Selection_Reason"] = (
        "length_ratio_above_threshold" if select_four else "length_ratio_at_or_below_threshold"
    )
    assign_low_quality(row, report_four if select_four else report_eleven)


def build_output_row(
    accession: str,
    report_four: ParsedCheckM2Report | None,
    report_eleven: ParsedCheckM2Report | None,
    *,
    warnings: list[str],
) -> dict[str, str]:
    """Create the final output row from two optional parsed reports."""
    row = empty_output_row(accession)

    if report_four is not None:
        row["Completeness_gcode4"] = format_metric(report_four.metrics["Completeness"])
        row["Contamination_gcode4"] = format_metric(report_four.metrics["Contamination"])
        row["Coding_Density_gcode4"] = format_metric(report_four.metrics["Coding_Density"])
        row["Average_Gene_Length_gcode4"] = format_metric(
            report_four.metrics["Average_Gene_Length"]
        )
        row["Total_Coding_Sequences_gcode4"] = format_metric(
            report_four.metrics["Total_Coding_Sequences"]
        )

    if report_eleven is not None:
        row["Completeness_gcode11"] = format_metric(
            report_eleven.metrics["Completeness"]
        )
        row["Contamination_gcode11"] = format_metric(
            report_eleven.metrics["Contamination"]
        )
        row["Coding_Density_gcode11"] = format_metric(
            report_eleven.metrics["Coding_Density"]
        )
        row["Average_Gene_Length_gcode11"] = format_metric(
            report_eleven.metrics["Average_Gene_Length"]
        )
        row["Total_Coding_Sequences_gcode11"] = format_metric(
            report_eleven.metrics["Total_Coding_Sequences"]
        )

    valid_pair = (
        report_four is not None
        and report_eleven is not None
        and report_four.translation_table == 4
        and report_eleven.translation_table == 11
        and reports_have_consistent_shared_stats(report_four, report_eleven)
    )
    if valid_pair:
        assert report_four is not None and report_eleven is not None
        assign_gcode_from_valid_pair(
            row,
            report_four,
            report_eleven,
        )
    else:
        row["Gcode"] = "NA"
        row["Low_quality"] = "NA"
        row["Gcode_Selection_Reason"] = "invalid_report_pair"

    failure_warnings = {
        "checkm2_gcode4_failed",
        "checkm2_gcode11_failed",
        "inconsistent_shared_stats",
    }
    row["checkm2_status"] = (
        "done"
        if valid_pair and not any(warning in failure_warnings for warning in warnings)
        else "failed"
    )
    row["warnings"] = ";".join(dict.fromkeys(warnings))
    return row


def write_output(path: Path, row: dict[str, str]) -> None:
    """Write the final single-row TSV output."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(OUTPUT_COLUMNS), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(row)


def run_summary(
    accession: str,
    gcode4_report: Path,
    gcode11_report: Path,
    output: Path,
) -> None:
    """Summarise paired CheckM2 reports into one per-sample TSV row."""
    warnings: list[str] = []

    try:
        report_four = parse_report(gcode4_report)
        if report_four.translation_table != 4:
            raise ValueError(f"Expected translation table 4 in {gcode4_report}.")
    except ValueError as error:
        LOGGER.warning(str(error))
        warnings.append("checkm2_gcode4_failed")
        report_four = None

    try:
        report_eleven = parse_report(gcode11_report)
        if report_eleven.translation_table != 11:
            raise ValueError(f"Expected translation table 11 in {gcode11_report}.")
    except ValueError as error:
        LOGGER.warning(str(error))
        warnings.append("checkm2_gcode11_failed")
        report_eleven = None

    if (
        report_four is not None
        and report_eleven is not None
        and not reports_have_consistent_shared_stats(report_four, report_eleven)
    ):
        warnings.append("inconsistent_shared_stats")

    row = build_output_row(
        accession,
        report_four,
        report_eleven,
        warnings=warnings,
    )
    write_output(output, row)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CheckM2 summarisation CLI."""
    args = parse_args(argv)
    configure_logging()
    run_summary(
        accession=args.accession,
        gcode4_report=args.gcode4_report,
        gcode11_report=args.gcode11_report,
        output=args.output,
    )
    LOGGER.info("Wrote CheckM2 summary for %s.", args.accession)
    return 0


if __name__ == "__main__":
    sys.exit(main())
