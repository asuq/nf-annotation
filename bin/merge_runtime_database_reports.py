#!/usr/bin/env python3
"""Merge runtime database preparation reports and render Nextflow arguments."""

from __future__ import annotations

import argparse
import csv
import shlex
import sys
from pathlib import Path
from typing import Sequence


REPORT_COLUMNS = ("component", "status", "source", "destination", "details")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Merge per-component runtime database reports and write a paste-ready "
            "Nextflow argument block."
        )
    )
    parser.add_argument(
        "--report",
        action="append",
        type=Path,
        required=True,
        help="One per-component TSV report. May be supplied multiple times.",
    )
    parser.add_argument(
        "--output-report",
        type=Path,
        required=True,
        help="Merged runtime database report path.",
    )
    parser.add_argument(
        "--output-args",
        type=Path,
        required=True,
        help="Paste-ready Nextflow argument file path.",
    )
    return parser.parse_args(argv)


def read_rows(path: Path) -> list[dict[str, str]]:
    """Read one TSV report file."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader)


def component_sort_key(row: dict[str, str]) -> tuple[int, str]:
    """Return a stable sort key for report rows."""
    ordering = {
        "taxdump": 0,
        "checkm2": 1,
        "codetta": 2,
        "busco_root": 3,
        "eggnog": 4,
    }
    return (ordering.get(row["component"], 99), row["component"])


def merge_rows(report_paths: Sequence[Path]) -> list[dict[str, str]]:
    """Merge and sort rows from all report files."""
    rows: list[dict[str, str]] = []
    for report_path in report_paths:
        rows.extend(read_rows(report_path))
    return sorted(rows, key=component_sort_key)


def write_merged_report(path: Path, rows: Sequence[dict[str, str]]) -> None:
    """Write the merged runtime database report."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_nextflow_args(rows: Sequence[dict[str, str]]) -> str:
    """Build the paste-ready Nextflow argument block."""
    flag_map = {
        "taxdump": "--taxdump",
        "checkm2": "--checkm2_db",
        "codetta": "--codetta_db",
        "busco_root": "--busco_db",
        "eggnog": "--eggnog_db",
        "cogclassifier": "--cogclassifier_db",
        "pfam": "--pfam_db",
        "kofam": "--kofam_db",
    }
    arguments = [
        (flag_map[row["component"]], row["destination"])
        for row in rows
        if row["component"] in flag_map
    ]
    lines = ["nextflow run . \\"]
    for index, (flag, value) in enumerate(arguments):
        suffix = " \\" if index < len(arguments) - 1 else ""
        lines.append(f"  {flag} {shlex.quote(value)}{suffix}")
    return "\n".join(lines) + "\n"


def write_nextflow_args(path: Path, rows: Sequence[dict[str, str]]) -> None:
    """Write the Nextflow argument block."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(build_nextflow_args(rows), encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    """Merge runtime database reports."""
    args = parse_args(argv)
    rows = merge_rows(args.report)
    write_merged_report(args.output_report, rows)
    write_nextflow_args(args.output_args, rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
