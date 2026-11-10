#!/usr/bin/env python3
"""Parse machine-readable BUSCO summaries into master-table-ready strings."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any, cast

LOGGER = logging.getLogger(__name__)

BUSCO_STRING_RE = re.compile(
    r"C:\s*\d+(?:\.\d+)?%\s*"
    r"\[\s*S:\s*\d+(?:\.\d+)?%,\s*D:\s*\d+(?:\.\d+)?%\s*\],\s*"
    r"F:\s*\d+(?:\.\d+)?%,\s*"
    r"M:\s*\d+(?:\.\d+)?%,\s*"
    r"n:\s*\d+",
    re.IGNORECASE,
)
PERCENT_ALIASES = {
    "C": (
        "completepercentage",
        "completepercent",
    ),
    "S": (
        "singlecopypercentage",
        "singlecopybuscospercentage",
        "singlepercentage",
    ),
    "D": (
        "duplicatedpercentage",
        "duplicatedbuscospercentage",
        "multicopypercentage",
    ),
    "F": (
        "fragmentedpercentage",
        "fragmentedbuscospercentage",
    ),
    "M": (
        "missingpercentage",
        "missingbuscospercentage",
    ),
}
COUNT_ALIASES = {
    "C": (
        "c",
        "complete",
        "completebuscos",
    ),
    "S": (
        "s",
        "singlecopy",
        "singlecopybuscos",
    ),
    "D": (
        "d",
        "duplicated",
        "duplicatedbuscos",
        "multicopy",
    ),
    "F": (
        "f",
        "fragmented",
        "fragmentedbuscos",
    ),
    "M": (
        "m",
        "missing",
        "missingbuscos",
    ),
    "n": (
        "n",
        "nmarkers",
        "nmarkerstotal",
        "totalbuscogroupssearched",
        "numberofbuscogroupssearched",
        "datasetnumberofbuscogenes",
    ),
}
COMPACT_PERCENTAGE_KEYS = ("C", "S", "D", "F", "M")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Summarise one BUSCO machine-readable summary file."
    )
    parser.add_argument(
        "--accession",
        required=True,
        help="Original sample accession for the output row.",
    )
    parser.add_argument(
        "--summary",
        required=True,
        type=Path,
        help="Path to a BUSCO machine-readable summary file.",
    )
    parser.add_argument(
        "--lineage",
        required=True,
        help="BUSCO lineage label for the summary being parsed.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Path to the summarised BUSCO TSV.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def busco_column_name(lineage: str) -> str:
    """Return the stable output column name for one BUSCO lineage."""
    return f"BUSCO_{lineage}"


def primary_busco_column(configured_lineages: Sequence[str]) -> str:
    """Return the BUSCO column name used as the primary ANI-scoring lineage."""
    lineages = tuple(
        lineage.strip() for lineage in configured_lineages if lineage.strip()
    )
    if not lineages:
        raise ValueError("At least one BUSCO lineage is required.")
    return busco_column_name(lineages[0])


def normalise_key(key: str) -> str:
    """Normalise a JSON key for alias matching."""
    return "".join(character for character in key.lower() if character.isalnum())


def collect_scalar_values(data: Any) -> dict[str, list[Any]]:
    """Collect scalar JSON values keyed by normalised field name."""
    collected: dict[str, list[Any]] = {}

    def _walk(node: Any) -> None:
        if isinstance(node, dict):
            for key, value in node.items():
                normalised = normalise_key(str(key))
                if isinstance(value, (dict, list)):
                    _walk(value)
                    continue
                collected.setdefault(normalised, []).append(value)
        elif isinstance(node, list):
            for item in node:
                _walk(item)

    _walk(data)
    return collected


def find_existing_busco_string(data: Any) -> str | None:
    """Find an existing BUSCO summary string anywhere in the JSON document."""
    if isinstance(data, str):
        match = BUSCO_STRING_RE.search(data)
        if match is None:
            return None
        return match.group(0)
    if isinstance(data, dict):
        for value in data.values():
            found = find_existing_busco_string(value)
            if found is not None:
                return found
    if isinstance(data, list):
        for item in data:
            found = find_existing_busco_string(item)
            if found is not None:
                return found
    return None


def coerce_float(value: Any) -> float | None:
    """Convert a scalar JSON value to float when possible."""
    if isinstance(value, bool):
        return None
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        stripped = value.strip().rstrip("%")
        if not stripped:
            return None
        try:
            return float(stripped)
        except ValueError:
            return None
    return None


def find_scalar_value(
    scalars: dict[str, list[Any]],
    aliases: Sequence[str],
) -> float | None:
    """Return the first scalar value matching a list of alias keys."""
    for alias in aliases:
        values = scalars.get(alias)
        if not values:
            continue
        for value in values:
            coerced = coerce_float(value)
            if coerced is not None:
                return coerced
    return None


def format_percent(value: float) -> str:
    """Format BUSCO percentages in a stable one-decimal representation."""
    return f"{value:.1f}"


def format_busco_string(percentages: dict[str, float], n_value: float) -> str:
    """Format BUSCO percentage metrics into the compact pipeline string."""
    return (
        "C:"
        f"{format_percent(percentages['C'])}%"
        "[S:"
        f"{format_percent(percentages['S'])}%,D:"
        f"{format_percent(percentages['D'])}%],F:"
        f"{format_percent(percentages['F'])}%,M:"
        f"{format_percent(percentages['M'])}%,n:{round(n_value)}"
    )


def build_compact_percentage_string(scalars: dict[str, list[Any]]) -> str | None:
    """Build a BUSCO string from direct C/S/D/F/M percentage fields when present."""
    compact_values = {
        key: find_scalar_value(scalars, (key.lower(),))
        for key in COMPACT_PERCENTAGE_KEYS
    }
    if any(value is None for value in compact_values.values()):
        return None

    n_value = find_scalar_value(scalars, COUNT_ALIASES["n"])
    if n_value is None:
        return None

    percentages = {
        key: value for key, value in compact_values.items() if value is not None
    }
    if any(value < 0 or value > 100 for value in percentages.values()):
        return None

    complete = percentages["C"]
    single = percentages["S"]
    duplicated = percentages["D"]
    fragmented = percentages["F"]
    missing = percentages["M"]
    if abs(complete - (single + duplicated)) > 0.05:
        return None
    if abs((complete + fragmented + missing) - 100.0) > 0.05:
        return None

    return format_busco_string(percentages, n_value)


def build_busco_string_from_scalars(scalars: dict[str, list[Any]]) -> str:
    """Construct a BUSCO summary string from numeric JSON fields."""
    compact_summary = build_compact_percentage_string(scalars)
    if compact_summary is not None:
        return compact_summary

    percentages = {
        key: find_scalar_value(scalars, aliases)
        for key, aliases in PERCENT_ALIASES.items()
    }
    if all(value is not None for value in percentages.values()):
        assert percentages["C"] is not None
        assert percentages["S"] is not None
        assert percentages["D"] is not None
        assert percentages["F"] is not None
        assert percentages["M"] is not None
        n_value = find_scalar_value(scalars, COUNT_ALIASES["n"])
        if n_value is None:
            raise ValueError("BUSCO JSON is missing the total marker count.")
        return format_busco_string(cast(dict[str, float], percentages), n_value)

    counts = {
        key: find_scalar_value(scalars, aliases)
        for key, aliases in COUNT_ALIASES.items()
    }
    if any(counts[key] is None for key in ("C", "S", "D", "F", "M", "n")):
        raise ValueError("BUSCO JSON is missing required summary fields.")

    total = counts["n"]
    assert total is not None
    if total <= 0:
        raise ValueError("BUSCO JSON reported a non-positive total marker count.")

    complete = counts["C"]
    single = counts["S"]
    duplicated = counts["D"]
    fragmented = counts["F"]
    missing = counts["M"]
    assert complete is not None
    assert single is not None
    assert duplicated is not None
    assert fragmented is not None
    assert missing is not None
    return (
        "C:"
        f"{format_percent(100 * complete / total)}%"
        "[S:"
        f"{format_percent(100 * single / total)}%,D:"
        f"{format_percent(100 * duplicated / total)}%"
        "],F:"
        f"{format_percent(100 * fragmented / total)}%,M:"
        f"{format_percent(100 * missing / total)}%,n:{round(total)}"
    )


def parse_summary(path: Path) -> str:
    """Parse a BUSCO JSON summary into the compact one-line summary string."""
    if not path.is_file():
        raise ValueError(f"Missing BUSCO summary file: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"Could not parse BUSCO JSON summary: {path}") from error

    existing_summary = find_existing_busco_string(payload)
    if existing_summary is not None:
        return existing_summary.replace(" ", "")

    scalars = collect_scalar_values(payload)
    return build_busco_string_from_scalars(scalars)


def write_output(
    accession: str,
    lineage: str,
    busco_summary: str,
    busco_status: str,
    warnings: str,
    output: Path,
) -> None:
    """Write the one-row BUSCO summary TSV."""
    busco_column = busco_column_name(lineage)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "accession",
                "lineage",
                busco_column,
                "busco_status",
                "warnings",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "accession": accession,
                "lineage": lineage,
                busco_column: busco_summary,
                "busco_status": busco_status,
                "warnings": warnings,
            }
        )


def main(argv: Sequence[str] | None = None) -> int:
    """Run the BUSCO summarisation CLI."""
    args = parse_args(argv)
    configure_logging()

    try:
        busco_summary = parse_summary(args.summary)
        busco_status = "done"
        warnings = ""
    except ValueError as error:
        LOGGER.warning(str(error))
        busco_summary = "NA"
        busco_status = "failed"
        warnings = "busco_summary_failed"

    write_output(
        accession=args.accession,
        lineage=args.lineage,
        busco_summary=busco_summary,
        busco_status=busco_status,
        warnings=warnings,
        output=args.output,
    )
    LOGGER.info("Wrote BUSCO summary for %s (%s).", args.accession, args.lineage)
    return 0


if __name__ == "__main__":
    sys.exit(main())
