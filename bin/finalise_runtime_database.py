#!/usr/bin/env python3
"""Validate runtime database directories and write ready markers."""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

from prepare_runtime_databases import (
    PreparationRecord,
    PrepareRuntimeDatabasesError,
    assess_destination,
    build_busco_validator,
    describe_validation,
    normalise_path,
    validate_checkm2,
    validate_codetta,
    validate_eggnog,
    validate_padloc,
    write_marker,
    write_report,
)

LOGGER = logging.getLogger(__name__)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Validate one prepared runtime database directory, optionally write "
            "its ready marker, and emit one TSV report row."
        )
    )
    parser.add_argument(
        "--component",
        choices=("checkm2", "busco_root", "codetta", "eggnog", "padloc"),
        required=True,
        help="Database component to validate.",
    )
    parser.add_argument(
        "--destination",
        type=Path,
        required=True,
        help="Database destination directory to validate.",
    )
    parser.add_argument(
        "--report",
        type=Path,
        required=True,
        help="One-row TSV report path.",
    )
    parser.add_argument(
        "--mode",
        choices=("reuse", "write", "prepared"),
        required=True,
        help="Reuse an existing ready marker or validate and write a new marker.",
    )
    parser.add_argument(
        "--source",
        default="NA",
        help="Source string recorded when writing a new marker.",
    )
    parser.add_argument(
        "--busco-lineage",
        action="append",
        default=[],
        help="Configured BUSCO lineage. May be supplied multiple times.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def build_validator(component: str, busco_lineages: Sequence[str]):
    """Return the validator for one runtime database component."""
    if component == "checkm2":
        return validate_checkm2
    if component == "busco_root":
        if not busco_lineages:
            raise PrepareRuntimeDatabasesError(
                "BUSCO finalisation requires at least one --busco-lineage value."
            )
        return build_busco_validator(busco_lineages)
    if component == "codetta":
        return validate_codetta
    if component == "eggnog":
        return validate_eggnog
    if component == "padloc":
        return validate_padloc
    raise PrepareRuntimeDatabasesError(f"Unsupported component: {component}")


def finalise_runtime_database(args: argparse.Namespace) -> PreparationRecord:
    """Validate one database destination and emit a report record."""
    destination = normalise_path(args.destination)
    validator = build_validator(args.component, args.busco_lineage)
    mode = "write" if args.mode == "prepared" else args.mode

    if mode == "reuse":
        _state, validation, marker = assess_destination(
            component=args.component,
            destination=destination,
            validator=validator,
        )
        if validation is None:
            raise PrepareRuntimeDatabasesError(
                "Reuse requires a valid prepared database"
            )
        return PreparationRecord(
            component=args.component,
            status="present",
            source=str(marker.get("source", "NA")) if marker is not None else "NA",
            destination=destination,
            details=describe_validation(validation),
        )

    validation = validator(destination)
    write_marker(
        component=args.component,
        source=args.source,
        destination=destination,
        validation=validation,
    )
    return PreparationRecord(
        component=args.component,
        status="prepared",
        source=args.source,
        destination=destination,
        details=describe_validation(validation, {"source_mode": "native_download"}),
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Run one runtime database finalisation step."""
    args = parse_args(argv)
    configure_logging()
    try:
        record = finalise_runtime_database(args)
        write_report(args.report, [record])
    except PrepareRuntimeDatabasesError as error:
        LOGGER.error("%s", error)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
