#!/usr/bin/env python3
"""Collect tool, container, and resource provenance into a final TSV report."""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Sequence


LOGGER = logging.getLogger(__name__)
OUTPUT_COLUMNS = ("component", "kind", "version", "image_or_path", "notes")
VERSION_KIND_KEYS = {
    "python": "runtime",
    "script": "pipeline",
    "transform": "pipeline",
}


class CollectVersionsError(RuntimeError):
    """Raised when version collection inputs are invalid."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Collect tool, container, and database provenance into a TSV."
    )
    parser.add_argument(
        "--version-file",
        action="append",
        default=[],
        type=Path,
        help="Path to one process versions.yml file. May be supplied multiple times.",
    )
    parser.add_argument(
        "--version-dir",
        action="append",
        default=[],
        type=Path,
        help="Directory containing staged process versions.yml files.",
    )
    parser.add_argument(
        "--inherited-version-dir",
        type=Path,
        help="Directory of provenance TSVs from explicitly reused published results.",
    )
    parser.add_argument(
        "--container-ref",
        action="append",
        default=[],
        help="Container reference as name=value.",
    )
    parser.add_argument(
        "--busco-lineage",
        action="append",
        default=[],
        help="Configured BUSCO lineage name. May be supplied multiple times.",
    )
    parser.add_argument(
        "--nextflow-version",
        default="NA",
        help="Nextflow version string.",
    )
    parser.add_argument(
        "--pipeline-version",
        default="NA",
        help="Pipeline manifest version string.",
    )
    parser.add_argument(
        "--git-commit",
        default="NA",
        help="Pipeline git commit hash or NA.",
    )
    parser.add_argument(
        "--container-engine",
        default="NA",
        help="Container engine name or NA.",
    )
    parser.add_argument(
        "--checkm2-db",
        default="NA",
        help="CheckM2 database path.",
    )
    parser.add_argument(
        "--checkm2-db-label",
        default="NA",
        help="Optional CheckM2 database version or label.",
    )
    parser.add_argument(
        "--taxdump",
        default="NA",
        help="Pinned taxdump path.",
    )
    parser.add_argument(
        "--taxdump-label",
        default="NA",
        help="Optional pinned taxdump label.",
    )
    parser.add_argument(
        "--codetta-db",
        default="NA",
        help="Codetta profile database path.",
    )
    parser.add_argument(
        "--codetta-db-label",
        default="NA",
        help="Optional Codetta profile database label.",
    )
    parser.add_argument(
        "--busco-db",
        default="NA",
        help="BUSCO dataset root path.",
    )
    parser.add_argument(
        "--eggnog-db",
        default="NA",
        help="eggNOG database path.",
    )
    parser.add_argument(
        "--eggnog-db-label",
        default="NA",
        help="Optional eggNOG database label.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Path to tool_and_db_versions.tsv.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def normalise_value(value: str | None) -> str:
    """Normalise blank values to NA."""
    if value is None:
        return "NA"
    stripped = value.strip()
    return stripped if stripped else "NA"


def build_row(
    component: str,
    kind: str,
    version: str | None = "NA",
    image_or_path: str | None = "NA",
    notes: str | None = "",
) -> dict[str, str]:
    """Build one normalized provenance row."""
    return {
        "component": normalise_value(component),
        "kind": normalise_value(kind),
        "version": normalise_value(version),
        "image_or_path": normalise_value(image_or_path),
        "notes": notes.strip() if notes else "",
    }


def parse_name_value(token: str, argument_name: str) -> tuple[str, str]:
    """Parse one name=value token."""
    if "=" not in token:
        raise CollectVersionsError(
            f"{argument_name} entry must use name=value syntax: {token!r}"
        )
    name, value = token.split("=", 1)
    name = name.strip()
    if not name:
        raise CollectVersionsError(
            f"{argument_name} entry is missing the name before '=': {token!r}"
        )
    return name, normalise_value(value)


def parse_versions_header(line: str, path: Path) -> str:
    """Parse one non-indented versions.yml header line."""
    stripped = line.strip()
    if not stripped.endswith(":"):
        raise CollectVersionsError(
            f"Malformed versions.yml header line in {path}: {stripped!r}"
        )
    return stripped[:-1].strip().strip("'\"")


def parse_versions_entry(line: str, path: Path) -> tuple[str, str]:
    """Parse one indented versions.yml entry line."""
    if not line.startswith("  ") or line.startswith("   ") or "\t" in line[:2]:
        raise CollectVersionsError(
            f"Malformed versions.yml indentation in {path}: {line.rstrip()!r}"
        )
    stripped = line.strip()
    if ":" not in stripped:
        raise CollectVersionsError(
            f"Malformed versions.yml entry in {path}: {stripped!r}"
        )
    name, value = stripped.split(":", 1)
    raw_value = value.strip()
    if not raw_value:
        raise CollectVersionsError(
            f"Missing versions.yml value in {path}: {stripped!r}"
        )
    if raw_value[:1] in {'"', "'"} and not raw_value.endswith(raw_value[:1]):
        raise CollectVersionsError(
            f"Multiline quoted values are not supported in {path}: {stripped!r}"
        )
    return name.strip(), raw_value.strip().strip("'\"")


def parse_versions_sections(path: Path) -> dict[str, list[dict[str, str]]]:
    """Parse one versions file into per-process row groups."""
    if not path.is_file():
        raise CollectVersionsError(f"Missing versions file: {path}")

    sections: dict[str, list[dict[str, str]]] = {}
    current_source: str | None = None
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            stripped = line.strip()
            if not stripped or stripped == "EOF":
                continue
            if not line.startswith((" ", "\t")):
                current_source = parse_versions_header(line, path)
                sections.setdefault(current_source, [])
                continue
            if current_source is None:
                raise CollectVersionsError(
                    f"Malformed versions.yml entry without a header in {path}: {line.rstrip()!r}"
                )
            clean_name, clean_value = parse_versions_entry(line, path)
            kind = VERSION_KIND_KEYS.get(clean_name, "tool")
            if clean_name == "script":
                sections[current_source].append(
                    build_row(
                        component=current_source,
                        kind=kind,
                        version="NA",
                        image_or_path=clean_value,
                        notes="script",
                    )
                )
                continue
            if clean_name == "transform":
                sections[current_source].append(
                    build_row(
                        component=current_source,
                        kind=kind,
                        version="NA",
                        image_or_path=clean_value,
                        notes="helper",
                    )
                )
                continue
            sections[current_source].append(
                build_row(
                    component=clean_name,
                    kind=kind,
                    version=clean_value,
                    notes=f"reported by {current_source}",
                )
            )
    return sections


def parse_versions_file(path: Path) -> list[dict[str, str]]:
    """Parse one simple Nextflow versions.yml file into row dictionaries."""
    rows: list[dict[str, str]] = []
    for section_rows in parse_versions_sections(path).values():
        rows.extend(section_rows)
    return rows


def build_runtime_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    """Build runtime and resource report rows from CLI arguments."""
    rows = [
        build_row("nextflow", "runtime", args.nextflow_version, notes="workflow"),
        build_row("pipeline", "pipeline", args.pipeline_version, notes="workflow manifest"),
        build_row("pipeline_git_commit", "pipeline", args.git_commit, notes="workflow revision"),
        build_row("container_engine", "runtime", args.container_engine, notes="workflow"),
        build_row(
            "checkm2_db",
            "database",
            args.checkm2_db_label,
            args.checkm2_db,
            "CheckM2 database",
        ),
        build_row("taxdump", "database", args.taxdump_label, args.taxdump, "Pinned taxdump"),
        build_row(
            "codetta_db",
            "database",
            args.codetta_db_label,
            args.codetta_db,
            "Codetta profile database",
        ),
        build_row(
            "busco_datasets",
            "database",
            ";".join(args.busco_lineage) if args.busco_lineage else "NA",
            args.busco_db,
            "Configured BUSCO lineages in order",
        ),
        build_row(
            "eggnog_db",
            "database",
            args.eggnog_db_label,
            args.eggnog_db,
            "eggNOG database",
        ),
    ]
    for lineage in args.busco_lineage:
        rows.append(
            build_row(
                component=lineage,
                kind="database",
                version="NA",
                image_or_path=(
                    str(Path(args.busco_db) / lineage)
                    if normalise_value(args.busco_db) != "NA"
                    else "NA"
                ),
                notes="BUSCO lineage dataset",
            )
        )
    for token in args.container_ref:
        name, value = parse_name_value(token, "--container-ref")
        rows.append(
            build_row(name, "container", "NA", value, "params container reference")
        )
    return rows


def discover_version_files(
    version_files: Sequence[Path], version_dirs: Sequence[Path]
) -> list[Path]:
    """Collect direct version files and staged directory contents in stable order."""
    discovered = list(version_files)
    for version_dir in version_dirs:
        if not version_dir.exists():
            raise CollectVersionsError(f"Missing versions directory: {version_dir}")
        if not version_dir.is_dir():
            raise CollectVersionsError(f"Versions path is not a directory: {version_dir}")
        discovered.extend(sorted(version_dir.glob("*.yml")))
    return discovered


def build_process_signature(rows: Sequence[dict[str, str]]) -> tuple[tuple[str, ...], ...]:
    """Build a comparable signature for one process row set."""
    return tuple(
        (
            row["component"],
            row["kind"],
            row["version"],
            row["image_or_path"],
            row["notes"],
        )
        for row in rows
    )


def collect_canonical_version_rows(version_files: Sequence[Path]) -> list[dict[str, str]]:
    """Keep one canonical parsed row set per process header."""
    canonical_rows: dict[str, list[dict[str, str]]] = {}
    canonical_sources: dict[str, Path] = {}
    canonical_signatures: dict[str, tuple[tuple[str, ...], ...]] = {}

    for version_file in version_files:
        for process_name, rows in parse_versions_sections(version_file).items():
            signature = build_process_signature(rows)
            if process_name not in canonical_signatures:
                canonical_signatures[process_name] = signature
                canonical_sources[process_name] = version_file
                canonical_rows[process_name] = rows
                continue
            if canonical_signatures[process_name] != signature:
                raise CollectVersionsError(
                    "Conflicting versions.yml entries for process "
                    f"{process_name!r}: {canonical_sources[process_name]} != {version_file}"
                )

    rows: list[dict[str, str]] = []
    for process_name in sorted(canonical_rows):
        rows.extend(canonical_rows[process_name])
    return rows


def deduplicate_rows(rows: Sequence[dict[str, str]]) -> list[dict[str, str]]:
    """Deduplicate rows while keeping a stable, category-aware order."""
    seen: set[tuple[str, str, str, str, str]] = set()
    ordered_rows: list[dict[str, str]] = []
    kind_order = {
        "pipeline": 0,
        "runtime": 0,
        "database": 1,
        "container": 3,
        "tool": 4,
    }
    for row in sorted(
        rows,
        key=lambda row: (
            kind_order.get(row["kind"], 99),
            row["component"],
            row["notes"],
            row["version"],
            row["image_or_path"],
        ),
    ):
        key = (
            row["component"],
            row["kind"],
            row["version"],
            row["image_or_path"],
            row["notes"],
        )
        if key in seen:
            continue
        seen.add(key)
        ordered_rows.append(row)
    return ordered_rows


def write_tsv(path: Path, rows: Sequence[dict[str, str]]) -> None:
    """Write the final versions TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(OUTPUT_COLUMNS), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_collect_versions(args: argparse.Namespace) -> None:
    """Build the final tool and database versions TSV."""
    rows = build_runtime_rows(args)
    version_files = discover_version_files(args.version_file, args.version_dir)
    rows.extend(collect_canonical_version_rows(version_files))
    if args.inherited_version_dir is not None:
        paths = sorted(args.inherited_version_dir.glob("*.tsv"))
        if not paths:
            raise CollectVersionsError(
                f"No inherited provenance reports in {args.inherited_version_dir}"
            )
        for path in paths:
            with path.open(encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames != list(OUTPUT_COLUMNS):
                    raise CollectVersionsError(
                        f"Invalid inherited provenance header: {path}"
                    )
                for row in reader:
                    if None in row or any(value is None for value in row.values()):
                        raise CollectVersionsError(
                            f"Malformed inherited provenance row: {path}"
                        )
                    rows.append(row)
    write_tsv(args.output, deduplicate_rows(rows))


def main(argv: Sequence[str] | None = None) -> int:
    """Run the versions collector CLI."""
    args = parse_args(argv)
    configure_logging()
    try:
        run_collect_versions(args)
    except CollectVersionsError as error:
        LOGGER.error(str(error))
        return 1

    LOGGER.info("Wrote tool and database provenance to %s.", args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
