#!/usr/bin/env python3
"""Expand Tax_ID lineage columns from a user-supplied pinned NCBI taxdump."""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

from build_master_table import (
    detect_key_column,
    find_column_by_normalised_name,
    read_table,
)

LOGGER = logging.getLogger(__name__)
LINEAGE_COLUMNS = (
    "Tax_ID",
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)
TARGET_RANKS = frozenset(LINEAGE_COLUMNS[1:])


class TaxonomyExpandError(RuntimeError):
    """Raised when taxonomy expansion inputs are invalid."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Expand metadata Tax_ID values using a pinned NCBI taxdump."
    )
    parser.add_argument(
        "--validated-samples",
        required=True,
        type=Path,
        help="Path to validated_samples.tsv.",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="Path to the metadata table (CSV or TSV).",
    )
    parser.add_argument(
        "--taxdump",
        required=True,
        type=Path,
        help="Directory containing names.dmp and nodes.dmp.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Path to taxonomy_expanded.tsv.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def split_dmp_fields(line: str) -> list[str]:
    """Split one NCBI `.dmp` line into stripped fields."""
    fields = [field.strip() for field in line.rstrip("\n").split("\t|\t")]
    if fields and fields[-1].endswith("\t|"):
        fields[-1] = fields[-1][:-2].strip()
    elif fields and fields[-1].endswith("|"):
        fields[-1] = fields[-1][:-1].strip()
    return fields


def load_requested_tax_ids(validated_samples: Path, metadata: Path) -> list[str]:
    """Return unique metadata Tax_ID values referenced by requested samples."""
    _sample_header, sample_rows = read_table(validated_samples, delimiter="\t")
    metadata_header, metadata_rows = read_table(metadata)
    metadata_key = detect_key_column(metadata_header, ("Accession", "accession"))
    tax_id_column = find_column_by_normalised_name(metadata_header, "tax_id")
    if tax_id_column is None:
        return []

    requested_accessions = [row.get("accession", "").strip() for row in sample_rows]
    metadata_index = {row[metadata_key].strip(): row for row in metadata_rows}

    requested_tax_ids: list[str] = []
    seen_tax_ids: set[str] = set()
    for accession in requested_accessions:
        metadata_row = metadata_index.get(accession)
        if metadata_row is None:
            continue
        tax_id = metadata_row.get(tax_id_column, "").strip()
        if not tax_id or tax_id.upper() == "NA" or tax_id in seen_tax_ids:
            continue
        requested_tax_ids.append(tax_id)
        seen_tax_ids.add(tax_id)
    return requested_tax_ids


def load_nodes(path: Path) -> dict[str, tuple[str, str]]:
    """Load taxon parent and rank records from `nodes.dmp`."""
    if not path.is_file():
        raise TaxonomyExpandError(f"Missing taxdump nodes file: {path}")

    nodes: dict[str, tuple[str, str]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = split_dmp_fields(line)
            if len(fields) < 3:
                raise TaxonomyExpandError(f"Malformed nodes.dmp line: {line.rstrip()}")
            tax_id, parent_tax_id, rank = fields[:3]
            nodes[tax_id] = (parent_tax_id, rank)
    return nodes


def load_scientific_names(path: Path) -> dict[str, str]:
    """Load scientific taxon names from `names.dmp`."""
    if not path.is_file():
        raise TaxonomyExpandError(f"Missing taxdump names file: {path}")

    names: dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = split_dmp_fields(line)
            if len(fields) < 4:
                raise TaxonomyExpandError(f"Malformed names.dmp line: {line.rstrip()}")
            tax_id, name_txt, _unique_name, name_class = fields[:4]
            if name_class == "scientific name" and tax_id not in names:
                names[tax_id] = name_txt
    return names


def expand_tax_id(
    tax_id: str,
    nodes: dict[str, tuple[str, str]],
    names: dict[str, str],
) -> dict[str, str] | None:
    """Expand one Tax_ID to the required lineage columns or return `None`."""
    if tax_id not in nodes:
        return None

    lineage = {column: "NA" for column in LINEAGE_COLUMNS}
    lineage["Tax_ID"] = tax_id
    current_tax_id = tax_id
    visited: set[str] = set()

    while current_tax_id in nodes and current_tax_id not in visited:
        visited.add(current_tax_id)
        parent_tax_id, rank = nodes[current_tax_id]
        if rank in TARGET_RANKS and lineage[rank] == "NA":
            lineage[rank] = names.get(current_tax_id, "NA")
        if parent_tax_id == current_tax_id:
            break
        current_tax_id = parent_tax_id

    return lineage


def write_tsv(path: Path, rows: Sequence[dict[str, str]]) -> None:
    """Write expanded taxonomy rows to a TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(LINEAGE_COLUMNS), delimiter="\t"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_taxonomy_expand(
    validated_samples: Path,
    metadata: Path,
    taxdump: Path,
    output: Path,
) -> None:
    """Expand requested metadata Tax_ID values into lineage columns."""
    requested_tax_ids = load_requested_tax_ids(validated_samples, metadata)
    nodes = load_nodes(taxdump / "nodes.dmp")
    names = load_scientific_names(taxdump / "names.dmp")

    expanded_rows: list[dict[str, str]] = []
    missing_tax_ids: list[str] = []
    for tax_id in requested_tax_ids:
        lineage = expand_tax_id(tax_id, nodes, names)
        if lineage is None:
            missing_tax_ids.append(tax_id)
            continue
        expanded_rows.append(lineage)

    if missing_tax_ids:
        LOGGER.warning(
            "Skipped %d metadata Tax_ID value(s) missing from the pinned taxdump.",
            len(missing_tax_ids),
        )

    write_tsv(output, expanded_rows)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the taxonomy expansion CLI."""
    args = parse_args(argv)
    configure_logging()
    try:
        run_taxonomy_expand(
            validated_samples=args.validated_samples,
            metadata=args.metadata,
            taxdump=args.taxdump,
            output=args.output,
        )
    except (TaxonomyExpandError, FileNotFoundError, OSError) as error:
        LOGGER.error(str(error))
        return 1

    LOGGER.info("Expanded taxonomy for requested metadata Tax_ID values.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
