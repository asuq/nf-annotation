#!/usr/bin/env python3
"""Select ANI representatives and emit stable ANI summary tables."""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Sequence

import master_table_contract
from ani_common import (
    AniInputError,
    is_empty_ani_cohort,
    load_ani_metadata,
    load_matrix,
    normalize_header,
    subset_matrix,
)
from ani_representatives import build_ani_summary_from_clusters


LOGGER = logging.getLogger(__name__)
ANI_SUMMARY_COLUMNS = ("Accession", *master_table_contract.ANI_COLUMNS)
ANI_REPRESENTATIVE_COLUMNS = (
    "Cluster_ID",
    "Representative_Accession",
    "Organism_Name",
    "CheckM2_Completeness",
    "CheckM2_Contamination",
    "BUSCO",
    "Assembly_Level",
    "N50",
    "Cluster_Size",
)


class RepresentativeSelectionError(RuntimeError):
    """Raise when ANI representative selection inputs are inconsistent."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments for ANI representative selection."""
    parser = argparse.ArgumentParser(
        description="Select ANI representatives from cluster memberships and FastANI outputs."
    )
    parser.add_argument(
        "--ani-clusters",
        required=True,
        type=Path,
        help="Path to cluster.tsv from cluster_ani.py.",
    )
    parser.add_argument(
        "--ani-metadata",
        required=True,
        type=Path,
        help="Path to ani_metadata.tsv from build_fastani_inputs.py.",
    )
    parser.add_argument(
        "--ani-matrix",
        required=True,
        type=Path,
        help="Path to the FastANI PHYLIP matrix.",
    )
    parser.add_argument(
        "--ani-score-profile",
        choices=("default", "isolate", "mag"),
        default="default",
        help="Representative-selection score profile. Default: default.",
    )
    parser.add_argument(
        "--ani-summary-output",
        required=True,
        type=Path,
        help="Output path for ani_summary.tsv.",
    )
    parser.add_argument(
        "--ani-representatives-output",
        required=True,
        type=Path,
        help="Output path for ani_representatives.tsv.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV file into a header and row dictionaries."""
    if not path.is_file():
        raise RepresentativeSelectionError(f"Missing input table: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = [row for row in reader if any(cell.strip() for cell in row)]

    if not rows:
        raise RepresentativeSelectionError(f"Input table is empty: {path}")

    header = [column.strip() for column in rows[0]]
    if any(not column for column in header):
        raise RepresentativeSelectionError(f"Input table has an empty header name: {path}")
    if len(set(header)) != len(header):
        raise RepresentativeSelectionError(f"Input table contains duplicate headers: {path}")

    records: list[dict[str, str]] = []
    for line_number, values in enumerate(rows[1:], start=2):
        if len(values) != len(header):
            raise RepresentativeSelectionError(
                f"Input table {path} has {len(values)} columns on line {line_number}, "
                f"expected {len(header)}."
            )
        records.append(
            {column: value.strip() for column, value in zip(header, values, strict=True)}
        )
    return header, records


def detect_key_column(header: Sequence[str], aliases: Sequence[str]) -> str:
    """Detect one join-key column from a set of allowed aliases."""
    alias_set = {normalize_header(alias) for alias in aliases}
    matching_columns = [
        column for column in header if normalize_header(column) in alias_set
    ]
    if len(matching_columns) > 1:
        raise RepresentativeSelectionError(
            f"Input table contains multiple possible key columns: {matching_columns}"
        )
    if not matching_columns:
        raise RepresentativeSelectionError(
            f"Input table is missing a required key column. Expected one of: {aliases}"
        )
    return matching_columns[0]


def find_required_column(header: Sequence[str], normalized_name: str) -> str:
    """Resolve one required column by normalized header name."""
    matches = [
        column for column in header if normalize_header(column) == normalize_header(normalized_name)
    ]
    if len(matches) > 1:
        raise RepresentativeSelectionError(
            f"Input table contains multiple possible {normalized_name!r} columns: {matches}"
        )
    if not matches:
        raise RepresentativeSelectionError(
            f"Input table is missing required column: {normalized_name}"
        )
    return matches[0]


def build_index(
    rows: Sequence[dict[str, str]],
    key_column: str,
    table_name: str,
) -> dict[str, dict[str, str]]:
    """Build a unique keyed row index and fail on duplicates."""
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        key = row.get(key_column, "").strip()
        if not key:
            raise RepresentativeSelectionError(f"{table_name} contains an empty key value.")
        if key in indexed:
            raise RepresentativeSelectionError(f"{table_name} contains duplicate key {key!r}.")
        indexed[key] = row
    return indexed


def load_cluster_members(path: Path) -> dict[str, list[str]]:
    """Load ANI cluster memberships keyed by Cluster_ID with accession members."""
    header, rows = read_tsv(path)
    accession_column = detect_key_column(header, ("accession", "Accession"))
    cluster_column = find_required_column(header, "Cluster_ID")
    cluster_index = build_index(rows, accession_column, "ANI cluster table")

    cluster_members: dict[str, list[str]] = {}
    for accession, row in cluster_index.items():
        cluster_id = row.get(cluster_column, "").strip()
        if not cluster_id:
            raise RepresentativeSelectionError(
                f"ANI cluster table contains an empty Cluster_ID for accession {accession!r}."
            )
        cluster_members.setdefault(cluster_id, []).append(accession)
    return cluster_members


def build_ani_outputs(
    *,
    ani_clusters: Path,
    ani_metadata: Path,
    ani_matrix: Path,
    score_profile: str,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Build ANI summary and representative tables from raw ANI artifacts."""
    cluster_members_by_accession = load_cluster_members(ani_clusters)

    try:
        if is_empty_ani_cohort(ani_metadata, ani_matrix, require_scoring=True):
            if cluster_members_by_accession:
                raise RepresentativeSelectionError(
                    "ANI clusters contain samples but eligibility metadata is empty."
                )
            return [], []
        names, ani, name_to_idx = load_matrix(ani_matrix)
        metadata_by_matrix_name, eligible_names = load_ani_metadata(
            ani_metadata=ani_metadata,
            matrix_names=names,
            busco_column=None,
            matrix_name_column="matrix_name",
            require_existing_paths=False,
        )
        names, ani, name_to_idx = subset_matrix(names, ani, eligible_names, name_to_idx)
        metadata_by_matrix_name = {
            matrix_name: metadata_by_matrix_name[matrix_name]
            for matrix_name in eligible_names
        }
    except AniInputError as error:
        raise RepresentativeSelectionError(str(error)) from error

    matrix_name_by_accession: dict[str, str] = {}
    for matrix_name, genome in metadata_by_matrix_name.items():
        if genome.Accession in matrix_name_by_accession:
            raise RepresentativeSelectionError(
                f"ANI metadata contains duplicate accession {genome.Accession!r}."
            )
        matrix_name_by_accession[genome.Accession] = matrix_name

    cluster_members: dict[str, list[str]] = {}
    for cluster_id, accessions in cluster_members_by_accession.items():
        cluster_members[cluster_id] = []
        for accession in accessions:
            if accession not in matrix_name_by_accession:
                raise RepresentativeSelectionError(
                    "ANI cluster table contains accession not present in ANI metadata: "
                    f"{accession!r}."
                )
            cluster_members[cluster_id].append(matrix_name_by_accession[accession])

    cluster_accessions = {
        accession
        for accessions in cluster_members_by_accession.values()
        for accession in accessions
    }
    missing_cluster_accessions = sorted(set(matrix_name_by_accession) - cluster_accessions)
    if missing_cluster_accessions:
        raise RepresentativeSelectionError(
            "ANI metadata contains eligible accession(s) missing from ANI clusters: "
            + ", ".join(missing_cluster_accessions)
        )

    try:
        summary_index, representative_rows = build_ani_summary_from_clusters(
            cluster_members=cluster_members,
            names=names,
            ani=ani,
            name_to_idx=name_to_idx,
            meta=metadata_by_matrix_name,
            score_profile=score_profile,
            threads=1,
        )
    except RuntimeError as error:
        raise RepresentativeSelectionError(str(error)) from error

    summary_rows = [
        {"Accession": accession, **summary_index[accession]}
        for accession in sorted(summary_index)
    ]
    return summary_rows, representative_rows


def write_tsv(path: Path, header: Sequence[str], rows: Sequence[dict[str, str]]) -> None:
    """Write a TSV with a fixed header."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(header), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run(args: argparse.Namespace) -> None:
    """Run ANI representative selection and write both output tables."""
    summary_rows, representative_rows = build_ani_outputs(
        ani_clusters=args.ani_clusters,
        ani_metadata=args.ani_metadata,
        ani_matrix=args.ani_matrix,
        score_profile=args.ani_score_profile,
    )
    write_tsv(args.ani_summary_output, ANI_SUMMARY_COLUMNS, summary_rows)
    write_tsv(
        args.ani_representatives_output,
        ANI_REPRESENTATIVE_COLUMNS,
        representative_rows,
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Run the CLI entrypoint."""
    args = parse_args(argv)
    configure_logging()
    try:
        run(args)
    except RepresentativeSelectionError as error:
        LOGGER.error(str(error))
        return 1

    LOGGER.info("Wrote ANI summary to %s.", args.ani_summary_output)
    LOGGER.info("Wrote ANI representatives to %s.", args.ani_representatives_output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
