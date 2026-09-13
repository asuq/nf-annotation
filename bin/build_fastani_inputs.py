#!/usr/bin/env python3
"""Build ANI eligibility tables and canonical FastANI input files."""

from __future__ import annotations

import argparse
import csv
import logging
import shutil
import sys
from pathlib import Path
from typing import Sequence

from ani_common import ANI_16S_POLICIES, derive_sixteen_s_ani_exclusion_reason
from atypical_warnings import classify_atypical_warnings
from build_master_table import (
    choose_assembly_level,
    detect_key_column,
    find_column_by_normalised_name,
    is_missing,
    load_assembly_stats_index,
    read_table,
    resolve_assembly_metric_value,
)


LOGGER = logging.getLogger(__name__)
ANI_METADATA_COLUMNS = (
    "accession",
    "matrix_name",
    "path",
    "assembly_level",
    "gcode",
    "checkm2_completeness",
    "checkm2_contamination",
    "n50",
    "scaffolds",
    "genome_size",
    "organism_name",
)
ANI_EXCLUSION_COLUMNS = (
    "accession",
    "internal_id",
    "ani_included",
    "ani_exclusion_reason",
)
MISSING_VALUE = "NA"


class FastAniInputError(RuntimeError):
    """Raised when ANI-input preparation cannot proceed."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build ANI eligibility metadata and canonical FastANI input files."
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
        "--staged-manifest",
        required=True,
        type=Path,
        help="Path to the staged genome manifest TSV.",
    )
    parser.add_argument(
        "--checkm2",
        required=True,
        type=Path,
        help="Path to the combined CheckM2 summary TSV.",
    )
    parser.add_argument(
        "--16s-status",
        dest="sixteen_s_status",
        required=True,
        type=Path,
        help="Path to the combined 16S status TSV.",
    )
    parser.add_argument(
        "--busco",
        action="append",
        required=True,
        type=Path,
        help="Path to a BUSCO summary TSV. May be supplied multiple times.",
    )
    parser.add_argument(
        "--primary-busco-column",
        required=True,
        help="Primary BUSCO column name, derived from the first configured lineage.",
    )
    parser.add_argument(
        "--assembly-stats",
        type=Path,
        help="Optional in-house assembly stats TSV keyed by accession.",
    )
    parser.add_argument(
        "--ani-16s-policy",
        choices=ANI_16S_POLICIES,
        default="complete",
        help="Require complete 16S, allow partial 16S, or ignore the 16S eligibility gate.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Directory for FastANI input outputs.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def join_tokens(tokens: Sequence[str]) -> str:
    """Join unique non-empty tokens in stable order."""
    return ";".join(dict.fromkeys(token for token in tokens if token))


def load_index(path: Path, key_aliases: Sequence[str]) -> dict[str, dict[str, str]]:
    """Load a TSV/CSV table into a unique accession-keyed index."""
    header, rows = read_table(path)
    key_column = detect_key_column(header, key_aliases)
    index: dict[str, dict[str, str]] = {}
    for row in rows:
        key = row.get(key_column, "").strip()
        if not key:
            raise FastAniInputError(f"Table {path} contains an empty key value.")
        if key in index:
            raise FastAniInputError(f"Table {path} contains duplicate key {key!r}.")
        index[key] = row
    return index


def load_staged_manifest(path: Path) -> dict[str, dict[str, str]]:
    """Load the staged genome manifest keyed by accession."""
    header, rows = read_table(path, delimiter="\t")
    required = ("accession", "internal_id", "staged_filename")
    missing = [column for column in required if column not in header]
    if missing:
        raise FastAniInputError(
            f"Staged genome manifest is missing required columns: {', '.join(missing)}"
        )

    manifest: dict[str, dict[str, str]] = {}
    for row in rows:
        accession = row.get("accession", "").strip()
        if not accession:
            raise FastAniInputError("Staged genome manifest contains an empty accession.")
        if accession in manifest:
            raise FastAniInputError(
                f"Staged genome manifest contains duplicate accession {accession!r}."
            )
        manifest[accession] = row
    return manifest


def load_busco_index(paths: Sequence[Path]) -> dict[str, dict[str, str]]:
    """Load one or more BUSCO summary TSVs keyed by accession."""
    index: dict[str, dict[str, str]] = {}
    seen_busco_columns: set[str] = set()

    for path in paths:
        header, rows = read_table(path, delimiter="\t")
        accession_column = detect_key_column(header, ("accession",))
        busco_columns = [column for column in header if column.startswith("BUSCO_")]
        if len(busco_columns) != 1:
            raise FastAniInputError(
                f"BUSCO table must contain exactly one BUSCO_<lineage> column: {path}"
            )
        busco_column = busco_columns[0]
        if busco_column in seen_busco_columns:
            raise FastAniInputError(f"Duplicate BUSCO lineage column supplied: {busco_column}")
        seen_busco_columns.add(busco_column)

        for row in rows:
            accession = row.get(accession_column, "").strip()
            if not accession:
                raise FastAniInputError(f"BUSCO table {path} contains an empty accession.")
            entry = index.setdefault(accession, {})
            value = row.get(busco_column, "").strip()
            entry[busco_column] = value if value else MISSING_VALUE
    return index


def load_metadata_index(path: Path) -> tuple[list[str], dict[str, dict[str, str]]]:
    """Load the metadata table keyed by accession."""
    header, rows = read_table(path)
    key_column = detect_key_column(header, ("Accession", "accession"))
    metadata_index: dict[str, dict[str, str]] = {}
    for row in rows:
        accession = row.get(key_column, "").strip()
        if not accession:
            raise FastAniInputError("Metadata table contains an empty accession.")
        if accession in metadata_index:
            raise FastAniInputError(
                f"Metadata table contains duplicate accession {accession!r}."
            )
        metadata_index[accession] = row
    return header, metadata_index


def detect_metadata_value(metadata_row: dict[str, str], column_name: str) -> str:
    """Return one metadata value by normalized column name or `NA`."""
    column = find_column_by_normalised_name(tuple(metadata_row), column_name)
    if column is None:
        return MISSING_VALUE
    value = metadata_row.get(column, "").strip()
    return MISSING_VALUE if is_missing(value) else value


def choose_checkm2_fields(sample_row: dict[str, str]) -> tuple[str, str, str]:
    """Return gcode-specific completeness and contamination values."""
    gcode = sample_row.get("Gcode", MISSING_VALUE) or MISSING_VALUE
    if gcode == "4":
        return (
            gcode,
            sample_row.get("Completeness_gcode4", MISSING_VALUE) or MISSING_VALUE,
            sample_row.get("Contamination_gcode4", MISSING_VALUE) or MISSING_VALUE,
        )
    if gcode == "11":
        return (
            gcode,
            sample_row.get("Completeness_gcode11", MISSING_VALUE) or MISSING_VALUE,
            sample_row.get("Contamination_gcode11", MISSING_VALUE) or MISSING_VALUE,
        )
    return gcode, MISSING_VALUE, MISSING_VALUE


def ensure_fastani_input(
    staged_filename: str,
    manifest_dir: Path,
    internal_id: str,
    output_dir: Path,
) -> tuple[str, str]:
    """Materialise one canonical FastANI input FASTA and return its path strings."""
    source = Path(staged_filename)
    if not source.is_absolute():
        source = manifest_dir / source
    if not source.is_file():
        raise FastAniInputError(f"Missing staged genome in build_fastani_inputs: {staged_filename}")

    output_dir.mkdir(parents=True, exist_ok=True)
    target_name = f"{internal_id}.fasta"
    target = output_dir / target_name
    if target.exists() or target.is_symlink():
        target.unlink()

    source_real = source.resolve()
    try:
        target.hardlink_to(source_real)
    except OSError:
        shutil.copy2(source_real, target)

    relative_path = str(Path(output_dir.name) / target_name)
    return relative_path, str(target)


def write_tsv(path: Path, header: Sequence[str], rows: Sequence[dict[str, str]]) -> None:
    """Write rows to a TSV file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(header),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_path_list(path: Path, paths: Sequence[str]) -> None:
    """Write one FastANI path per line."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(paths) + ("\n" if paths else ""), encoding="utf-8")


def run_build_fastani_inputs(
    validated_samples: Path,
    metadata: Path,
    staged_manifest: Path,
    checkm2: Path,
    sixteen_s_status: Path,
    busco: Sequence[Path],
    primary_busco_column: str,
    assembly_stats: Path | None,
    outdir: Path,
    ani_16s_policy: str = "complete",
) -> None:
    """Build FastANI path lists, ANI metadata, and ANI exclusion rows."""
    validated_header, validated_rows = read_table(validated_samples, delimiter="\t")
    if "accession" not in validated_header or "internal_id" not in validated_header:
        raise FastAniInputError(
            "validated_samples.tsv is missing accession/internal_id columns."
        )

    metadata_header, metadata_index = load_metadata_index(metadata)
    checkm2_index = load_index(checkm2, ("accession",))
    sixteen_s_index = load_index(sixteen_s_status, ("accession",))
    busco_index = load_busco_index(busco)
    staged_index = load_staged_manifest(staged_manifest)
    assembly_stats_index = load_assembly_stats_index(
        assembly_stats,
        {row["accession"] for row in validated_rows},
    )

    if primary_busco_column not in {
        column for busco_row in busco_index.values() for column in busco_row
    }:
        raise FastAniInputError(
            f"Primary BUSCO column {primary_busco_column!r} is not present in BUSCO summaries."
        )

    fastani_inputs_dir = outdir / "fastani_inputs"
    fastani_inputs_dir.mkdir(parents=True, exist_ok=True)
    metadata_rows: list[dict[str, str]] = []
    exclusion_rows: list[dict[str, str]] = []
    fastani_paths: list[str] = []

    for sample_row in validated_rows:
        accession = sample_row.get("accession", "").strip()
        internal_id = sample_row.get("internal_id", "").strip()
        metadata_row = metadata_index.get(accession, {})
        checkm2_row = checkm2_index.get(accession, {})
        sixteen_s_row = sixteen_s_index.get(accession, {})
        busco_row = busco_index.get(accession, {})
        staged_row = staged_index.get(accession, {})
        assembly_stats_row = assembly_stats_index.get(accession)

        reasons: list[str] = []
        gcode, checkm2_completeness, checkm2_contamination = choose_checkm2_fields(
            checkm2_row
        )
        if gcode not in {"4", "11"}:
            reasons.append("gcode_na")

        low_quality = checkm2_row.get("Low_quality", MISSING_VALUE) or MISSING_VALUE
        if low_quality == "true":
            reasons.append("low_quality")
        elif low_quality == "NA" and gcode in {"4", "11"}:
            reasons.append("missing_low_quality")

        sixteen_s_reason = derive_sixteen_s_ani_exclusion_reason(
            sixteen_s_row.get("16S", MISSING_VALUE) or MISSING_VALUE,
            policy=ani_16s_policy,
        )
        if sixteen_s_reason is not None:
            reasons.append(sixteen_s_reason)

        atypical_value = (
            detect_metadata_value(metadata_row, "Atypical_Warnings")
            if metadata_row
            else MISSING_VALUE
        )
        is_atypical, is_exception = classify_atypical_warnings(atypical_value)
        if is_atypical and not is_exception:
            reasons.append("atypical")

        primary_busco_value = busco_row.get(primary_busco_column, MISSING_VALUE)
        if is_missing(primary_busco_value):
            reasons.append("missing_primary_busco")

        assembly_level = choose_assembly_level(sample_row, metadata_row)
        if is_missing(assembly_level):
            reasons.append("missing_assembly_level")

        n50 = (
            resolve_assembly_metric_value(metadata_row, assembly_stats_row, "N50")
            if metadata_row or assembly_stats_row
            else MISSING_VALUE
        )
        if is_missing(n50):
            reasons.append("missing_n50")

        scaffolds = (
            resolve_assembly_metric_value(metadata_row, assembly_stats_row, "Scaffolds")
            if metadata_row or assembly_stats_row
            else MISSING_VALUE
        )
        if is_missing(scaffolds):
            reasons.append("missing_scaffolds")

        genome_size = (
            resolve_assembly_metric_value(
                metadata_row, assembly_stats_row, "Genome_Size"
            )
            if metadata_row or assembly_stats_row
            else MISSING_VALUE
        )
        if is_missing(genome_size):
            reasons.append("missing_genome_size")

        if is_missing(checkm2_completeness):
            reasons.append("missing_checkm2_completeness")
        if is_missing(checkm2_contamination):
            reasons.append("missing_checkm2_contamination")

        staged_filename = staged_row.get("staged_filename", "").strip()
        if not staged_filename:
            reasons.append("missing_staged_genome")

        inclusion = "false"
        relative_path = MISSING_VALUE
        if not reasons:
            relative_path, _absolute_path = ensure_fastani_input(
                staged_filename=staged_filename,
                manifest_dir=staged_manifest.parent,
                internal_id=internal_id,
                output_dir=fastani_inputs_dir,
            )
            metadata_record = {
                "accession": accession,
                "matrix_name": relative_path,
                "path": relative_path,
                "assembly_level": assembly_level,
                "gcode": gcode,
                "checkm2_completeness": checkm2_completeness,
                "checkm2_contamination": checkm2_contamination,
                "n50": n50,
                "scaffolds": scaffolds,
                "genome_size": genome_size,
                "organism_name": (
                    detect_metadata_value(metadata_row, "Organism_Name")
                    if metadata_row
                    else MISSING_VALUE
                ),
                primary_busco_column: primary_busco_value,
            }
            metadata_rows.append(metadata_record)
            fastani_paths.append(relative_path)
            inclusion = "true"

        exclusion_rows.append(
            {
                "accession": accession,
                "internal_id": internal_id,
                "ani_included": inclusion,
                "ani_exclusion_reason": join_tokens(reasons),
            }
        )

    metadata_header_out = [*ANI_METADATA_COLUMNS, primary_busco_column]
    write_path_list(outdir / "fastani_paths.txt", fastani_paths)
    write_tsv(outdir / "ani_metadata.tsv", metadata_header_out, metadata_rows)
    write_tsv(outdir / "ani_exclusions.tsv", ANI_EXCLUSION_COLUMNS, exclusion_rows)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the FastANI-input preparation CLI."""
    args = parse_args(argv)
    configure_logging()
    try:
        run_build_fastani_inputs(
            validated_samples=args.validated_samples,
            metadata=args.metadata,
            staged_manifest=args.staged_manifest,
            checkm2=args.checkm2,
            sixteen_s_status=args.sixteen_s_status,
            busco=args.busco,
            primary_busco_column=args.primary_busco_column,
            assembly_stats=args.assembly_stats,
            outdir=args.outdir,
            ani_16s_policy=args.ani_16s_policy,
        )
    except (FastAniInputError, FileNotFoundError, OSError) as error:
        LOGGER.error(str(error))
        return 1

    LOGGER.info("Built FastANI input paths and ANI-ready metadata.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
