#!/usr/bin/env python3
"""Shared ANI parsing helpers for clustering and representative selection."""

from __future__ import annotations

import csv
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence


REQUIRED_SCORING_METADATA_COLS: set[str] = {
    "accession",
    "path",
    "assembly_level",
    "gcode",
    "checkm2_completeness",
    "checkm2_contamination",
    "n50",
    "scaffolds",
    "genome_size",
}
BUSCO_RE = re.compile(r"C:(?P<C>\d+(?:\.\d+)?)%.*?M:(?P<M>\d+(?:\.\d+)?)%", re.IGNORECASE)
ASSEMBLY_LEVEL_MAP: dict[str, str] = {
    "complete genome": "Complete Genome",
    "complete genomes": "Complete Genome",
    "chromosome": "Chromosome",
    "chromosomes": "Chromosome",
    "scaffold": "Scaffold",
    "scaffolds": "Scaffold",
    "contig": "Contig",
    "contigs": "Contig",
}
ASSEMBLY_RANK: dict[str, int] = {
    "Contig": 0,
    "Scaffold": 1,
    "Chromosome": 2,
    "Complete Genome": 3,
}
SIXTEEN_S_STATUS_MAP: dict[str, str] = {
    "yes": "Yes",
    "no": "No",
    "partial": "partial",
}


class AniInputError(RuntimeError):
    """Raised when ANI matrix or ANI metadata inputs are invalid."""


@dataclass(slots=True)
class Genome:
    """Store ANI representative-selection metadata for one eligible genome."""

    Accession: str
    Organism_Name: str
    Gcode: int
    CheckM2_Completeness: float
    CheckM2_Contamination: float
    N50: int
    Scaffolds: int
    Genome_Size: int
    BUSCO_str: str
    BUSCO_C: float
    BUSCO_M: float
    Assembly_Level: str
    Assembly_Rank: int
    Path: str
    Score: float | None = None


def normalize_header(name: str) -> str:
    """Normalize a column name into a lowercase underscore token."""
    return re.sub(r"[^a-z0-9]+", "_", name.strip().casefold()).strip("_")


def is_empty_ani_cohort(
    metadata_path: Path,
    matrix_path: Path,
    *,
    matrix_name_column: str = "matrix_name",
    require_scoring: bool = False,
) -> bool:
    """Recognise an explicitly empty cohort without accepting missing inputs.

    A header-only eligibility table and a zero-byte matrix are the pipeline's
    empty-cohort contract. Matrix data without eligible metadata is an error.
    """
    if not metadata_path.is_file() or not matrix_path.is_file():
        raise AniInputError("ANI metadata and matrix files must both exist.")
    with metadata_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = [normalize_header(column) for column in next(reader, [])]
        required = {"accession", normalize_header(matrix_name_column)}
        if require_scoring:
            required.update(
                normalize_header(column) for column in REQUIRED_SCORING_METADATA_COLS
            )
        if (
            not header
            or any(not column for column in header)
            or len(header) != len(set(header))
            or not required <= set(header)
            or (
                require_scoring
                and sum(column.startswith("busco_") for column in header) != 1
            )
        ):
            raise AniInputError(f"Invalid ANI metadata header: {metadata_path}")
        if any(row for row in reader):
            return False
    if matrix_path.stat().st_size:
        raise AniInputError(
            "ANI matrix is non-empty but eligibility metadata has no samples."
        )
    return True


def is_missing_value(value: Any) -> bool:
    """Return True when a scalar should be treated as missing."""
    return str(value).strip().upper() in {"", "NA"}


def normalise_sixteen_s_status(value: Any) -> str:
    """Normalise one 16S status value to the public reporting vocabulary."""
    text = str(value or "").strip()
    if is_missing_value(text):
        return "NA"
    return SIXTEEN_S_STATUS_MAP.get(text.casefold(), text)


def derive_sixteen_s_ani_exclusion_reason(
    value: Any,
    *,
    allow_incomplete: bool = False,
) -> str | None:
    """Return the ANI exclusion reason for one 16S status, if any."""
    status = normalise_sixteen_s_status(value)
    if status == "Yes":
        return None
    if status == "partial":
        return None if allow_incomplete else "partial_16s"
    if status == "No":
        return None if allow_incomplete else "no_16s"
    return "16s_na"


def try_parse_int_like(value: Any) -> int | None:
    """Parse an integer or integer-like value and return None on failure."""
    text = str(value).strip()
    if is_missing_value(text):
        return None
    if re.fullmatch(r"[+-]?\d+", text):
        return int(text)
    try:
        parsed = float(text)
    except Exception:
        return None
    if parsed.is_integer():
        return int(parsed)
    return None


def try_parse_float_like(value: Any) -> float | None:
    """Parse a float-like value and return None on failure."""
    text = str(value).strip()
    if is_missing_value(text):
        return None
    try:
        return float(text)
    except Exception:
        return None


def try_parse_busco(busco_str: Any) -> tuple[float, float] | None:
    """Extract BUSCO complete and missing percentages from the compact summary."""
    if not isinstance(busco_str, str) or is_missing_value(busco_str):
        return None
    match = BUSCO_RE.search(busco_str)
    if not match:
        return None
    return float(match.group("C")), float(match.group("M"))


def load_phylip_lower_triangular(path: Path) -> tuple[list[str], list[list[str]]]:
    """Load a PHYLIP lower-triangular matrix with strict structure checks."""
    if not path.is_file():
        raise AniInputError(f"ANI matrix file not found: {path}")

    lines: list[str] = []
    with path.open("rt", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                lines.append(stripped)

    if not lines:
        raise AniInputError("ANI matrix: empty file or missing taxon count.")
    try:
        sample_count = int(lines[0])
    except Exception as error:
        raise AniInputError(
            f"First ANI matrix line must be an integer taxon count: {lines[0]!r}"
        ) from error

    if sample_count == 0:
        raise AniInputError("ANI matrix does not contain any sample.")
    if len(lines) - 1 != sample_count:
        raise AniInputError(
            f"Expected {sample_count} ANI matrix rows, found {len(lines) - 1}."
        )

    names: list[str] = []
    rows: list[list[str]] = []
    for row_index in range(1, sample_count + 1):
        raw = lines[row_index]
        expected_values = row_index - 1
        parts = raw.rsplit(maxsplit=expected_values)
        if len(parts) != expected_values + 1:
            raise AniInputError(
                f"ANI matrix row {row_index} expected {expected_values} values: {raw!r}"
            )
        name = parts[0].strip()
        values = [part.strip() for part in parts[1:]]
        for column_index, token in enumerate(values, start=1):
            if token == "NA":
                continue
            try:
                parsed = float(token)
            except Exception as error:
                raise AniInputError(
                    f"ANI matrix token is not numeric or NA at row {name!r}, "
                    f"column {column_index}: {token!r}"
                ) from error
            if not 0.0 <= parsed <= 100.0:
                raise AniInputError(
                    f"ANI value out of range [0,100] at row {name!r}, "
                    f"column {column_index}: {parsed}"
                )
        names.append(name)
        rows.append(values)
    return names, rows


def load_matrix(path: Path) -> tuple[list[str], "np.ndarray", dict[str, int]]:
    """Convert a PHYLIP lower-triangular ANI file into a full symmetric matrix."""
    import numpy as np

    names, rows = load_phylip_lower_triangular(path)
    sample_count = len(names)
    ani = np.full((sample_count, sample_count), np.nan, dtype=np.float64)
    name_to_idx: dict[str, int] = {}

    for index, name in enumerate(names):
        if name in name_to_idx:
            raise AniInputError(f"Duplicate name in ANI matrix: {name!r}")
        name_to_idx[name] = index
        ani[index, index] = 100.0

    for row_index, values in enumerate(rows):
        for column_index, token in enumerate(values):
            if token == "NA":
                continue
            value = float(token)
            ani[row_index, column_index] = value
            ani[column_index, row_index] = value

    return names, ani, name_to_idx


def resolve_busco_column(header_map: dict[str, str], busco_column: str | None) -> str:
    """Resolve the normalized BUSCO column name from the ANI metadata header."""
    if busco_column is not None:
        normalized = normalize_header(busco_column)
        if normalized not in header_map:
            raise AniInputError(
                "ANI metadata TSV is missing the requested BUSCO column "
                f"{busco_column!r}."
            )
        return normalized

    busco_columns = [column for column in header_map if column.startswith("busco_")]
    if len(busco_columns) != 1:
        raise AniInputError(
            "ANI metadata TSV must contain exactly one BUSCO_<lineage> column "
            "when no BUSCO column is specified."
        )
    return busco_columns[0]


def build_genome_from_row(
    row: dict[str, str],
    *,
    acc: str,
    busco_column: str,
    require_existing_path: bool = True,
) -> tuple[Genome | None, list[str]]:
    """Convert one ANI metadata row into a Genome record or exclusion reasons."""
    reasons: list[str] = []

    path_str = row["path"].strip()
    if is_missing_value(path_str):
        reasons.append("missing_path")
    elif require_existing_path and not Path(path_str).exists():
        reasons.append("path_not_found")

    assembly_key = row["assembly_level"].strip().casefold()
    if assembly_key in ASSEMBLY_LEVEL_MAP:
        assembly_level = ASSEMBLY_LEVEL_MAP[assembly_key]
    else:
        assembly_level = ""
        reasons.append("invalid_assembly_level")

    gcode = try_parse_int_like(row["gcode"])
    if gcode not in (4, 11):
        reasons.append("invalid_gcode")

    checkm2_completeness = try_parse_float_like(row["checkm2_completeness"])
    if checkm2_completeness is None:
        reasons.append("missing_checkm2_completeness")

    checkm2_contamination = try_parse_float_like(row["checkm2_contamination"])
    if checkm2_contamination is None:
        reasons.append("missing_checkm2_contamination")

    n50 = try_parse_int_like(row["n50"])
    if n50 is None:
        reasons.append("missing_n50")

    scaffolds = try_parse_int_like(row["scaffolds"])
    if scaffolds is None:
        reasons.append("missing_scaffolds")

    genome_size = try_parse_int_like(row["genome_size"])
    if genome_size is None:
        reasons.append("missing_genome_size")

    busco_str = row[busco_column].strip()
    busco_parsed = try_parse_busco(busco_str)
    if busco_parsed is None:
        reasons.append(f"missing_{busco_column}")

    if reasons:
        return None, reasons

    assert gcode is not None
    assert checkm2_completeness is not None
    assert checkm2_contamination is not None
    assert n50 is not None
    assert scaffolds is not None
    assert genome_size is not None
    busco_c, busco_m = busco_parsed
    return (
        Genome(
            Accession=acc,
            Organism_Name=row.get("organism_name", "").strip(),
            Gcode=gcode,
            CheckM2_Completeness=checkm2_completeness,
            CheckM2_Contamination=checkm2_contamination,
            N50=n50,
            Scaffolds=scaffolds,
            Genome_Size=genome_size,
            BUSCO_str=busco_str,
            BUSCO_C=busco_c,
            BUSCO_M=busco_m,
            Assembly_Level=assembly_level,
            Assembly_Rank=ASSEMBLY_RANK[assembly_level],
            Path=path_str,
        ),
        [],
    )


def load_ani_metadata(
    ani_metadata: Path,
    matrix_names: Sequence[str],
    *,
    busco_column: str | None,
    matrix_name_column: str,
    require_existing_paths: bool = True,
) -> tuple[dict[str, Genome], list[str]]:
    """Load ANI scoring metadata and exclude rows missing representative inputs."""
    if not ani_metadata.is_file():
        raise AniInputError(f"ANI metadata TSV not found: {ani_metadata}")

    with ani_metadata.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AniInputError(f"ANI metadata TSV is empty or missing a header: {ani_metadata}")

        header_map: dict[str, str] = {}
        for original in reader.fieldnames:
            normalized = normalize_header(original)
            if normalized in header_map:
                raise AniInputError(
                    "ANI metadata TSV contains duplicate normalized column names: "
                    f"{original!r} and {header_map[normalized]!r}"
                )
            header_map[normalized] = original

        normalized_busco_column = resolve_busco_column(header_map, busco_column)
        normalized_matrix_name_column = normalize_header(matrix_name_column)
        missing_columns = sorted(
            (
                REQUIRED_SCORING_METADATA_COLS
                | {normalized_busco_column, normalized_matrix_name_column}
            )
            - set(header_map)
        )
        if missing_columns:
            raise AniInputError(
                "ANI metadata TSV is missing required column(s): "
                f"{missing_columns}. Requested matrix-name column: {matrix_name_column!r}."
            )

        metadata: dict[str, Genome] = {}
        excluded: dict[str, list[str]] = {}
        seen_accessions: set[str] = set()
        for row_number, raw_row in enumerate(reader, start=2):
            row = {
                normalize_header(key): str(value or "").strip()
                for key, value in raw_row.items()
                if key is not None
            }
            accession = row.get("accession", "").strip()
            if is_missing_value(accession):
                raise AniInputError(f"ANI metadata TSV has an empty accession at line {row_number}.")
            matrix_name = row.get(normalized_matrix_name_column, "").strip()
            if is_missing_value(matrix_name):
                raise AniInputError(
                    "ANI metadata TSV has an empty matrix-name value at line "
                    f"{row_number}. Requested matrix column: {matrix_name_column!r}."
                )
            if matrix_name in metadata or matrix_name in excluded:
                raise AniInputError(
                    f"Duplicate matrix-name value {matrix_name!r} in ANI metadata TSV."
                )
            if accession in seen_accessions:
                raise AniInputError(f"Duplicate accession {accession!r} in ANI metadata TSV.")
            seen_accessions.add(accession)

            genome, reasons = build_genome_from_row(
                row,
                acc=accession,
                busco_column=normalized_busco_column,
                require_existing_path=require_existing_paths,
            )
            if genome is None:
                excluded[matrix_name] = reasons
                logging.warning(
                    "Excluding accession %r from ANI representative selection due to: %s",
                    accession,
                    ";".join(reasons),
                )
                continue
            metadata[matrix_name] = genome

    matrix_name_set = set(matrix_names)
    metadata_names = set(metadata)
    excluded_names = set(excluded)

    missing_rows = sorted(matrix_name_set - metadata_names - excluded_names)
    if missing_rows:
        raise AniInputError(
            "ANI matrix contains matrix names missing from ANI metadata TSV: "
            f"{missing_rows[:20]}"
        )

    extra_rows = sorted((metadata_names | excluded_names) - matrix_name_set)
    if extra_rows:
        logging.warning(
            "Ignoring %d ANI metadata row(s) not present in the ANI matrix (first 20): %s",
            len(extra_rows),
            extra_rows[:20],
        )

    eligible_names = [name for name in matrix_names if name in metadata]
    if not eligible_names:
        raise AniInputError("No ANI-eligible samples remain after loading ANI metadata.")
    return metadata, eligible_names


def load_cluster_metadata(
    ani_metadata: Path,
    matrix_names: Sequence[str],
    *,
    matrix_name_column: str,
) -> tuple[dict[str, str], list[str]]:
    """Load only accession-to-matrix-name metadata needed for ANI clustering."""
    if not ani_metadata.is_file():
        raise AniInputError(f"ANI metadata TSV not found: {ani_metadata}")

    with ani_metadata.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AniInputError(f"ANI metadata TSV is empty or missing a header: {ani_metadata}")

        header_map: dict[str, str] = {}
        for original in reader.fieldnames:
            normalized = normalize_header(original)
            if normalized in header_map:
                raise AniInputError(
                    "ANI metadata TSV contains duplicate normalized column names: "
                    f"{original!r} and {header_map[normalized]!r}"
                )
            header_map[normalized] = original
        normalized_matrix_name_column = normalize_header(matrix_name_column)
        missing_columns = sorted(
            {"accession", normalized_matrix_name_column} - set(header_map)
        )
        if missing_columns:
            raise AniInputError(
                "ANI metadata TSV is missing required cluster column(s): "
                f"{missing_columns}."
            )

        accession_by_matrix_name: dict[str, str] = {}
        seen_accessions: set[str] = set()
        for row_number, raw_row in enumerate(reader, start=2):
            row = {
                normalize_header(key): str(value or "").strip()
                for key, value in raw_row.items()
                if key is not None
            }
            accession = row.get("accession", "").strip()
            matrix_name = row.get(normalized_matrix_name_column, "").strip()
            if is_missing_value(accession):
                raise AniInputError(f"ANI metadata TSV has an empty accession at line {row_number}.")
            if is_missing_value(matrix_name):
                raise AniInputError(
                    "ANI metadata TSV has an empty matrix-name value at line "
                    f"{row_number}. Requested matrix column: {matrix_name_column!r}."
                )
            if matrix_name in accession_by_matrix_name:
                raise AniInputError(
                    f"Duplicate matrix-name value {matrix_name!r} in ANI metadata TSV."
                )
            if accession in seen_accessions:
                raise AniInputError(f"Duplicate accession {accession!r} in ANI metadata TSV.")
            seen_accessions.add(accession)
            accession_by_matrix_name[matrix_name] = accession

    matrix_name_set = set(matrix_names)
    missing_rows = sorted(matrix_name_set - set(accession_by_matrix_name))
    if missing_rows:
        raise AniInputError(
            "ANI matrix contains matrix names missing from ANI metadata TSV: "
            f"{missing_rows[:20]}"
        )

    extra_rows = sorted(set(accession_by_matrix_name) - matrix_name_set)
    if extra_rows:
        logging.warning(
            "Ignoring %d ANI metadata row(s) not present in the ANI matrix (first 20): %s",
            len(extra_rows),
            extra_rows[:20],
        )

    eligible_names = [name for name in matrix_names if name in accession_by_matrix_name]
    if not eligible_names:
        raise AniInputError("No ANI samples remain after loading clustering metadata.")
    return accession_by_matrix_name, eligible_names


def subset_matrix(
    names: list[str],
    ani: "np.ndarray",
    eligible_names: Sequence[str],
    name_to_idx: dict[str, int],
) -> tuple[list[str], "np.ndarray", dict[str, int]]:
    """Restrict the ANI matrix to the ordered subset of eligible matrix names."""
    import numpy as np

    eligible_name_list = list(eligible_names)
    if len(eligible_name_list) == len(names):
        return names, ani, name_to_idx

    indices = [name_to_idx[name] for name in eligible_name_list]
    subset = ani[np.ix_(indices, indices)]
    subset_name_to_idx = {name: index for index, name in enumerate(eligible_name_list)}
    logging.info(
        "Proceeding with %d ANI sample(s) after metadata filtering.",
        len(eligible_name_list),
    )
    return eligible_name_list, subset, subset_name_to_idx
