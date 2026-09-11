"""Define the locked master-table column contract for v1."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence


DEFAULT_BUSCO_LINEAGES = ("bacillota_odb12", "mycoplasmatota_odb12")
TAXONOMY_COLUMNS = (
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)
MASTER_TABLE_PREFIX_COLUMNS = ("is_new",)
CHECKM2_COLUMNS = (
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
)
ASSEMBLY_DERIVED_COLUMNS = ("GC_Content",)
CODETTA_COLUMNS = ("Codetta_Genetic_Code", "Codetta_NCBI_Table_Candidates")
GCODE_PROVENANCE_COLUMNS = (
    "Gcode_Rule",
    "Gcode_Length_Ratio",
    "Gcode_Length_Ratio_Threshold",
    "Gcode_Selection_Reason",
)
GCODE_QC_COLUMNS = (
    "Gcode", *GCODE_PROVENANCE_COLUMNS, *CODETTA_COLUMNS, "Low_quality", "16S"
)
CRISPR_COLUMNS = ("CRISPRS", "SPACERS_SUM", "CRISPR_FRAC")
ANI_COLUMNS = ("Cluster_ID", "Is_Representative", "ANI_to_Representative", "Score")
DEFAULT_APPEND_COLUMNS_ASSET = (
    Path(__file__).resolve().parents[1]
    / "assets"
    / "tables"
    / "contracts"
    / "master_table_append_columns.txt"
)
SAMPLE_STATUS_PREFIX_COLUMNS = (
    "accession",
    "internal_id",
    "is_new",
    "validation_status",
    "taxonomy_status",
    "barrnap_status",
    "checkm2_gcode4_status",
    "checkm2_gcode11_status",
    "gcode_status",
    "gcode",
    "low_quality",
)
SAMPLE_STATUS_SUFFIX_COLUMNS = (
    "codetta_status",
    "prokka_status",
    "ccfinder_status",
    "padloc_status",
    "eggnog_status",
    "ani_included",
    "ani_exclusion_reason",
    "warnings",
    "notes",
)
DEFAULT_SAMPLE_STATUS_COLUMNS_ASSET = (
    Path(__file__).resolve().parents[1]
    / "assets"
    / "tables"
    / "contracts"
    / "sample_status_columns.txt"
)


def normalise_busco_lineages(
    busco_lineages: Sequence[str] | None = None,
) -> tuple[str, ...]:
    """Validate BUSCO lineage names and preserve their order."""
    raw_lineages = busco_lineages or DEFAULT_BUSCO_LINEAGES
    normalised = tuple(lineage.strip() for lineage in raw_lineages if lineage.strip())
    if not normalised:
        raise ValueError("At least one BUSCO lineage is required.")
    if len(set(normalised)) != len(normalised):
        raise ValueError("BUSCO lineages must be unique.")
    return normalised


def build_append_columns(busco_lineages: Sequence[str] | None = None) -> list[str]:
    """Return the appended derived columns in locked output order."""
    busco_columns = [
        f"BUSCO_{lineage}" for lineage in normalise_busco_lineages(busco_lineages)
    ]
    return [
        *MASTER_TABLE_PREFIX_COLUMNS,
        *TAXONOMY_COLUMNS,
        *CHECKM2_COLUMNS,
        *ASSEMBLY_DERIVED_COLUMNS,
        *GCODE_QC_COLUMNS,
        *busco_columns,
        *CRISPR_COLUMNS,
        *ANI_COLUMNS,
    ]


def build_master_table_columns(
    metadata_columns: Sequence[str],
    busco_lineages: Sequence[str] | None = None,
) -> list[str]:
    """Return the full master-table header from metadata and derived contracts."""
    cleaned_metadata_columns = [column.strip() for column in metadata_columns]
    if not cleaned_metadata_columns:
        raise ValueError("Metadata columns cannot be empty.")
    if any(not column for column in cleaned_metadata_columns):
        raise ValueError("Metadata columns cannot contain empty names.")
    if len(set(cleaned_metadata_columns)) != len(cleaned_metadata_columns):
        raise ValueError("Metadata columns must be unique.")

    append_columns = build_append_columns(busco_lineages)
    overlap = sorted(set(cleaned_metadata_columns) & set(append_columns))
    if overlap:
        raise ValueError(
            "Metadata columns overlap with derived master-table columns: "
            + ", ".join(overlap)
        )
    return [*cleaned_metadata_columns, *append_columns]


def read_append_columns_asset(path: Path = DEFAULT_APPEND_COLUMNS_ASSET) -> list[str]:
    """Read the repository append-column asset file."""
    return [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def build_sample_status_columns(
    busco_lineages: Sequence[str] | None = None,
) -> list[str]:
    """Return the final sample-status columns in locked output order."""
    busco_columns = [
        f"busco_{lineage}_status" for lineage in normalise_busco_lineages(busco_lineages)
    ]
    return [
        *SAMPLE_STATUS_PREFIX_COLUMNS,
        *busco_columns,
        *SAMPLE_STATUS_SUFFIX_COLUMNS,
    ]


def read_sample_status_columns_asset(
    path: Path = DEFAULT_SAMPLE_STATUS_COLUMNS_ASSET,
) -> list[str]:
    """Read the repository sample-status column asset file."""
    return [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def extract_busco_lineages_from_append_columns(
    append_columns: Sequence[str],
) -> tuple[str, ...]:
    """Return BUSCO lineage names encoded in one append-column contract."""
    prefix_columns = [
        *MASTER_TABLE_PREFIX_COLUMNS,
        *TAXONOMY_COLUMNS,
        *CHECKM2_COLUMNS,
        *ASSEMBLY_DERIVED_COLUMNS,
        *GCODE_QC_COLUMNS,
    ]
    suffix_columns = [*CRISPR_COLUMNS, *ANI_COLUMNS]
    append_list = [column.strip() for column in append_columns]

    if append_list[: len(prefix_columns)] != prefix_columns:
        raise ValueError("Append-column contract has an unexpected non-BUSCO prefix.")
    if append_list[-len(suffix_columns) :] != suffix_columns:
        raise ValueError("Append-column contract has an unexpected non-BUSCO suffix.")

    busco_columns = append_list[len(prefix_columns) : len(append_list) - len(suffix_columns)]
    if not busco_columns:
        raise ValueError("Append-column contract must contain at least one BUSCO column.")
    if any(not column.startswith("BUSCO_") for column in busco_columns):
        raise ValueError("Append-column contract contains a non-BUSCO column in the BUSCO block.")

    busco_lineages = normalise_busco_lineages(
        [column.removeprefix("BUSCO_") for column in busco_columns]
    )
    expected_columns = build_append_columns(busco_lineages)
    if append_list != expected_columns:
        raise ValueError("Append-column contract does not match the BUSCO-aware contract.")
    return busco_lineages


def extract_busco_lineages_from_sample_status_columns(
    sample_status_columns: Sequence[str],
) -> tuple[str, ...]:
    """Return BUSCO lineage names encoded in one sample-status contract."""
    status_list = [column.strip() for column in sample_status_columns]

    if status_list[: len(SAMPLE_STATUS_PREFIX_COLUMNS)] != list(SAMPLE_STATUS_PREFIX_COLUMNS):
        raise ValueError("Sample-status contract has an unexpected non-BUSCO prefix.")
    if status_list[-len(SAMPLE_STATUS_SUFFIX_COLUMNS) :] != list(SAMPLE_STATUS_SUFFIX_COLUMNS):
        raise ValueError("Sample-status contract has an unexpected non-BUSCO suffix.")

    busco_columns = status_list[
        len(SAMPLE_STATUS_PREFIX_COLUMNS) : len(status_list) - len(SAMPLE_STATUS_SUFFIX_COLUMNS)
    ]
    if not busco_columns:
        raise ValueError("Sample-status contract must contain at least one BUSCO status column.")
    if any(
        not column.startswith("busco_") or not column.endswith("_status")
        for column in busco_columns
    ):
        raise ValueError(
            "Sample-status contract contains a non-BUSCO status column in the BUSCO block."
        )

    busco_lineages = normalise_busco_lineages(
        [
            column.removeprefix("busco_").removesuffix("_status")
            for column in busco_columns
        ]
    )
    expected_columns = build_sample_status_columns(busco_lineages)
    if status_list != expected_columns:
        raise ValueError("Sample-status contract does not match the BUSCO-aware contract.")
    return busco_lineages


def resolve_append_columns(
    path: Path | None = None,
    busco_lineages: Sequence[str] | None = None,
) -> tuple[list[str], tuple[str, ...]]:
    """Resolve one append-column contract from an override, lineages, or defaults."""
    if path is not None:
        if not path.is_file():
            raise ValueError(f"Missing append-column contract: {path}")
        append_columns = read_append_columns_asset(path)
        return append_columns, extract_busco_lineages_from_append_columns(append_columns)

    if busco_lineages:
        resolved_lineages = normalise_busco_lineages(busco_lineages)
        return build_append_columns(resolved_lineages), resolved_lineages

    validate_default_append_columns_asset()
    return build_append_columns(DEFAULT_BUSCO_LINEAGES), DEFAULT_BUSCO_LINEAGES


def resolve_sample_status_columns(
    path: Path | None = None,
    busco_lineages: Sequence[str] | None = None,
) -> tuple[list[str], tuple[str, ...]]:
    """Resolve one sample-status contract from an override, lineages, or defaults."""
    if path is not None:
        if not path.is_file():
            raise ValueError(f"Missing sample-status column asset: {path}")
        sample_status_columns = read_sample_status_columns_asset(path)
        return (
            sample_status_columns,
            extract_busco_lineages_from_sample_status_columns(sample_status_columns),
        )

    if busco_lineages:
        resolved_lineages = normalise_busco_lineages(busco_lineages)
        return build_sample_status_columns(resolved_lineages), resolved_lineages

    validate_default_sample_status_columns_asset()
    return build_sample_status_columns(DEFAULT_BUSCO_LINEAGES), DEFAULT_BUSCO_LINEAGES


def validate_default_append_columns_asset(
    path: Path = DEFAULT_APPEND_COLUMNS_ASSET,
) -> None:
    """Raise when the default asset file drifts from the code-defined contract."""
    asset_columns = read_append_columns_asset(path)
    expected_columns = build_append_columns(DEFAULT_BUSCO_LINEAGES)
    if asset_columns != expected_columns:
        raise ValueError(
            f"Append-column asset {path} does not match the default v1 contract."
        )


def validate_default_sample_status_columns_asset(
    path: Path = DEFAULT_SAMPLE_STATUS_COLUMNS_ASSET,
) -> None:
    """Raise when the default sample-status asset drifts from the code contract."""
    asset_columns = read_sample_status_columns_asset(path)
    expected_columns = build_sample_status_columns(DEFAULT_BUSCO_LINEAGES)
    if asset_columns != expected_columns:
        raise ValueError(
            f"Sample-status column asset {path} does not match the default v1 contract."
        )
