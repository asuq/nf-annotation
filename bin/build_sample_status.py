#!/usr/bin/env python3
"""Build the authoritative final sample-status table."""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from ani_common import ANI_16S_POLICIES, derive_sixteen_s_ani_exclusion_reason
from atypical_warnings import classify_atypical_warnings
import build_master_table as table_helpers
import master_table_contract
import validate_inputs


LOGGER = logging.getLogger(__name__)
PROKKA_MANIFEST_COLUMNS = ("exit_code", "gff_size", "faa_size")
CODETTA_SUMMARY_COLUMNS = (
    "Codetta_Genetic_Code",
    "Codetta_NCBI_Table_Candidates",
    "codetta_status",
    "warnings",
)
STATUS_COLUMNS_WITH_NA_DEFAULT = {"gcode", "low_quality"}


@dataclass(frozen=True)
class EffectiveAssemblyMetrics:
    """Store effective assembly metrics for ANI eligibility."""

    n50: str
    scaffolds: str
    genome_size: str


class SampleStatusError(RuntimeError):
    """Raise when the final sample-status table cannot be built safely."""


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Build the authoritative final sample_status.tsv audit table."
    )
    parser.add_argument(
        "--validated-samples",
        required=True,
        type=Path,
        help="Path to validated_samples.tsv.",
    )
    parser.add_argument(
        "--initial-status",
        required=True,
        type=Path,
        help="Path to the validation-time sample_status.tsv seed table.",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        type=Path,
        help="Path to the metadata table (CSV or TSV).",
    )
    parser.add_argument(
        "--taxonomy",
        type=Path,
        help="Optional taxonomy-expansion TSV keyed by Tax_ID.",
    )
    parser.add_argument(
        "--checkm2",
        type=Path,
        help="Optional per-sample CheckM2 summary TSV.",
    )
    parser.add_argument(
        "--16s-status",
        dest="sixteen_s_status",
        type=Path,
        help="Optional per-sample 16S status TSV.",
    )
    parser.add_argument(
        "--busco",
        action="append",
        default=[],
        type=Path,
        help="Optional per-lineage BUSCO summary TSV. May be supplied multiple times.",
    )
    parser.add_argument(
        "--codetta-summary",
        dest="codetta_summary",
        type=Path,
        help="Optional per-sample Codetta summary TSV.",
    )
    parser.add_argument(
        "--ccfinder-strains",
        type=Path,
        help="Optional CRISPRCasFinder strain summary TSV.",
    )
    parser.add_argument(
        "--prokka-manifest",
        type=Path,
        help="Optional Prokka status manifest TSV.",
    )
    parser.add_argument(
        "--ani",
        type=Path,
        help="Optional ANI summary TSV from select_ani_representatives.py.",
    )
    parser.add_argument(
        "--assembly-stats",
        type=Path,
        help="Optional in-house assembly stats TSV keyed by accession.",
    )
    parser.add_argument(
        "--primary-busco-column",
        type=str,
        help="BUSCO_<lineage> column used for ANI eligibility decisions.",
    )
    parser.add_argument(
        "--ani-16s-policy",
        choices=ANI_16S_POLICIES,
        default="complete",
        help="Require complete 16S, allow partial 16S, or ignore the 16S eligibility gate.",
    )
    parser.add_argument(
        "--busco-lineage",
        action="append",
        help="Configured BUSCO lineage name. May be supplied multiple times.",
    )
    parser.add_argument(
        "--columns",
        type=Path,
        help="Optional path to the ordered sample-status column asset override.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Path to the final sample_status.tsv output.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def split_notes(value: str | None) -> list[str]:
    """Split a semicolon-delimited notes field into stripped note fragments."""
    if value is None:
        return []
    return [note.strip() for note in value.split(";") if note.strip()]


def write_tsv(path: Path, header: Sequence[str], rows: Sequence[dict[str, str]]) -> None:
    """Write a TSV with a fixed header."""
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


def resolve_output_columns(
    path: Path | None = None,
    busco_lineages: Sequence[str] | None = None,
) -> tuple[list[str], tuple[str, ...]]:
    """Resolve the sample-status contract from an override, lineages, or defaults."""
    try:
        return master_table_contract.resolve_sample_status_columns(
            path=path,
            busco_lineages=busco_lineages,
        )
    except ValueError as error:
        raise SampleStatusError(str(error)) from error


def load_initial_status(
    path: Path,
    *,
    output_columns: Sequence[str],
    validated_samples: Sequence[dict[str, str]],
) -> dict[str, dict[str, str]]:
    """Load the seed sample-status table and validate its identity columns."""
    try:
        header, rows = table_helpers.read_table(path, delimiter="\t")
        table_helpers.require_columns(header, output_columns, "initial sample-status table")
        accession_column = table_helpers.detect_key_column(header, ("accession",))
        status_index = table_helpers.build_index(rows, accession_column, "initial sample-status table")
    except table_helpers.MasterTableError as error:
        raise SampleStatusError(str(error)) from error

    expected_accessions = {row["accession"] for row in validated_samples}
    actual_accessions = set(status_index)
    if actual_accessions != expected_accessions:
        missing = sorted(expected_accessions - actual_accessions)
        extra = sorted(actual_accessions - expected_accessions)
        details: list[str] = []
        if missing:
            details.append("missing: " + ", ".join(missing))
        if extra:
            details.append("extra: " + ", ".join(extra))
        raise SampleStatusError(
            "initial sample-status table does not match validated_samples.tsv: "
            + "; ".join(details)
        )

    for sample_row in validated_samples:
        accession = sample_row["accession"]
        status_row = status_index[accession]
        if status_row.get("internal_id", "") != sample_row["internal_id"]:
            raise SampleStatusError(
                f"initial sample-status table disagrees on internal_id for accession {accession!r}."
            )
        status_is_new = normalise_is_new_value(
            status_row.get("is_new", ""),
            accession=accession,
            table_name="initial sample-status table",
        )
        sample_is_new = normalise_is_new_value(
            sample_row["is_new"],
            accession=accession,
            table_name="validated_samples.tsv",
        )
        if status_is_new != sample_is_new:
            raise SampleStatusError(
                f"initial sample-status table disagrees on is_new for accession {accession!r}."
            )
    return status_index


def normalise_is_new_value(
    value: str,
    *,
    accession: str,
    table_name: str,
) -> str:
    """Normalise one `is_new` value or raise a table-specific status error."""
    try:
        return validate_inputs.normalise_boolean(value)
    except validate_inputs.ValidationError as error:
        raise SampleStatusError(
            f"{table_name} contains invalid is_new value {value!r} for accession {accession!r}."
        ) from error


def load_annotation_manifest(
    path: Path | None,
    required_columns: Sequence[str],
    table_name: str,
    validated_accessions: set[str],
) -> dict[str, dict[str, str]]:
    """Load one optional annotation manifest keyed by accession."""
    try:
        return table_helpers.load_optional_accession_index(
            path,
            required_columns,
            table_name,
            validated_accessions,
        )
    except table_helpers.MasterTableError as error:
        raise SampleStatusError(str(error)) from error


def parse_nonnegative_int(
    value: str | None,
    *,
    table_name: str,
    accession: str,
    column_name: str,
) -> int:
    """Parse a non-negative integer manifest field."""
    if value is None or not value.strip():
        raise SampleStatusError(
            f"{table_name} is missing {column_name} for accession {accession!r}."
        )
    try:
        parsed = int(value)
    except ValueError as error:
        raise SampleStatusError(
            f"{table_name} has a non-integer {column_name} for accession {accession!r}: {value!r}."
        ) from error
    if parsed < 0:
        raise SampleStatusError(
            f"{table_name} has a negative {column_name} for accession {accession!r}: {value!r}."
        )
    return parsed


def has_non_na_values(row: dict[str, str], columns: Sequence[str]) -> bool:
    """Return True when any requested column contains a non-missing value."""
    return any(not table_helpers.is_missing(row.get(column)) for column in columns)


def derive_taxonomy_status(
    metadata_row: dict[str, str],
    taxonomy_index: dict[str, dict[str, str]],
    taxonomy_requested: bool,
) -> tuple[str, list[str]]:
    """Return taxonomy status and warning tokens for one sample."""
    if not taxonomy_requested:
        return "na", []
    tax_id_column = table_helpers.find_column_by_normalised_name(tuple(metadata_row), "tax_id")
    if tax_id_column is None:
        return "na", []
    tax_id = metadata_row.get(tax_id_column, "NA")
    if table_helpers.is_missing(tax_id):
        return "na", []
    if tax_id in taxonomy_index:
        return "done", []
    return "failed", ["taxonomy_missing"]


def derive_barrnap_status(
    accession: str,
    sixteen_s_index: dict[str, dict[str, str]],
    sixteen_s_requested: bool,
) -> tuple[str, list[str]]:
    """Return Barrnap status and warning tokens for one sample."""
    if not sixteen_s_requested:
        return "na", []
    row = sixteen_s_index.get(accession)
    if row is None:
        return "failed", ["missing_16s_summary"]
    warnings = table_helpers.split_tokens(row.get("warnings", ""))
    status_value = row.get("16S", "NA") or "NA"
    if status_value == "NA":
        return "failed", warnings or ["invalid_barrnap_output"]
    return "done", warnings


def derive_checkm2_statuses(
    accession: str,
    checkm2_index: dict[str, dict[str, str]],
    checkm2_requested: bool,
) -> tuple[str, str, str, list[str], str, str]:
    """Return CheckM2-related statuses, warnings, and chosen QC values."""
    if not checkm2_requested:
        return "na", "na", "na", [], "NA", "NA"
    row = checkm2_index.get(accession)
    if row is None:
        return "failed", "failed", "failed", ["missing_checkm2_summary"], "NA", "NA"

    warnings = table_helpers.split_tokens(row.get("warnings", ""))
    gcode4_columns = [column for column in table_helpers.CHECKM2_COLUMNS if column.endswith("_gcode4")]
    gcode11_columns = [column for column in table_helpers.CHECKM2_COLUMNS if column.endswith("_gcode11")]
    has_gcode4 = has_non_na_values(row, gcode4_columns)
    has_gcode11 = has_non_na_values(row, gcode11_columns)

    gcode4_status = "failed" if "checkm2_gcode4_failed" in warnings or not has_gcode4 else "done"
    gcode11_status = "failed" if "checkm2_gcode11_failed" in warnings or not has_gcode11 else "done"

    gcode_value = row.get("Gcode", "NA") or "NA"
    low_quality_value = row.get("Low_quality", "NA") or "NA"
    gcode_status = "done" if gcode_value in {"4", "11"} else "failed"
    return gcode4_status, gcode11_status, gcode_status, warnings, gcode_value, low_quality_value


def derive_busco_statuses(
    accession: str,
    busco_index: dict[str, dict[str, table_helpers.BuscoLineageSummary]],
    provided_busco_columns: set[str],
    output_columns: Sequence[str],
) -> tuple[dict[str, str], list[str], dict[str, str]]:
    """Return BUSCO status columns, warning tokens, and compact BUSCO values."""
    status_values = {
        column: "na"
        for column in output_columns
        if column.startswith("busco_") and column.endswith("_status")
    }
    busco_values: dict[str, str] = {}
    warnings: list[str] = []
    sample_busco = busco_index.get(accession, {})

    for status_column in status_values:
        lineage = status_column.removeprefix("busco_").removesuffix("_status")
        busco_column = f"BUSCO_{lineage}"
        if busco_column not in provided_busco_columns:
            continue
        summary = sample_busco.get(busco_column)
        if summary is None:
            status_values[status_column] = "failed"
            warnings.append(f"missing_{busco_column.lower()}_summary")
            busco_values[busco_column] = "NA"
            continue
        status_values[status_column] = summary.status if summary.status else "done"
        warnings.extend(table_helpers.split_tokens(summary.warnings))
        busco_values[busco_column] = summary.value or "NA"
    return status_values, warnings, busco_values


def derive_ccfinder_status(
    accession: str,
    ccfinder_index: dict[str, dict[str, str]],
    ccfinder_requested: bool,
    gcode_value: str,
) -> tuple[str, list[str]]:
    """Return CRISPRCasFinder status and warning tokens for one sample."""
    row = ccfinder_index.get(accession)
    if row is not None:
        warnings = table_helpers.split_tokens(row.get("warnings", ""))
        status = row.get("ccfinder_status", "") or "done"
        return status, warnings
    if gcode_value == "NA":
        return ("skipped", []) if ccfinder_requested else ("na", [])
    if not ccfinder_requested:
        return "na", []
    return "failed", ["missing_ccfinder_summary"]


def derive_codetta_status(
    accession: str,
    codetta_index: dict[str, dict[str, str]],
    codetta_requested: bool,
) -> tuple[str, list[str]]:
    """Return Codetta status and warning tokens for one sample."""
    row = codetta_index.get(accession)
    if row is not None:
        warnings = table_helpers.split_tokens(row.get("warnings", ""))
        status = row.get("codetta_status", "") or "done"
        return status, warnings
    if not codetta_requested:
        return "na", []
    return "failed", ["missing_codetta_summary"]


def derive_ani_decision(
    accession: str,
    *,
    gcode_value: str,
    low_quality_value: str,
    sixteen_s_value: str,
    assembly_level_value: str,
    metadata_row: dict[str, str],
    assembly_metrics: EffectiveAssemblyMetrics,
    ani_index: dict[str, dict[str, str]],
    ani_requested: bool,
    primary_busco_value: str,
    ani_16s_policy: str = "complete",
) -> tuple[str, str]:
    """Return ANI inclusion status and exclusion reasons for one sample."""
    if not ani_requested:
        return "na", ""

    exclusion_reasons: list[str] = []
    atypical_value = table_helpers.detect_metadata_value(
        metadata_row, "Atypical_Warnings"
    )
    is_atypical, is_exception = classify_atypical_warnings(atypical_value)

    if gcode_value == "NA":
        exclusion_reasons.append("gcode_na")
    if low_quality_value == "true":
        exclusion_reasons.append("low_quality")
    elif low_quality_value == "NA" and gcode_value in {"4", "11"}:
        exclusion_reasons.append("missing_low_quality")

    sixteen_s_reason = derive_sixteen_s_ani_exclusion_reason(
        sixteen_s_value,
        policy=ani_16s_policy,
    )
    if sixteen_s_reason is not None:
        exclusion_reasons.append(sixteen_s_reason)

    if is_atypical and not is_exception:
        exclusion_reasons.append("atypical")

    if table_helpers.is_missing(primary_busco_value):
        exclusion_reasons.append("missing_primary_busco")
    if table_helpers.is_missing(assembly_level_value):
        exclusion_reasons.append("missing_assembly_level")
    if table_helpers.is_missing(assembly_metrics.n50):
        exclusion_reasons.append("missing_n50")
    if table_helpers.is_missing(assembly_metrics.scaffolds):
        exclusion_reasons.append("missing_scaffolds")
    if table_helpers.is_missing(assembly_metrics.genome_size):
        exclusion_reasons.append("missing_genome_size")

    if accession in ani_index:
        if exclusion_reasons:
            raise SampleStatusError(
                f"ANI summary contains ineligible accession {accession!r}: "
                + table_helpers.join_tokens(exclusion_reasons)
            )
        return "true", ""

    if exclusion_reasons:
        return "false", table_helpers.join_tokens(exclusion_reasons)

    raise SampleStatusError(f"ANI summary is missing eligible accession {accession!r}.")


def derive_annotation_status(
    accession: str,
    *,
    manifest_index: dict[str, dict[str, str]],
    manifest_requested: bool,
    gcode_value: str,
    table_name: str,
    failed_warning: str,
    missing_warning: str,
    has_outputs,
    missing_outputs_warning: str | None = None,
) -> tuple[str, list[str]]:
    """Return one annotation-tool status from its per-sample manifest."""
    if not manifest_requested:
        return "na", []
    if gcode_value == "NA":
        return "skipped", []

    row = manifest_index.get(accession)
    if row is None:
        return "failed", [missing_warning]

    exit_code = parse_nonnegative_int(
        row.get("exit_code"),
        table_name=table_name,
        accession=accession,
        column_name="exit_code",
    )
    if exit_code != 0:
        return "failed", [failed_warning]

    if not has_outputs(row, accession):
        warning = missing_outputs_warning or failed_warning
        return "failed", [warning]

    return "done", []


def prokka_outputs_present(row: dict[str, str], accession: str) -> bool:
    """Return True when Prokka emitted non-empty GFF and FAA outputs."""
    gff_size = parse_nonnegative_int(
        row.get("gff_size"),
        table_name="Prokka manifest",
        accession=accession,
        column_name="gff_size",
    )
    faa_size = parse_nonnegative_int(
        row.get("faa_size"),
        table_name="Prokka manifest",
        accession=accession,
        column_name="faa_size",
    )
    return gff_size > 0 and faa_size > 0


def build_status_row(
    sample_row: dict[str, str],
    *,
    initial_row: dict[str, str],
    output_columns: Sequence[str],
    metadata_row: dict[str, str],
    assembly_stats_row: dict[str, str] | None,
    taxonomy_index: dict[str, dict[str, str]],
    taxonomy_requested: bool,
    checkm2_index: dict[str, dict[str, str]],
    checkm2_requested: bool,
    sixteen_s_index: dict[str, dict[str, str]],
    sixteen_s_requested: bool,
    busco_index: dict[str, dict[str, table_helpers.BuscoLineageSummary]],
    provided_busco_columns: set[str],
    codetta_index: dict[str, dict[str, str]],
    codetta_requested: bool,
    ccfinder_index: dict[str, dict[str, str]],
    ccfinder_requested: bool,
    prokka_index: dict[str, dict[str, str]],
    prokka_requested: bool,
    ani_index: dict[str, dict[str, str]],
    ani_requested: bool,
    primary_busco_column: str | None,
    ani_16s_policy: str,
) -> dict[str, str]:
    """Build one final status row by overlaying derived statuses on the seed row."""
    row = {column: initial_row.get(column, "") for column in output_columns}
    for column in output_columns:
        if column not in row or row[column] == "":
            if column.endswith("_status") or column == "ani_included":
                row[column] = "na"
            elif column in STATUS_COLUMNS_WITH_NA_DEFAULT:
                row[column] = "NA"
            else:
                row[column] = ""

    accession = sample_row["accession"]
    row["accession"] = accession
    row["internal_id"] = sample_row["internal_id"]
    row["is_new"] = sample_row["is_new"]

    warnings = table_helpers.split_tokens(initial_row.get("warnings", ""))
    notes = split_notes(initial_row.get("notes", ""))

    taxonomy_status, taxonomy_warnings = derive_taxonomy_status(
        metadata_row=metadata_row,
        taxonomy_index=taxonomy_index,
        taxonomy_requested=taxonomy_requested,
    )
    row["taxonomy_status"] = taxonomy_status
    warnings.extend(taxonomy_warnings)

    barrnap_status, barrnap_warnings = derive_barrnap_status(
        accession=accession,
        sixteen_s_index=sixteen_s_index,
        sixteen_s_requested=sixteen_s_requested,
    )
    row["barrnap_status"] = barrnap_status
    warnings.extend(barrnap_warnings)

    (
        row["checkm2_gcode4_status"],
        row["checkm2_gcode11_status"],
        row["gcode_status"],
        checkm2_warnings,
        gcode_value,
        low_quality_value,
    ) = derive_checkm2_statuses(
        accession=accession,
        checkm2_index=checkm2_index,
        checkm2_requested=checkm2_requested,
    )
    row["gcode"] = gcode_value
    row["low_quality"] = low_quality_value
    warnings.extend(checkm2_warnings)

    busco_status_values, busco_warnings, busco_values = derive_busco_statuses(
        accession=accession,
        busco_index=busco_index,
        provided_busco_columns=provided_busco_columns,
        output_columns=output_columns,
    )
    row.update(busco_status_values)
    warnings.extend(busco_warnings)

    codetta_status, codetta_warnings = derive_codetta_status(
        accession=accession,
        codetta_index=codetta_index,
        codetta_requested=codetta_requested,
    )
    row["codetta_status"] = codetta_status
    warnings.extend(codetta_warnings)

    ccfinder_status, ccfinder_warnings = derive_ccfinder_status(
        accession=accession,
        ccfinder_index=ccfinder_index,
        ccfinder_requested=ccfinder_requested,
        gcode_value=gcode_value,
    )
    row["ccfinder_status"] = ccfinder_status
    warnings.extend(ccfinder_warnings)

    prokka_status, prokka_warnings = derive_annotation_status(
        accession,
        manifest_index=prokka_index,
        manifest_requested=prokka_requested,
        gcode_value=gcode_value,
        table_name="Prokka manifest",
        failed_warning="prokka_failed",
        missing_warning="missing_prokka_result",
        has_outputs=prokka_outputs_present,
        missing_outputs_warning="missing_prokka_outputs",
    )
    row["prokka_status"] = prokka_status
    warnings.extend(prokka_warnings)

    sixteen_s_value = sixteen_s_index.get(accession, {}).get("16S", "NA") or "NA"
    primary_busco_value = "NA"
    if primary_busco_column is not None:
        primary_busco_value = busco_values.get(primary_busco_column, "NA")
    assembly_level_value = table_helpers.choose_assembly_level(sample_row, metadata_row)
    assembly_metrics = EffectiveAssemblyMetrics(
        n50=table_helpers.resolve_assembly_metric_value(
            metadata_row,
            assembly_stats_row,
            "N50",
        ),
        scaffolds=table_helpers.resolve_assembly_metric_value(
            metadata_row,
            assembly_stats_row,
            "Scaffolds",
        ),
        genome_size=table_helpers.resolve_assembly_metric_value(
            metadata_row,
            assembly_stats_row,
            "Genome_Size",
        ),
    )

    row["ani_included"], row["ani_exclusion_reason"] = derive_ani_decision(
        accession,
        gcode_value=gcode_value,
        low_quality_value=low_quality_value,
        sixteen_s_value=sixteen_s_value,
        assembly_level_value=assembly_level_value,
        metadata_row=metadata_row,
        assembly_metrics=assembly_metrics,
        ani_index=ani_index,
        ani_requested=ani_requested,
        primary_busco_value=primary_busco_value,
        ani_16s_policy=ani_16s_policy,
    )

    row["warnings"] = table_helpers.join_tokens(warnings)
    row["notes"] = table_helpers.join_notes(notes)
    return row


def run_build(args: argparse.Namespace) -> None:
    """Build the final sample-status table from validation and downstream summaries."""
    if args.ani is not None and not args.primary_busco_column:
        raise SampleStatusError(
            "--primary-busco-column is required when --ani is supplied."
        )

    try:
        validated_samples = table_helpers.load_validated_samples(args.validated_samples)
        validated_accessions = {row["accession"] for row in validated_samples}
        output_columns, contract_busco_lineages = resolve_output_columns(
            args.columns,
            busco_lineages=args.busco_lineage,
        )
        expected_busco_columns = {
            f"BUSCO_{lineage}" for lineage in contract_busco_lineages
        }
        initial_status_index = load_initial_status(
            args.initial_status,
            output_columns=output_columns,
            validated_samples=validated_samples,
        )
        metadata_header, metadata_key_column, metadata_index = (
            table_helpers.load_metadata(args.metadata)
        )
        assembly_stats_index = table_helpers.load_assembly_stats_index(
            args.assembly_stats,
            validated_accessions,
        )
        taxonomy_index = table_helpers.load_taxonomy_index(args.taxonomy)
        checkm2_index = table_helpers.load_optional_accession_index(
            args.checkm2,
            table_helpers.CHECKM2_COLUMNS,
            "CheckM2 summary table",
            validated_accessions,
        )
        sixteen_s_index = table_helpers.load_optional_accession_index(
            args.sixteen_s_status,
            table_helpers.SIXTEEN_S_COLUMNS,
            "16S status table",
            validated_accessions,
        )
        codetta_index = table_helpers.load_optional_accession_index(
            args.codetta_summary,
            CODETTA_SUMMARY_COLUMNS,
            "Codetta summary table",
            validated_accessions,
        )
        ccfinder_index = table_helpers.load_optional_accession_index(
            args.ccfinder_strains,
            table_helpers.CRISPR_COLUMNS,
            "CRISPR strain summary table",
            validated_accessions,
        )
        prokka_index = load_annotation_manifest(
            args.prokka_manifest,
            PROKKA_MANIFEST_COLUMNS,
            "Prokka manifest",
            validated_accessions,
        )
        ani_index = table_helpers.load_optional_accession_index(
            args.ani,
            table_helpers.ANI_COLUMNS,
            "ANI summary table",
            validated_accessions,
        )
        busco_index, provided_busco_columns = table_helpers.load_busco_index(
            args.busco,
            validated_accessions,
        )
    except table_helpers.MasterTableError as error:
        raise SampleStatusError(str(error)) from error
    except ValueError as error:
        raise SampleStatusError(str(error)) from error

    unexpected_busco_columns = sorted(provided_busco_columns - expected_busco_columns)
    if unexpected_busco_columns:
        raise SampleStatusError(
            "BUSCO summaries contain lineage columns not present in the sample-status "
            "contract: " + ", ".join(unexpected_busco_columns)
        )
    if (
        args.primary_busco_column is not None
        and args.primary_busco_column not in expected_busco_columns
    ):
        raise SampleStatusError(
            "--primary-busco-column is not present in the sample-status BUSCO contract: "
            + args.primary_busco_column
        )

    rows: list[dict[str, str]] = []
    for sample_row in validated_samples:
        accession = sample_row["accession"]
        assembly_stats_row = assembly_stats_index.get(accession)
        metadata_row = table_helpers.build_metadata_row(
            sample_row=sample_row,
            metadata_header=metadata_header,
            metadata_key_column=metadata_key_column,
            metadata_index=metadata_index,
            assembly_stats_row=assembly_stats_row,
        )
        rows.append(
            build_status_row(
                sample_row,
                initial_row=initial_status_index[accession],
                output_columns=output_columns,
                metadata_row=metadata_row,
                assembly_stats_row=assembly_stats_row,
                taxonomy_index=taxonomy_index,
                taxonomy_requested=args.taxonomy is not None,
                checkm2_index=checkm2_index,
                checkm2_requested=args.checkm2 is not None,
                sixteen_s_index=sixteen_s_index,
                sixteen_s_requested=args.sixteen_s_status is not None,
                busco_index=busco_index,
                provided_busco_columns=provided_busco_columns,
                codetta_index=codetta_index,
                codetta_requested=args.codetta_summary is not None,
                ccfinder_index=ccfinder_index,
                ccfinder_requested=args.ccfinder_strains is not None,
                prokka_index=prokka_index,
                prokka_requested=args.prokka_manifest is not None,
                ani_index=ani_index,
                ani_requested=args.ani is not None,
                primary_busco_column=args.primary_busco_column,
                ani_16s_policy=args.ani_16s_policy,
            )
        )

    write_tsv(args.output, output_columns, rows)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the sample-status writer CLI."""
    args = parse_args(argv)
    configure_logging()
    try:
        run_build(args)
    except SampleStatusError as error:
        LOGGER.error(str(error))
        return 1

    LOGGER.info("Wrote sample status to %s.", args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
