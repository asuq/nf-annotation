#!/usr/bin/env python3
"""Run layered acceptance tests for local and SLURM executions."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import logging
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence
from urllib.request import urlopen

import master_table_contract
import validate_inputs


LOGGER = logging.getLogger(__name__)
ROOT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_WORK_ROOT = ROOT_DIR / "assets" / "fixtures" / "local" / "acceptance"
DEFAULT_SOURCE_CATALOG = ROOT_DIR / "assets" / "tables" / "acceptance" / "source_catalog.tsv"
DEFAULT_COHORT_PLAN = ROOT_DIR / "assets" / "tables" / "acceptance" / "cohort_plan.tsv"
DEFAULT_MEDIUM_SOURCE_CATALOG = (
    ROOT_DIR / "assets" / "tables" / "medium" / "source_catalog.tsv"
)
DEFAULT_MEDIUM_COHORT_PLAN = ROOT_DIR / "assets" / "tables" / "medium" / "cohort_plan.tsv"
DEFAULT_LOCAL_PROFILE = "debug,local,docker"
DEFAULT_SLURM_PROFILE = "debug,slurm,singularity"
DEFAULT_DBPREP_PROFILE = "slurm,singularity"
REAL_RUN_NOTE = (
    "CRISPRCasFinder uses params.ccfinder_container from pipeline config; "
    "the acceptance harness does not override it."
)
DBPREP_RUN_NOTE = (
    "This mode runs prepare_databases.nf on SLURM, downloads the runtime "
    "databases into explicit tool directories, and validates the prepared tree."
)
DOWNLOAD_COLUMNS = ("source_accession", "sha256", "fasta_path")
SOURCE_STATS_COLUMNS = (
    "source_accession",
    "organism_name",
    "tax_id",
    "assembly_level",
    "scaffolds",
    "n50",
    "genome_size",
    "fasta_path",
    "sha256",
)
METADATA_COLUMNS = (
    "Accession",
    "Tax_ID",
    "Organism_Name",
    "Assembly_Level",
    "N50",
    "Scaffolds",
    "Genome_Size",
    "Atypical_Warnings",
)
SAMPLE_COLUMNS = ("accession", "is_new", "assembly_level", "genome_fasta")
FINAL_TABLES = ("master_table.tsv", "sample_status.tsv", "tool_and_db_versions.tsv")
DBPREP_OUTPUTS = (
    "runtime_database_report.tsv",
    "nextflow_args.txt",
    "pipeline_info/trace.tsv",
    "pipeline_info/report.html",
    "pipeline_info/timeline.html",
    "pipeline_info/dag.html",
)
DEFAULT_DBPREP_BUSCO_LINEAGES = ("bacillota_odb12", "mycoplasmatota_odb12")
SAMPLE_STATUS_COLUMNS_ASSET = (
    ROOT_DIR / "assets" / "tables" / "contracts" / "sample_status_columns.txt"
)
REQUIRED_ROLE_TAGS = {
    "gcode4_candidate",
    "gcode11_candidate",
    "crispr_positive_candidate",
    "crispr_negative_candidate",
    "eggnog_smoke_candidate",
    "ani_cluster_candidate",
    "atypical_excluded_candidate",
    "atypical_exception_candidate",
    "atypical_mixed_candidate",
    "missing_metadata_case",
    "collision_candidate",
}
TRUE_TOKENS = {"true", "t", "yes", "y", "1"}
FALSE_TOKENS = {"false", "f", "no", "n", "0"}


class AcceptanceTestError(RuntimeError):
    """Raised when acceptance-test preparation or assertions fail."""


@dataclass(frozen=True)
class SourceRecord:
    """Describe one tracked source genome download."""

    source_accession: str
    organism_name: str
    tax_id: str
    assembly_level: str
    source_url: str


@dataclass(frozen=True)
class CohortRecord:
    """Describe one generated cohort sample."""

    accession: str
    source_accession: str
    is_new: str
    include_metadata: bool
    atypical_warnings: str
    role_tags: tuple[str, ...]


@dataclass(frozen=True)
class SourceStats:
    """Store derived FASTA statistics for one downloaded source genome."""

    source_accession: str
    organism_name: str
    tax_id: str
    assembly_level: str
    scaffolds: int
    n50: int
    genome_size: int
    fasta_path: Path
    sha256: str


@dataclass(frozen=True)
class PreparedCohort:
    """Describe the generated acceptance cohort files and metadata."""

    work_root: Path
    sample_csv: Path
    metadata_tsv: Path
    source_stats_tsv: Path
    checksums_tsv: Path
    cohort_plan: tuple[CohortRecord, ...]
    source_stats: dict[str, SourceStats]


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def normalise_boolean(value: str, *, field_name: str) -> bool:
    """Normalise one boolean-like value."""
    token = value.strip().lower()
    if token in TRUE_TOKENS:
        return True
    if token in FALSE_TOKENS:
        return False
    raise AcceptanceTestError(f"Invalid boolean value for {field_name}: {value!r}")


def normalise_optional(value: str | None) -> str:
    """Normalise blank optional values to `NA`."""
    if value is None:
        return "NA"
    cleaned = value.strip()
    return cleaned if cleaned else "NA"


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read one TSV into a header and row dictionaries."""
    if not path.is_file():
        raise AcceptanceTestError(f"Missing TSV file: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise AcceptanceTestError(f"TSV file is empty: {path}")
        return reader.fieldnames, list(reader)


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read one CSV into a header and row dictionaries."""
    if not path.is_file():
        raise AcceptanceTestError(f"Missing CSV file: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise AcceptanceTestError(f"CSV file is empty: {path}")
        return reader.fieldnames, list(reader)


def write_tsv(path: Path, header: Sequence[str], rows: Sequence[dict[str, str]]) -> None:
    """Write rows to a TSV file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(header), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_csv(path: Path, header: Sequence[str], rows: Sequence[dict[str, str]]) -> None:
    """Write rows to a CSV file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(header))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def load_source_catalog(path: Path) -> dict[str, SourceRecord]:
    """Load the tracked source catalog."""
    header, rows = read_tsv(path)
    required = ("source_accession", "organism_name", "tax_id", "assembly_level", "source_url")
    missing = [column for column in required if column not in header]
    if missing:
        raise AcceptanceTestError(
            f"Source catalog is missing required columns: {', '.join(missing)}"
        )

    catalog: dict[str, SourceRecord] = {}
    for row in rows:
        source_accession = row["source_accession"].strip()
        if not source_accession:
            raise AcceptanceTestError("Source catalog contains an empty source_accession.")
        if source_accession in catalog:
            raise AcceptanceTestError(
                f"Source catalog contains duplicate source_accession {source_accession!r}."
            )
        catalog[source_accession] = SourceRecord(
            source_accession=source_accession,
            organism_name=row["organism_name"].strip(),
            tax_id=row["tax_id"].strip(),
            assembly_level=row["assembly_level"].strip(),
            source_url=row["source_url"].strip(),
        )
    return catalog


def parse_role_tags(token: str) -> tuple[str, ...]:
    """Parse a semicolon-delimited role-tag string."""
    role_tags = tuple(
        dict.fromkeys(tag.strip() for tag in token.split(";") if tag.strip())
    )
    if not role_tags:
        raise AcceptanceTestError("Cohort plan row is missing role_tags.")
    return role_tags


def load_cohort_plan(path: Path, catalog: dict[str, SourceRecord]) -> tuple[CohortRecord, ...]:
    """Load the tracked cohort plan."""
    header, rows = read_tsv(path)
    required = (
        "accession",
        "source_accession",
        "is_new",
        "include_metadata",
        "atypical_warnings",
        "role_tags",
    )
    missing = [column for column in required if column not in header]
    if missing:
        raise AcceptanceTestError(
            f"Cohort plan is missing required columns: {', '.join(missing)}"
        )

    plan: list[CohortRecord] = []
    seen_accessions: set[str] = set()
    seen_required_roles: set[str] = set()
    for row in rows:
        accession = row["accession"].strip()
        if not accession:
            raise AcceptanceTestError("Cohort plan contains an empty accession.")
        if accession in seen_accessions:
            raise AcceptanceTestError(
                f"Cohort plan contains duplicate accession {accession!r}."
            )
        seen_accessions.add(accession)
        validate_inputs.ensure_path_safe_accession(accession)

        source_accession = row["source_accession"].strip()
        if source_accession not in catalog:
            raise AcceptanceTestError(
                f"Cohort plan accession {accession!r} references unknown source_accession "
                f"{source_accession!r}."
            )

        is_new = "true" if normalise_boolean(row["is_new"], field_name="is_new") else "false"
        include_metadata = normalise_boolean(
            row["include_metadata"],
            field_name="include_metadata",
        )
        if is_new == "false" and not include_metadata:
            raise AcceptanceTestError(
                f"Cohort plan accession {accession!r} cannot be is_new=false without metadata."
            )

        role_tags = parse_role_tags(row["role_tags"])
        seen_required_roles.update(set(role_tags) & REQUIRED_ROLE_TAGS)
        plan.append(
            CohortRecord(
                accession=accession,
                source_accession=source_accession,
                is_new=is_new,
                include_metadata=include_metadata,
                atypical_warnings=normalise_optional(row["atypical_warnings"]),
                role_tags=role_tags,
            )
        )

    missing_roles = sorted(REQUIRED_ROLE_TAGS - seen_required_roles)
    if missing_roles:
        raise AcceptanceTestError(
            "Cohort plan is missing required role tags: " + ", ".join(missing_roles)
        )
    eggnog_smoke_candidates = [record for record in plan if "eggnog_smoke_candidate" in record.role_tags]
    if len(eggnog_smoke_candidates) != 1:
        raise AcceptanceTestError("Cohort plan must contain exactly one eggnog_smoke_candidate.")
    return tuple(plan)


def stream_copy(source, handle) -> None:
    """Copy one binary stream to a writable handle."""
    while True:
        chunk = source.read(1024 * 1024)
        if not chunk:
            return
        handle.write(chunk)


def download_source_fasta(
    record: SourceRecord,
    *,
    downloads_dir: Path,
    offline: bool,
) -> Path:
    """Download one source FASTA if it is not already cached."""
    destination = downloads_dir / f"{record.source_accession}_genomic.fna"
    if destination.is_file():
        return destination
    if offline:
        raise AcceptanceTestError(
            f"Missing cached FASTA for {record.source_accession} while offline: {destination}"
        )

    downloads_dir.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(".fna.tmp")
    LOGGER.info("Downloading %s to %s", record.source_accession, destination)
    try:
        with urlopen(record.source_url, timeout=300) as response:
            with temporary.open("wb") as handle:
                if record.source_url.endswith(".gz"):
                    with gzip.GzipFile(fileobj=response) as gz_handle:
                        stream_copy(gz_handle, handle)
                else:
                    stream_copy(response, handle)
        temporary.replace(destination)
    except Exception as error:
        temporary.unlink(missing_ok=True)
        raise AcceptanceTestError(
            f"Failed to download {record.source_accession} from {record.source_url}: {error}"
        ) from error
    return destination


def calculate_fasta_stats(path: Path) -> tuple[int, int, int]:
    """Return scaffolds, N50, and genome size for one FASTA file."""
    if not path.is_file():
        raise AcceptanceTestError(f"Missing FASTA file: {path}")

    sequence_lengths: list[int] = []
    current_length = 0
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_length:
                    sequence_lengths.append(current_length)
                current_length = 0
                continue
            current_length += len(line)
    if current_length:
        sequence_lengths.append(current_length)

    if not sequence_lengths:
        raise AcceptanceTestError(f"FASTA file contains no sequences: {path}")

    genome_size = sum(sequence_lengths)
    sorted_lengths = sorted(sequence_lengths, reverse=True)
    half_total = genome_size / 2
    running_total = 0
    n50 = 0
    for length in sorted_lengths:
        running_total += length
        if running_total >= half_total:
            n50 = length
            break
    return len(sequence_lengths), n50, genome_size


def sha256_file(path: Path) -> str:
    """Calculate the SHA-256 digest for one file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                return digest.hexdigest()
            digest.update(chunk)


def load_existing_checksums(path: Path) -> dict[str, str]:
    """Load any previously recorded download checksums."""
    if not path.is_file():
        return {}
    _header, rows = read_tsv(path)
    checksums: dict[str, str] = {}
    for row in rows:
        checksums[row["source_accession"]] = row["sha256"]
    return checksums


def record_and_validate_checksums(
    source_stats: Sequence[SourceStats],
    path: Path,
) -> None:
    """Persist checksums and fail if any cached download has drifted."""
    existing = load_existing_checksums(path)
    rows: list[dict[str, str]] = []
    for stats in sorted(source_stats, key=lambda item: item.source_accession):
        recorded = existing.get(stats.source_accession)
        if recorded is not None and recorded != stats.sha256:
            raise AcceptanceTestError(
                f"download_checksum_drift:{stats.source_accession}"
            )
        rows.append(
            {
                "source_accession": stats.source_accession,
                "sha256": stats.sha256,
                "fasta_path": str(stats.fasta_path.resolve()),
            }
        )
    write_tsv(path, DOWNLOAD_COLUMNS, rows)


def build_source_stats(
    catalog: dict[str, SourceRecord],
    plan: Sequence[CohortRecord],
    *,
    downloads_dir: Path,
    offline: bool,
) -> dict[str, SourceStats]:
    """Download and summarise all source genomes referenced by the cohort plan."""
    source_stats: dict[str, SourceStats] = {}
    for source_accession in sorted({record.source_accession for record in plan}):
        source = catalog[source_accession]
        fasta_path = download_source_fasta(
            source,
            downloads_dir=downloads_dir,
            offline=offline,
        )
        scaffolds, n50, genome_size = calculate_fasta_stats(fasta_path)
        source_stats[source_accession] = SourceStats(
            source_accession=source_accession,
            organism_name=source.organism_name,
            tax_id=source.tax_id,
            assembly_level=source.assembly_level,
            scaffolds=scaffolds,
            n50=n50,
            genome_size=genome_size,
            fasta_path=fasta_path.resolve(),
            sha256=sha256_file(fasta_path),
        )
    return source_stats


def build_sample_rows(
    plan: Sequence[CohortRecord],
    source_stats: dict[str, SourceStats],
) -> list[dict[str, str]]:
    """Build the generated sample manifest rows."""
    rows: list[dict[str, str]] = []
    for record in plan:
        source = source_stats[record.source_accession]
        rows.append(
            {
                "accession": record.accession,
                "is_new": record.is_new,
                "assembly_level": source.assembly_level,
                "genome_fasta": str(source.fasta_path),
            }
        )
    return rows


def build_metadata_rows(
    plan: Sequence[CohortRecord],
    source_stats: dict[str, SourceStats],
) -> list[dict[str, str]]:
    """Build the generated metadata table rows."""
    rows: list[dict[str, str]] = []
    for record in plan:
        if not record.include_metadata:
            continue
        source = source_stats[record.source_accession]
        rows.append(
            {
                "Accession": record.accession,
                "Tax_ID": source.tax_id,
                "Organism_Name": source.organism_name,
                "Assembly_Level": source.assembly_level,
                "N50": str(source.n50),
                "Scaffolds": str(source.scaffolds),
                "Genome_Size": str(source.genome_size),
                "Atypical_Warnings": record.atypical_warnings,
            }
        )
    return rows


def write_source_stats(path: Path, stats: dict[str, SourceStats]) -> None:
    """Write derived source-genome statistics for inspection."""
    rows = [
        {
            "source_accession": record.source_accession,
            "organism_name": record.organism_name,
            "tax_id": record.tax_id,
            "assembly_level": record.assembly_level,
            "scaffolds": str(record.scaffolds),
            "n50": str(record.n50),
            "genome_size": str(record.genome_size),
            "fasta_path": str(record.fasta_path),
            "sha256": record.sha256,
        }
        for record in sorted(stats.values(), key=lambda item: item.source_accession)
    ]
    write_tsv(path, SOURCE_STATS_COLUMNS, rows)


def prepare_cohort(
    *,
    work_root: Path,
    offline: bool,
    source_catalog_path: Path = DEFAULT_SOURCE_CATALOG,
    cohort_plan_path: Path = DEFAULT_COHORT_PLAN,
) -> PreparedCohort:
    """Prepare the generated acceptance cohort and its cached downloads."""
    downloads_dir = work_root / "downloads"
    generated_dir = work_root / "generated"
    checksums_tsv = work_root / "download_checksums.tsv"
    sample_csv = generated_dir / "sample_sheet.csv"
    metadata_tsv = generated_dir / "metadata.tsv"
    source_stats_tsv = generated_dir / "source_stats.tsv"

    catalog = load_source_catalog(source_catalog_path)
    plan = load_cohort_plan(cohort_plan_path, catalog)
    source_stats = build_source_stats(
        catalog,
        plan,
        downloads_dir=downloads_dir,
        offline=offline,
    )
    record_and_validate_checksums(tuple(source_stats.values()), checksums_tsv)
    write_csv(sample_csv, SAMPLE_COLUMNS, build_sample_rows(plan, source_stats))
    write_tsv(metadata_tsv, METADATA_COLUMNS, build_metadata_rows(plan, source_stats))
    write_source_stats(source_stats_tsv, source_stats)
    validate_prepared_cohort(plan, sample_csv, metadata_tsv, checksums_tsv)
    LOGGER.info("Prepared acceptance cohort under %s", generated_dir)
    return PreparedCohort(
        work_root=work_root,
        sample_csv=sample_csv,
        metadata_tsv=metadata_tsv,
        source_stats_tsv=source_stats_tsv,
        checksums_tsv=checksums_tsv,
        cohort_plan=tuple(plan),
        source_stats=source_stats,
    )


def prepare_medium_cohort(
    *,
    work_root: Path,
    offline: bool,
) -> PreparedCohort:
    """Prepare the fixed medium Mycoplasmatota/Bacillota cohort."""
    return prepare_cohort(
        work_root=work_root,
        offline=offline,
        source_catalog_path=DEFAULT_MEDIUM_SOURCE_CATALOG,
        cohort_plan_path=DEFAULT_MEDIUM_COHORT_PLAN,
    )


def validate_prepared_cohort(
    plan: Sequence[CohortRecord],
    sample_csv: Path,
    metadata_tsv: Path,
    checksums_tsv: Path,
) -> None:
    """Validate the generated cohort artefacts before pipeline execution."""
    sample_header, sample_rows = read_csv(sample_csv)
    if tuple(sample_header) != SAMPLE_COLUMNS:
        raise AcceptanceTestError("generated_sample_sheet_header_invalid")
    metadata_header, metadata_rows = read_tsv(metadata_tsv)
    if tuple(metadata_header) != METADATA_COLUMNS:
        raise AcceptanceTestError("generated_metadata_header_invalid")
    checksum_header, checksum_rows = read_tsv(checksums_tsv)
    if tuple(checksum_header) != DOWNLOAD_COLUMNS:
        raise AcceptanceTestError("download_checksums_header_invalid")

    expected_accessions = [record.accession for record in plan]
    if [row["accession"] for row in sample_rows] != expected_accessions:
        raise AcceptanceTestError("generated_sample_sheet_accessions_mismatch")

    metadata_accessions = {row["Accession"] for row in metadata_rows}
    for record in plan:
        has_metadata = record.accession in metadata_accessions
        if record.include_metadata != has_metadata:
            raise AcceptanceTestError(
                f"generated_metadata_presence_mismatch:{record.accession}"
            )
    if not checksum_rows:
        raise AcceptanceTestError("download_checksums_empty")


def require_command(name: str) -> None:
    """Raise when one required external command is unavailable."""
    if shutil.which(name) is None:
        raise AcceptanceTestError(f"Required command is not available: {name}")


def run_command(
    command: Sequence[str],
    *,
    cwd: Path,
    env: dict[str, str] | None = None,
) -> None:
    """Run one subprocess command and raise on failure."""
    display = " ".join(command)
    LOGGER.info("Running command: %s", display)
    completed = subprocess.run(command, cwd=cwd, env=env, check=False)
    if completed.returncode != 0:
        raise AcceptanceTestError(
            f"Command failed with exit code {completed.returncode}: {display}"
        )


def maybe_add_parameter(command: list[str], name: str, value: object | None) -> None:
    """Append one Nextflow-style parameter when it has a value."""
    if value is None or value == "":
        return
    text_value = str(value)
    if text_value.startswith("-"):
        command.append(f"{name}={text_value}")
        return
    command.extend([name, text_value])


def validate_real_run_args(args: argparse.Namespace) -> None:
    """Validate the required resource arguments for real-data runs."""
    missing: list[str] = []
    if not args.taxdump:
        missing.append("--taxdump")
    if not args.checkm2_db:
        missing.append("--checkm2-db")
    if not args.codetta_db:
        missing.append("--codetta-db")
    if not args.eggnog_db:
        missing.append("--eggnog-db")
    if not args.prepare_busco_datasets and not args.busco_db:
        missing.append("--busco-db or --prepare-busco-datasets")
    if missing:
        raise AcceptanceTestError(
            "Missing required arguments for real-data runs: " + ", ".join(missing)
        )


def validate_dbprep_run_args(args: argparse.Namespace) -> None:
    """Validate the required destination arguments for dbprep-slurm."""
    missing: list[str] = []
    if not args.taxdump:
        missing.append("--taxdump")
    if not args.checkm2_db:
        missing.append("--checkm2-db")
    if not args.codetta_db:
        missing.append("--codetta-db")
    if not args.busco_db:
        missing.append("--busco-db")
    if not args.eggnog_db:
        missing.append("--eggnog-db")
    if missing:
        raise AcceptanceTestError(
            "Missing required arguments for dbprep-slurm: " + ", ".join(missing)
        )


def build_nextflow_command(
    *,
    profile: str,
    work_dir: Path,
    outdir: Path,
    cohort: PreparedCohort,
    args: argparse.Namespace,
) -> list[str]:
    """Build the Nextflow command for one pipeline execution."""
    command = [
        "nextflow",
        "run",
        ".",
        "-profile",
        profile,
        "-work-dir",
        str(work_dir),
        "--sample_csv",
        str(cohort.sample_csv),
        "--metadata",
        str(cohort.metadata_tsv),
        "--taxdump",
        str(Path(args.taxdump).resolve()),
        "--checkm2_db",
        str(Path(args.checkm2_db).resolve()),
        "--codetta_db",
        str(Path(args.codetta_db).resolve()),
        "--eggnog_db",
        str(Path(args.eggnog_db).resolve()),
        "--outdir",
        str(outdir),
    ]
    if args.resume:
        command.append("-resume")
    if args.prepare_busco_datasets:
        command.extend(["--prepare_busco_datasets", "true"])
    if args.busco_db:
        command.extend(["--busco_db", str(Path(args.busco_db).resolve())])
    maybe_add_parameter(command, "--slurm_queue", args.slurm_queue)
    maybe_add_parameter(command, "--slurm_qos", args.slurm_qos)
    maybe_add_parameter(command, "--slurm_cluster_options", args.slurm_cluster_options)
    maybe_add_parameter(command, "--singularity_cache_dir", args.singularity_cache_dir)
    maybe_add_parameter(command, "--singularity_run_options", args.singularity_run_options)
    return command


def build_dbprep_command(
    *,
    profile: str,
    work_dir: Path,
    outdir: Path,
    args: argparse.Namespace,
) -> list[str]:
    """Build the Nextflow command for one runtime database prep execution."""
    command = [
        "nextflow",
        "run",
        "prepare_databases.nf",
        "-profile",
        profile,
        "-work-dir",
        str(work_dir),
        "--taxdump",
        str(Path(args.taxdump).resolve()),
        "--checkm2_db",
        str(Path(args.checkm2_db).resolve()),
        "--codetta_db",
        str(Path(args.codetta_db).resolve()),
        "--busco_db",
        str(Path(args.busco_db).resolve()),
        "--eggnog_db",
        str(Path(args.eggnog_db).resolve()),
        "--download_missing_databases",
        "true",
        "--outdir",
        str(outdir),
    ]
    if args.resume:
        command.append("-resume")
    if args.force_runtime_database_rebuild:
        command.extend(["--force_runtime_database_rebuild", "true"])
    maybe_add_parameter(command, "--slurm_queue", args.slurm_queue)
    maybe_add_parameter(command, "--slurm_qos", args.slurm_qos)
    maybe_add_parameter(command, "--slurm_cluster_options", args.slurm_cluster_options)
    maybe_add_parameter(command, "--singularity_cache_dir", args.singularity_cache_dir)
    maybe_add_parameter(command, "--singularity_run_options", args.singularity_run_options)
    return command


def read_sample_status_columns() -> list[str]:
    """Read the locked sample-status column order."""
    return [
        line.strip()
        for line in SAMPLE_STATUS_COLUMNS_ASSET.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def require_output_tables(outdir: Path) -> dict[str, Path]:
    """Require the published final tables for one pipeline run."""
    table_dir = outdir / "tables"
    missing = [name for name in FINAL_TABLES if not (table_dir / name).is_file()]
    if missing:
        raise AcceptanceTestError(
            "missing_output_tables:" + ",".join(sorted(missing))
        )
    return {name: table_dir / name for name in FINAL_TABLES}


def require_dbprep_outputs(outdir: Path) -> dict[str, Path]:
    """Require the published artefacts for one successful database prep run."""
    missing = [name for name in DBPREP_OUTPUTS if not (outdir / name).is_file()]
    if missing:
        raise AcceptanceTestError("missing_dbprep_outputs:" + ",".join(sorted(missing)))
    return {name: outdir / name for name in DBPREP_OUTPUTS}


def assert_dbprep_database_tree(
    taxdump_dir: Path,
    checkm2_dir: Path,
    codetta_dir: Path,
    busco_root: Path,
    eggnog_dir: Path,
) -> None:
    """Assert that the prepared runtime database tree is complete."""
    required_files = (
        taxdump_dir / "names.dmp",
        taxdump_dir / "nodes.dmp",
        eggnog_dir / "eggnog.db",
        eggnog_dir / "eggnog_proteins.dmnd",
        taxdump_dir / ".nf_myco_ready.json",
        checkm2_dir / ".nf_myco_ready.json",
        codetta_dir / "Pfam-A_enone.hmm",
        codetta_dir / "Pfam-A_enone.hmm.h3f",
        codetta_dir / "Pfam-A_enone.hmm.h3i",
        codetta_dir / "Pfam-A_enone.hmm.h3m",
        codetta_dir / "Pfam-A_enone.hmm.h3p",
        codetta_dir / ".nf_myco_ready.json",
        busco_root / ".nf_myco_ready.json",
        eggnog_dir / ".nf_myco_ready.json",
    )
    missing = [str(path) for path in required_files if not path.is_file()]
    if missing:
        raise AcceptanceTestError("missing_dbprep_paths:" + ",".join(sorted(missing)))

    dmnd_candidates = sorted(path for path in checkm2_dir.glob("*.dmnd") if path.is_file())
    if len(dmnd_candidates) != 1:
        raise AcceptanceTestError("invalid_dbprep_checkm2_dmnd_count")

    for lineage in DEFAULT_DBPREP_BUSCO_LINEAGES:
        lineage_dir = busco_root / lineage
        dataset_cfg = lineage_dir / "dataset.cfg"
        if not dataset_cfg.is_file():
            raise AcceptanceTestError(f"missing_dbprep_busco_lineage:{lineage}")


def assert_dbprep_report_contract(report_path: Path) -> None:
    """Assert that the merged prep report includes all required components."""
    _header, rows = read_tsv(report_path)
    seen_components = {row["component"] for row in rows}
    required_components = {
        "taxdump",
        "checkm2",
        "codetta",
        "busco_root",
        "eggnog",
    }
    missing = sorted(required_components - seen_components)
    if missing:
        raise AcceptanceTestError("missing_dbprep_report_rows:" + ",".join(missing))


def read_rows_indexed_by(
    path: Path,
    *,
    delimiter: str,
    key_column: str,
) -> tuple[list[str], dict[str, dict[str, str]]]:
    """Read one keyed table and index it by one unique column."""
    reader_fn = read_tsv if delimiter == "\t" else read_csv
    header, rows = reader_fn(path)
    index: dict[str, dict[str, str]] = {}
    for row in rows:
        key = row.get(key_column, "").strip()
        if not key:
            raise AcceptanceTestError(f"Table {path} contains an empty {key_column} value.")
        if key in index:
            raise AcceptanceTestError(f"Table {path} contains duplicate key {key!r}.")
        index[key] = row
    return header, index


def assert_metadata_contract(master_table: Path, metadata_tsv: Path) -> None:
    """Assert that the final master table preserves metadata and append contracts."""
    master_header, _master_rows = read_tsv(master_table)
    metadata_header, _metadata_rows = read_tsv(metadata_tsv)
    append_columns = master_table_contract.read_append_columns_asset()
    metadata_width = len(metadata_header)
    if master_header[:metadata_width] != metadata_header:
        raise AcceptanceTestError("metadata_block_order_mismatch")
    if master_header[metadata_width:] != append_columns:
        raise AcceptanceTestError("append_columns_mismatch")


def assert_sample_status_contract(sample_status: Path) -> None:
    """Assert that the final sample-status table matches the locked contract."""
    header, _rows = read_tsv(sample_status)
    if header != read_sample_status_columns():
        raise AcceptanceTestError("sample_status_columns_mismatch")


def parse_int(value: str) -> int | None:
    """Parse one optional integer-like value."""
    token = value.strip()
    if token in {"", "NA"}:
        return None
    try:
        return int(float(token))
    except ValueError:
        return None


def status_contains_token(value: str, token: str) -> bool:
    """Return True when one semicolon-delimited field contains a token."""
    return token in {part.strip() for part in value.split(";") if part.strip()}


def records_with_role(plan: Sequence[CohortRecord], role_tag: str) -> list[CohortRecord]:
    """Return all cohort records carrying one role tag."""
    return [record for record in plan if role_tag in record.role_tags]


def read_16s_manifest_accessions(path: Path, missing_column_error: str) -> set[str]:
    """Read accessions from one cohort 16S manifest."""
    header, rows = read_tsv(path)
    if "accession" not in header:
        raise AcceptanceTestError(missing_column_error)
    return {
        row.get("accession", "").strip()
        for row in rows
        if row.get("accession", "").strip()
    }


def assert_role_coverage(
    *,
    plan: Sequence[CohortRecord],
    master_rows: dict[str, dict[str, str]],
    status_rows: dict[str, dict[str, str]],
    intact_manifest_accessions: set[str],
    standard_manifest_accessions: set[str],
    low_quality_manifest_accessions: set[str],
) -> None:
    """Assert that the positive cohort proved all required acceptance roles."""
    gcode4_candidates = records_with_role(plan, "gcode4_candidate")
    if not any(master_rows[record.accession]["Gcode"] == "4" for record in gcode4_candidates):
        raise AcceptanceTestError("missing_gcode4_sample")

    gcode11_candidates = records_with_role(plan, "gcode11_candidate")
    if not any(master_rows[record.accession]["Gcode"] == "11" for record in gcode11_candidates):
        raise AcceptanceTestError("missing_gcode11_sample")

    crispr_positive_candidates = records_with_role(plan, "crispr_positive_candidate")
    if not any(
        (parse_int(master_rows[record.accession]["CRISPRS"]) or 0) > 0
        for record in crispr_positive_candidates
    ):
        raise AcceptanceTestError("missing_crispr_positive_sample")

    crispr_negative_candidates = records_with_role(plan, "crispr_negative_candidate")
    if not any(master_rows[record.accession]["CRISPRS"] == "0" for record in crispr_negative_candidates):
        raise AcceptanceTestError("missing_crispr_negative_sample")

    ani_cluster_candidates = records_with_role(plan, "ani_cluster_candidate")
    cluster_counts: dict[str, int] = {}
    for record in ani_cluster_candidates:
        cluster_id = master_rows[record.accession]["Cluster_ID"]
        if cluster_id == "NA":
            continue
        cluster_counts[cluster_id] = cluster_counts.get(cluster_id, 0) + 1
    if not any(count >= 2 for count in cluster_counts.values()):
        raise AcceptanceTestError("missing_ani_cluster_pair")

    low_quality_candidates = [
        record
        for record in plan
        if status_contains_token(
            status_rows[record.accession]["ani_exclusion_reason"],
            "low_quality",
        )
    ]
    if any(
        record.accession in standard_manifest_accessions
        for record in low_quality_candidates
    ):
        raise AcceptanceTestError("unexpected_low_quality_16s_inclusion")
    low_quality_16s_candidates = [
        record
        for record in low_quality_candidates
        if master_rows[record.accession].get("16S") in {"Yes", "partial"}
    ]
    if low_quality_16s_candidates and not any(
        record.accession in low_quality_manifest_accessions
        for record in low_quality_16s_candidates
    ):
        raise AcceptanceTestError("missing_low_quality_16s_manifest_inclusion")

    atypical_excluded_candidates = records_with_role(plan, "atypical_excluded_candidate")
    if not any(
        status_rows[record.accession]["ani_included"] == "false"
        and "atypical" in status_rows[record.accession]["ani_exclusion_reason"]
        for record in atypical_excluded_candidates
    ):
        raise AcceptanceTestError("missing_atypical_exclusion")
    if any(
        record.accession in intact_manifest_accessions
        for record in atypical_excluded_candidates
    ):
        raise AcceptanceTestError("unexpected_atypical_excluded_16s_inclusion")

    atypical_exception_candidates = records_with_role(plan, "atypical_exception_candidate")
    if not any(
        status_rows[record.accession]["ani_included"] == "true"
        for record in atypical_exception_candidates
    ):
        raise AcceptanceTestError("missing_atypical_exception_inclusion")
    if not any(
        record.accession in intact_manifest_accessions
        for record in atypical_exception_candidates
    ):
        raise AcceptanceTestError("missing_atypical_exception_16s_inclusion")

    atypical_mixed_candidates = records_with_role(plan, "atypical_mixed_candidate")
    if not any(
        status_rows[record.accession]["ani_included"] == "false"
        and "atypical" in status_rows[record.accession]["ani_exclusion_reason"]
        for record in atypical_mixed_candidates
    ):
        raise AcceptanceTestError("missing_atypical_mixed_exclusion")
    if any(
        record.accession in intact_manifest_accessions
        for record in atypical_mixed_candidates
    ):
        raise AcceptanceTestError("unexpected_atypical_mixed_16s_inclusion")

    missing_metadata_candidates = records_with_role(plan, "missing_metadata_case")
    if not any(
        status_contains_token(status_rows[record.accession]["warnings"], "missing_metadata_for_new_sample")
        and master_rows[record.accession]["Tax_ID"] == "NA"
        for record in missing_metadata_candidates
    ):
        raise AcceptanceTestError("missing_missing_metadata_case")

    collision_candidates = records_with_role(plan, "collision_candidate")
    collision_rows = [status_rows[record.accession] for record in collision_candidates]
    if len(collision_rows) < 2:
        raise AcceptanceTestError("collision_candidate_underpopulated")
    internal_ids = {row["internal_id"] for row in collision_rows}
    bases = {
        validate_inputs.sanitise_accession(row["accession"])
        for row in collision_rows
    }
    if len(internal_ids) != len(collision_rows):
        raise AcceptanceTestError("collision_candidate_internal_ids_not_unique")
    if len(bases) != 1:
        raise AcceptanceTestError("collision_candidate_sanitised_bases_mismatch")
    if not all(
        status_contains_token(row["warnings"], "internal_id_collision_resolved")
        for row in collision_rows
    ):
        raise AcceptanceTestError("missing_collision_warning")


def assert_run_outputs(outdir: Path, plan: Sequence[CohortRecord], metadata_tsv: Path) -> dict[str, Path]:
    """Assert the published outputs for one successful pipeline run."""
    tables = require_output_tables(outdir)
    assert_metadata_contract(tables["master_table.tsv"], metadata_tsv)
    assert_sample_status_contract(tables["sample_status.tsv"])

    _master_header, master_rows = read_rows_indexed_by(
        tables["master_table.tsv"],
        delimiter="\t",
        key_column="Accession",
    )
    _status_header, status_rows = read_rows_indexed_by(
        tables["sample_status.tsv"],
        delimiter="\t",
        key_column="accession",
    )
    if len(master_rows) != len(plan):
        raise AcceptanceTestError("master_table_row_count_mismatch")
    if len(status_rows) != len(plan):
        raise AcceptanceTestError("sample_status_row_count_mismatch")

    missing_sample_dirs = [
        record.accession
        for record in plan
        if not (outdir / "samples" / record.accession).is_dir()
    ]
    if missing_sample_dirs:
        raise AcceptanceTestError(
            "missing_sample_folders:" + ",".join(missing_sample_dirs)
        )

    cohort_16s_dir = outdir / "cohort" / "16s"
    intact_manifest_accessions = read_16s_manifest_accessions(
        cohort_16s_dir / "all_best_16S_manifest.tsv",
        "missing_intact_16s_manifest_accession_column",
    )
    partial_manifest_accessions = read_16s_manifest_accessions(
        cohort_16s_dir / "all_partial_16S_manifest.tsv",
        "missing_partial_16s_manifest_accession_column",
    )
    low_quality_best_manifest_accessions = read_16s_manifest_accessions(
        cohort_16s_dir / "low_quality_best_16S_manifest.tsv",
        "missing_low_quality_best_16s_manifest_accession_column",
    )
    low_quality_partial_manifest_accessions = read_16s_manifest_accessions(
        cohort_16s_dir / "low_quality_partial_16S_manifest.tsv",
        "missing_low_quality_partial_16s_manifest_accession_column",
    )

    assert_role_coverage(
        plan=plan,
        master_rows=master_rows,
        status_rows=status_rows,
        intact_manifest_accessions=intact_manifest_accessions,
        standard_manifest_accessions=intact_manifest_accessions | partial_manifest_accessions,
        low_quality_manifest_accessions=(
            low_quality_best_manifest_accessions | low_quality_partial_manifest_accessions
        ),
    )
    return tables


def compare_files_exact(local_path: Path, slurm_path: Path, error_code: str) -> None:
    """Compare two files byte-for-byte."""
    if local_path.read_bytes() != slurm_path.read_bytes():
        raise AcceptanceTestError(error_code)


def load_versions_rows(path: Path) -> list[dict[str, str]]:
    """Load one versions TSV."""
    _header, rows = read_tsv(path)
    return rows


def compare_versions_logically(local_path: Path, slurm_path: Path) -> None:
    """Compare stable provenance rows while ignoring runtime-only differences."""
    local_rows = {
        (
            row["component"],
            row["kind"],
            row["version"],
            row["image_or_path"],
            row["notes"],
        )
        for row in load_versions_rows(local_path)
        if row["kind"] in {"pipeline", "database", "container"}
    }
    slurm_rows = {
        (
            row["component"],
            row["kind"],
            row["version"],
            row["image_or_path"],
            row["notes"],
        )
        for row in load_versions_rows(slurm_path)
        if row["kind"] in {"pipeline", "database", "container"}
    }
    if local_rows != slurm_rows:
        raise AcceptanceTestError("slurm_versions_logical_differs_from_local")


def run_unit_tests() -> None:
    """Run the repository's Python unit-test suite."""
    run_command(
        [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-p", "test_*.py"],
        cwd=ROOT_DIR,
    )


def run_stub_smoke(args: argparse.Namespace) -> None:
    """Run the full-pipeline stub smoke test."""
    require_command("nextflow")
    run_dir = args.work_root / "runs" / "stub"
    outdir = run_dir / "results"
    work_dir = run_dir / "work"
    command = [
        "nextflow",
        "run",
        ".",
        "-profile",
        "test",
        "-stub-run",
        "-work-dir",
        str(work_dir),
        "--outdir",
        str(outdir),
    ]
    if args.resume:
        command.append("-resume")
    run_command(command, cwd=ROOT_DIR)
    require_output_tables(outdir)


def run_real_pipeline(args: argparse.Namespace, *, stage_name: str, profile: str) -> dict[str, Path]:
    """Prepare the cohort, run the real pipeline, and assert outputs."""
    validate_real_run_args(args)
    require_command("nextflow")
    if stage_name == "slurm":
        require_command("sbatch")

    cohort = prepare_cohort(
        work_root=args.work_root,
        offline=args.offline,
        source_catalog_path=args.source_catalog,
        cohort_plan_path=args.cohort_plan,
    )
    run_dir = cohort.work_root / "runs" / stage_name
    outdir = run_dir / "results"
    work_dir = run_dir / "work"
    command = build_nextflow_command(
        profile=profile,
        work_dir=work_dir,
        outdir=outdir,
        cohort=cohort,
        args=args,
    )
    run_command(command, cwd=ROOT_DIR)
    return assert_run_outputs(outdir, cohort.cohort_plan, cohort.metadata_tsv)


def run_slurm_with_comparison(args: argparse.Namespace) -> None:
    """Run the SLURM pipeline and compare stable outputs with the local baseline."""
    local_tables_dir = args.work_root / "runs" / "local" / "results" / "tables"
    if not local_tables_dir.is_dir():
        raise AcceptanceTestError("missing_local_baseline")
    slurm_tables = run_real_pipeline(args, stage_name="slurm", profile=args.slurm_profile)
    compare_files_exact(
        local_tables_dir / "master_table.tsv",
        slurm_tables["master_table.tsv"],
        "slurm_master_table_differs_from_local",
    )
    compare_files_exact(
        local_tables_dir / "sample_status.tsv",
        slurm_tables["sample_status.tsv"],
        "slurm_sample_status_differs_from_local",
    )
    compare_versions_logically(
        local_tables_dir / "tool_and_db_versions.tsv",
        slurm_tables["tool_and_db_versions.tsv"],
    )


def run_dbprep_slurm(args: argparse.Namespace) -> None:
    """Run the runtime database prep workflow on SLURM and validate outputs."""
    validate_dbprep_run_args(args)
    require_command("nextflow")
    require_command("sbatch")

    run_dir = args.work_root / "runs" / "dbprep-slurm"
    outdir = run_dir / "results"
    work_dir = run_dir / "work"

    command = build_dbprep_command(
        profile=args.dbprep_profile,
        work_dir=work_dir,
        outdir=outdir,
        args=args,
    )
    run_command(command, cwd=ROOT_DIR)
    outputs = require_dbprep_outputs(outdir)
    assert_dbprep_database_tree(
        Path(args.taxdump).resolve(),
        Path(args.checkm2_db).resolve(),
        Path(args.codetta_db).resolve(),
        Path(args.busco_db).resolve(),
        Path(args.eggnog_db).resolve(),
    )
    assert_dbprep_report_contract(outputs["runtime_database_report.tsv"])


def build_common_parser() -> argparse.ArgumentParser:
    """Build the parent parser for shared arguments."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "--work-root",
        type=Path,
        default=DEFAULT_WORK_ROOT,
        help="Directory for generated acceptance inputs and run artefacts.",
    )
    parser.add_argument(
        "--offline",
        action="store_true",
        help="Do not download missing source genomes during prepare.",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Pass -resume to Nextflow runs.",
    )
    parser.add_argument(
        "--source-catalog",
        type=Path,
        default=DEFAULT_SOURCE_CATALOG,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--cohort-plan",
        type=Path,
        default=DEFAULT_COHORT_PLAN,
        help=argparse.SUPPRESS,
    )
    return parser


def build_real_run_parser() -> argparse.ArgumentParser:
    """Build the parent parser for real-data run arguments."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--taxdump", type=Path, default=None, help="Pinned taxdump directory.")
    parser.add_argument("--checkm2-db", type=Path, default=None, help="CheckM2 database path.")
    parser.add_argument(
        "--codetta-db",
        type=Path,
        default=None,
        help="Codetta profile database path.",
    )
    parser.add_argument(
        "--busco-db",
        type=Path,
        default=None,
        help="Pre-downloaded BUSCO lineage directory root.",
    )
    parser.add_argument(
        "--prepare-busco-datasets",
        action="store_true",
        help="Allow the pipeline to prepare BUSCO lineage datasets itself.",
    )
    parser.add_argument("--eggnog-db", type=Path, default=None, help="eggNOG database path.")
    parser.add_argument(
        "--local-profile",
        default=DEFAULT_LOCAL_PROFILE,
        help="Nextflow profile string for local real-data runs.",
    )
    parser.add_argument(
        "--slurm-profile",
        default=DEFAULT_SLURM_PROFILE,
        help="Nextflow profile string for SLURM real-data runs.",
    )
    parser.add_argument("--slurm-queue", default=None, help="Optional SLURM queue.")
    parser.add_argument("--slurm-qos", default=None, help="Optional SLURM QoS.")
    parser.add_argument(
        "--slurm-cluster-options",
        default=None,
        help="Optional extra SLURM cluster options.",
    )
    parser.add_argument(
        "--singularity-cache-dir",
        default=None,
        help="Optional Singularity cache directory.",
    )
    parser.add_argument(
        "--singularity-run-options",
        default=None,
        help="Optional extra Singularity run options.",
    )
    return parser


def build_dbprep_run_parser() -> argparse.ArgumentParser:
    """Build the parent parser for SLURM database prep arguments."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "--dbprep-profile",
        default=DEFAULT_DBPREP_PROFILE,
        help="Nextflow profile string for SLURM database prep runs.",
    )
    parser.add_argument("--taxdump", type=Path, default=None, help="Pinned taxdump directory.")
    parser.add_argument("--checkm2-db", type=Path, default=None, help="CheckM2 database path.")
    parser.add_argument("--codetta-db", type=Path, default=None, help="Codetta profile database path.")
    parser.add_argument("--busco-db", type=Path, default=None, help="BUSCO database root.")
    parser.add_argument("--eggnog-db", type=Path, default=None, help="eggNOG database path.")
    parser.add_argument(
        "--force-runtime-database-rebuild",
        action="store_true",
        help="Rebuild incomplete runtime database destinations in place.",
    )
    parser.add_argument("--slurm-queue", default=None, help="Optional SLURM queue.")
    parser.add_argument("--slurm-qos", default=None, help="Optional SLURM QoS.")
    parser.add_argument(
        "--slurm-cluster-options",
        default=None,
        help="Optional extra SLURM cluster options.",
    )
    parser.add_argument(
        "--singularity-cache-dir",
        default=None,
        help="Optional Singularity cache directory.",
    )
    parser.add_argument(
        "--singularity-run-options",
        default=None,
        help="Optional extra Singularity run options.",
    )
    return parser


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    common_parser = build_common_parser()
    real_parser = build_real_run_parser()
    parser = argparse.ArgumentParser(
        description="Run layered acceptance tests for local and SLURM execution paths."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    dbprep_parser = build_dbprep_run_parser()
    subparsers.add_parser("prepare", parents=[common_parser], help="Prepare the acceptance cohort.")
    subparsers.add_parser("unit", parents=[common_parser], help="Run the Python unit-test layer.")
    subparsers.add_parser("stub", parents=[common_parser], help="Run the stub Nextflow smoke test.")
    subparsers.add_parser(
        "local",
        parents=[common_parser, real_parser],
        help="Run the real-data local acceptance cohort.",
        description=(
            "Run the real-data local acceptance cohort. "
            f"{REAL_RUN_NOTE}"
        ),
    )
    subparsers.add_parser(
        "slurm",
        parents=[common_parser, real_parser],
        help="Run the real-data SLURM acceptance cohort.",
        description=(
            "Run the real-data SLURM acceptance cohort. "
            f"{REAL_RUN_NOTE}"
        ),
    )
    subparsers.add_parser(
        "dbprep-slurm",
        parents=[common_parser, dbprep_parser],
        help="Run the runtime database prep workflow on SLURM.",
        description=DBPREP_RUN_NOTE,
    )
    subparsers.add_parser(
        "all",
        parents=[common_parser, real_parser],
        help="Run prepare, unit, stub, local, and SLURM layers in sequence.",
        description=(
            "Run prepare, unit, stub, local, and SLURM layers in sequence. "
            f"{REAL_RUN_NOTE}"
        ),
    )
    return parser.parse_args(argv)


def dispatch(args: argparse.Namespace) -> None:
    """Dispatch the requested acceptance-test subcommand."""
    if args.command == "prepare":
        prepare_cohort(
            work_root=args.work_root,
            offline=args.offline,
            source_catalog_path=args.source_catalog,
            cohort_plan_path=args.cohort_plan,
        )
        return
    if args.command == "unit":
        run_unit_tests()
        return
    if args.command == "stub":
        run_stub_smoke(args)
        return
    if args.command == "local":
        run_real_pipeline(args, stage_name="local", profile=args.local_profile)
        return
    if args.command == "slurm":
        run_slurm_with_comparison(args)
        return
    if args.command == "dbprep-slurm":
        run_dbprep_slurm(args)
        return
    if args.command == "all":
        prepare_cohort(
            work_root=args.work_root,
            offline=args.offline,
            source_catalog_path=args.source_catalog,
            cohort_plan_path=args.cohort_plan,
        )
        run_unit_tests()
        run_stub_smoke(args)
        run_real_pipeline(args, stage_name="local", profile=args.local_profile)
        run_slurm_with_comparison(args)
        return
    raise AcceptanceTestError(f"Unsupported command: {args.command}")


def main(argv: Sequence[str] | None = None) -> int:
    """Run the layered acceptance-test harness."""
    args = parse_args(argv)
    configure_logging()
    try:
        dispatch(args)
    except AcceptanceTestError as error:
        LOGGER.error(str(error))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
