#!/usr/bin/env python3
"""Prepare runtime databases in operator-managed destinations."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import os
import shlex
import shutil
import subprocess
import sys
import tarfile
import zipfile
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Callable, Iterator, Sequence
from urllib.parse import urlparse

from annotation_resources import (
    RESOURCE_TOOLS,
    AnnotationResourceError,
    prepare_resource,
    validate_resource,
)
from download_annotation_resources import ResourceDownloadError, acquire_component

LOGGER = logging.getLogger(__name__)
MARKER_FILE_NAME = ".nf_myco_ready.json"
LOCK_FILE_SUFFIX = ".nf_myco_prepare.lock"
DEFAULT_BUSCO_LINEAGES = ("bacillota_odb12", "mycoplasmatota_odb12")
REPORT_COLUMNS = ("component", "status", "source", "destination", "details")
LINK_MODES = ("copy", "symlink", "hardlink")
MAX_MD5_DOWNLOAD_ATTEMPTS = 3
DEFAULT_REMOTE_SOURCE_MANIFEST = Path(
    os.environ.get(
        "NF_MYCO_RUNTIME_DB_SOURCE_MANIFEST",
        Path(__file__).resolve().parents[1]
        / "assets"
        / "runtime"
        / "runtime_database_sources.json",
    )
)


class PrepareRuntimeDatabasesError(RuntimeError):
    """Raised when one runtime database cannot be prepared safely."""


class ChecksumMismatchError(PrepareRuntimeDatabasesError):
    """Raised when one downloaded file checksum does not match."""

    def __init__(
        self,
        *,
        path: Path,
        checksum_type: str,
        expected: str,
        observed: str,
    ) -> None:
        """Initialise one checksum-mismatch error."""
        self.path = path
        self.checksum_type = checksum_type
        self.expected = expected
        self.observed = observed
        super().__init__(
            f"Checksum mismatch for {path.name}: expected {expected}, observed {observed}"
        )


@dataclass(frozen=True)
class ValidationResult:
    """Describe one validated runtime database layout."""

    required_paths: tuple[str, ...]
    details: dict[str, str]


@dataclass(frozen=True)
class PreparationRecord:
    """Summarise one prepared or reused database destination."""

    component: str
    status: str
    source: str | None
    destination: Path
    details: str


Validator = Callable[[Path], ValidationResult]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Prepare runtime databases in operator-managed destinations without "
            "relying on system tmp."
        )
    )
    parser.add_argument("--taxdump-source", type=Path, help="Pinned taxdump source path.")
    parser.add_argument("--taxdump-dest", type=Path, help="Final taxdump destination.")
    parser.add_argument("--taxdump-version", default=None, help="Pinned remote taxdump version.")
    parser.add_argument("--checkm2-source", type=Path, help="CheckM2 database source path.")
    parser.add_argument("--checkm2-dest", type=Path, help="Final CheckM2 destination.")
    parser.add_argument("--checkm2-version", default=None, help="Pinned remote CheckM2 version.")
    parser.add_argument("--codetta-source", type=Path, help="Codetta profile database source path.")
    parser.add_argument("--codetta-dest", type=Path, help="Final Codetta profile database destination.")
    parser.add_argument("--codetta-version", default=None, help="Pinned remote Codetta profile database version.")
    parser.add_argument(
        "--busco-lineage-source",
        action="append",
        default=[],
        help="BUSCO lineage mapping as LINEAGE=PATH. May be supplied multiple times.",
    )
    parser.add_argument(
        "--busco-dest-root",
        type=Path,
        help="Final BUSCO lineage root directory.",
    )
    parser.add_argument("--busco-version", default=None, help="Pinned remote BUSCO version.")
    parser.add_argument("--eggnog-source", type=Path, help="eggNOG database source path.")
    parser.add_argument("--eggnog-dest", type=Path, help="Final eggNOG destination.")
    parser.add_argument("--eggnog-version", default=None, help="Pinned remote eggNOG version.")
    for component in ("cogclassifier", "pfam", "kofam"):
        parser.add_argument(f"--{component}-source", type=Path, help=f"Acquired or prepared {component} resource directory.")
        parser.add_argument(f"--{component}-dest", type=Path, help=f"Immutable {component} resource destination.")
        parser.add_argument(f"--{component}-version", default=None, help=f"Pinned {component} source version.")
    parser.add_argument("--padloc-source", type=Path, help="PADLOC database source path.")
    parser.add_argument("--padloc-dest", type=Path, help="Final PADLOC destination.")
    parser.add_argument("--padloc-version", default=None, help="Pinned remote PADLOC version.")
    parser.add_argument(
        "--download",
        action="store_true",
        help="Allow missing sources to be downloaded from curated remote sources.",
    )
    parser.add_argument(
        "--remote-source-manifest",
        type=Path,
        default=DEFAULT_REMOTE_SOURCE_MANIFEST,
        help="JSON manifest describing curated remote database sources.",
    )
    parser.add_argument(
        "--scratch-root",
        type=Path,
        default=None,
        help=(
            "Optional scratch root for downloads and archive extraction. Defaults "
            "to the destination parent directory."
        ),
    )
    parser.add_argument(
        "--link-mode",
        choices=LINK_MODES,
        default="copy",
        help="How directory sources are materialised into final destinations.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace invalid or incomplete destinations.",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=None,
        help="Optional TSV report path describing the prepared destinations.",
    )
    return parser.parse_args(argv)


def configure_logging() -> None:
    """Configure process logging."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def normalise_path(path: Path) -> Path:
    """Expand one user path without requiring it to exist yet."""
    return path.expanduser().resolve(strict=False)


def timestamp_token() -> str:
    """Return one UTC timestamp token suitable for file names."""
    return datetime.now(UTC).strftime("%Y%m%dT%H%M%SZ")


def marker_path(destination: Path) -> Path:
    """Return the ready-marker path for one prepared destination."""
    return destination / MARKER_FILE_NAME


def lock_path(destination: Path) -> Path:
    """Return the sidecar lock path for one destination."""
    return destination.parent / f".{destination.name}{LOCK_FILE_SUFFIX}"


def report_row(record: PreparationRecord) -> dict[str, str]:
    """Build one TSV report row."""
    return {
        "component": record.component,
        "status": record.status,
        "source": record.source if record.source else "NA",
        "destination": str(record.destination),
        "details": record.details,
    }


def parse_name_path(token: str, argument_name: str) -> tuple[str, Path]:
    """Parse one NAME=PATH token."""
    if "=" not in token:
        raise PrepareRuntimeDatabasesError(
            f"{argument_name} entries must use LINEAGE=PATH syntax: {token!r}"
        )
    name, raw_path = token.split("=", 1)
    clean_name = name.strip()
    if not clean_name:
        raise PrepareRuntimeDatabasesError(
            f"{argument_name} entry is missing the lineage before '=': {token!r}"
        )
    clean_path = raw_path.strip()
    if not clean_path:
        raise PrepareRuntimeDatabasesError(
            f"{argument_name} entry is missing the source path after '=': {token!r}"
        )
    return clean_name, normalise_path(Path(clean_path))


def ensure_parent_directory(path: Path) -> None:
    """Create the parent directory for one destination or report."""
    path.parent.mkdir(parents=True, exist_ok=True)


def remove_path(path: Path) -> None:
    """Remove one filesystem path regardless of type."""
    if not path.exists() and not path.is_symlink():
        return
    if path.is_symlink() or path.is_file():
        path.unlink()
        return
    shutil.rmtree(path)


def destination_is_empty(path: Path) -> bool:
    """Return whether one existing destination directory is empty."""
    return not any(path.iterdir())


def load_marker(destination: Path, component: str) -> dict[str, object] | None:
    """Load and validate one ready marker."""
    path = marker_path(destination)
    if not path.is_file():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise PrepareRuntimeDatabasesError(
            f"Marker file is not valid JSON for {component}: {path}"
        ) from error
    if not isinstance(data, dict):
        raise PrepareRuntimeDatabasesError(
            f"Marker file must contain a JSON object for {component}: {path}"
        )
    if data.get("component") != component:
        raise PrepareRuntimeDatabasesError(
            f"Marker component mismatch for {component}: {path}"
        )
    return data


def marker_exists(destination: Path) -> bool:
    """Return whether one destination already has a ready marker file."""
    return marker_path(destination).is_file()


def write_marker(
    *,
    component: str,
    source: str,
    destination: Path,
    validation: ValidationResult,
    source_metadata: dict[str, str] | None = None,
) -> None:
    """Write one ready marker after validation has succeeded."""
    payload = {
        "component": component,
        "source": source,
        "destination": str(destination),
        "prepared_at": datetime.now(UTC).isoformat(),
        "required_paths": list(validation.required_paths),
        "details": validation.details,
        "source_metadata": source_metadata or {},
    }
    marker_path(destination).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


@contextmanager
def preparation_lock(destination: Path, *, force: bool) -> Iterator[Path]:
    """Guard one destination against concurrent preparation."""
    sidecar = lock_path(destination)
    ensure_parent_directory(sidecar)
    if sidecar.exists():
        if not force:
            raise PrepareRuntimeDatabasesError(
                f"Preparation lock already exists for {destination}: {sidecar}"
            )
        LOGGER.warning("Removing existing preparation lock for %s", destination)
        remove_path(sidecar)
    payload = {
        "destination": str(destination),
        "pid": os.getpid(),
        "created_at": datetime.now(UTC).isoformat(),
    }
    sidecar.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    try:
        yield sidecar
    finally:
        sidecar.unlink(missing_ok=True)


def validate_taxdump(path: Path) -> ValidationResult:
    """Validate one taxdump directory."""
    required = ("names.dmp", "nodes.dmp")
    missing = [name for name in required if not (path / name).is_file()]
    if missing:
        raise PrepareRuntimeDatabasesError(
            f"Taxdump directory is missing required files in {path}: {', '.join(missing)}"
        )
    return ValidationResult(required_paths=required, details={"files": str(len(required))})


def validate_checkm2(path: Path) -> ValidationResult:
    """Validate one CheckM2 database directory."""
    candidates = sorted(file.name for file in path.glob("*.dmnd") if file.is_file())
    if len(candidates) != 1:
        raise PrepareRuntimeDatabasesError(
            f"CheckM2 destination must contain exactly one top-level .dmnd file: {path}"
        )
    return ValidationResult(required_paths=(candidates[0],), details={"dmnd": candidates[0]})


def validate_eggnog(path: Path) -> ValidationResult:
    """Validate a prepared v3-compatible eggNOG resource and content identity."""
    return validate_annotation_resource(path, "eggnog")


def validate_annotation_resource(path: Path, component: str) -> ValidationResult:
    """Expose the annotation content contract to the canonical ready-marker API."""
    try:
        resource = validate_resource(path, component)
    except AnnotationResourceError as error:
        raise PrepareRuntimeDatabasesError(str(error)) from error
    return ValidationResult(
        required_paths=tuple(entry["path"] for entry in resource["contract"]["files"]),
        details={"resource_id": resource["resource_id"], "version": resource["contract"]["version"]},
    )


def validate_codetta(path: Path) -> ValidationResult:
    """Validate one Codetta profile database directory."""
    required = (
        "Pfam-A_enone.hmm",
        "Pfam-A_enone.hmm.h3f",
        "Pfam-A_enone.hmm.h3i",
        "Pfam-A_enone.hmm.h3m",
        "Pfam-A_enone.hmm.h3p",
    )
    missing = [name for name in required if not (path / name).is_file()]
    if missing:
        raise PrepareRuntimeDatabasesError(
            f"Codetta destination is missing required files in {path}: {', '.join(missing)}"
        )
    return ValidationResult(
        required_paths=required,
        details={"profile": "Pfam-A_enone.hmm", "files": str(len(required))},
    )


def validate_padloc(path: Path) -> ValidationResult:
    """Validate one PADLOC data directory."""
    required = ("hmm/padlocdb.hmm",)
    if not (path / required[0]).is_file():
        raise PrepareRuntimeDatabasesError(
            f"PADLOC destination must contain hmm/padlocdb.hmm: {path}"
        )
    return ValidationResult(required_paths=required, details={"hmm": "hmm/padlocdb.hmm"})


def build_busco_lineage_validator(lineage: str) -> Validator:
    """Build a BUSCO lineage validator for one lineage directory."""

    def validate_busco_lineage(path: Path) -> ValidationResult:
        """Validate one BUSCO lineage directory."""
        required = ("dataset.cfg",)
        if not (path / required[0]).is_file():
            raise PrepareRuntimeDatabasesError(
                f"BUSCO lineage {lineage} is missing dataset.cfg in {path}"
            )
        return ValidationResult(required_paths=required, details={"lineage": lineage})

    return validate_busco_lineage


def build_busco_validator(lineages: Sequence[str]) -> Validator:
    """Build a BUSCO root validator for one lineage set."""

    def validate_busco_root(path: Path) -> ValidationResult:
        """Validate one BUSCO lineage root directory."""
        required = tuple(f"{lineage}/dataset.cfg" for lineage in lineages)
        missing = [entry for entry in required if not (path / entry).is_file()]
        if missing:
            raise PrepareRuntimeDatabasesError(
                f"BUSCO destination is missing required lineage files in {path}: "
                + ", ".join(missing)
            )
        return ValidationResult(
            required_paths=required,
            details={"lineages": ";".join(lineages)},
        )

    return validate_busco_root


def resolve_directory_source(path: Path, validator: Validator, component: str) -> Path:
    """Resolve one directory source to the single valid content root."""
    candidates = [path]
    candidates.extend(sorted(child for child in path.iterdir() if child.is_dir()))
    matches: list[Path] = []
    for candidate in candidates:
        try:
            validator(candidate)
        except PrepareRuntimeDatabasesError:
            continue
        matches.append(candidate)
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise PrepareRuntimeDatabasesError(
            f"Could not find a valid {component} root inside {path}"
        )
    raise PrepareRuntimeDatabasesError(
        f"Found multiple valid {component} roots inside {path}; choose a more specific source path."
    )


def copy_directory_tree(source: Path, destination: Path) -> None:
    """Copy one directory tree into an existing destination directory."""
    for child in source.iterdir():
        target = destination / child.name
        if child.is_dir():
            shutil.copytree(child, target)
        else:
            shutil.copy2(child, target)


def hardlink_directory_tree(source: Path, destination: Path) -> None:
    """Materialise one directory tree with hardlinked files."""
    for child in source.iterdir():
        target = destination / child.name
        if child.is_dir():
            target.mkdir(parents=True, exist_ok=False)
            hardlink_directory_tree(child, target)
            continue
        os.link(child, target)


def symlink_directory_tree(source: Path, destination: Path) -> None:
    """Materialise one directory tree with top-level symlinks."""
    for child in source.iterdir():
        target = destination / child.name
        os.symlink(child, target, target_is_directory=child.is_dir())


def materialise_directory(source: Path, destination: Path, link_mode: str) -> None:
    """Materialise one directory source into the final destination."""
    destination.mkdir(parents=True, exist_ok=True)
    if not destination_is_empty(destination):
        raise PrepareRuntimeDatabasesError(
            f"Destination must be absent or empty before preparation: {destination}"
        )
    if link_mode == "copy":
        copy_directory_tree(source, destination)
        return
    if link_mode == "hardlink":
        hardlink_directory_tree(source, destination)
        return
    symlink_directory_tree(source, destination)


def materialise_checkm2_file(source: Path, destination: Path, link_mode: str) -> None:
    """Materialise one raw CheckM2 .dmnd source into the final destination."""
    destination.mkdir(parents=True, exist_ok=True)
    if not destination_is_empty(destination):
        raise PrepareRuntimeDatabasesError(
            f"Destination must be absent or empty before preparation: {destination}"
        )
    target = destination / source.name
    if link_mode == "copy":
        shutil.copy2(source, target)
        return
    if link_mode == "hardlink":
        os.link(source, target)
        return
    os.symlink(source, target)


def ensure_same_filesystem(scratch_root: Path, destination_parent: Path, component: str) -> None:
    """Require archive scratch and destination to share one filesystem."""
    if scratch_root.stat().st_dev != destination_parent.stat().st_dev:
        raise PrepareRuntimeDatabasesError(
            f"{component} archive scratch root must be on the same filesystem as "
            f"the destination parent for rename-based finalisation: "
            f"{scratch_root} -> {destination_parent}"
        )


def extract_zip(archive_path: Path, destination: Path) -> None:
    """Extract one zip archive safely."""
    destination_root = destination.resolve()
    with zipfile.ZipFile(archive_path) as archive:
        for member in archive.infolist():
            target = (destination / member.filename).resolve()
            if not target.is_relative_to(destination_root):
                raise PrepareRuntimeDatabasesError(
                    f"Zip archive contains an unsafe path: {member.filename!r}"
                )
        archive.extractall(destination)


def extract_tar(archive_path: Path, destination: Path) -> None:
    """Extract one tar archive safely."""
    with tarfile.open(archive_path) as archive:
        archive.extractall(destination, filter="data")


def extract_archive(archive_path: Path, scratch_dir: Path) -> None:
    """Extract one supported archive into scratch space."""
    if zipfile.is_zipfile(archive_path):
        extract_zip(archive_path, scratch_dir)
        return
    if tarfile.is_tarfile(archive_path):
        extract_tar(archive_path, scratch_dir)
        return
    raise PrepareRuntimeDatabasesError(f"Unsupported archive format: {archive_path}")


def prepare_from_archive(
    *,
    component: str,
    source: Path,
    destination: Path,
    validator: Validator,
    scratch_root: Path | None,
) -> ValidationResult:
    """Prepare one database from an archive beside the final destination."""
    destination_parent = destination.parent
    destination_parent.mkdir(parents=True, exist_ok=True)
    scratch_base = normalise_path(scratch_root) if scratch_root else destination_parent
    scratch_base.mkdir(parents=True, exist_ok=True)
    ensure_same_filesystem(scratch_base, destination_parent, component)

    scratch_dir = scratch_base / f".prepare-{destination.name}-{timestamp_token()}"
    if scratch_dir.exists():
        raise PrepareRuntimeDatabasesError(
            f"Archive scratch directory already exists: {scratch_dir}"
        )

    try:
        scratch_dir.mkdir(parents=True, exist_ok=False)
        extract_archive(source, scratch_dir)
        resolved_root = resolve_directory_source(scratch_dir, validator, component)
        validation = validator(resolved_root)
        if destination.exists():
            if destination.is_dir() and destination_is_empty(destination):
                destination.rmdir()
            else:
                raise PrepareRuntimeDatabasesError(
                    f"Archive destination must be absent or empty before finalisation: "
                    f"{destination}"
                )
        resolved_root.rename(destination)
        if scratch_dir.exists():
            shutil.rmtree(scratch_dir)
        return validation
    except Exception:
        shutil.rmtree(scratch_dir, ignore_errors=True)
        raise


def describe_validation(
    validation: ValidationResult,
    extras: dict[str, str] | None = None,
) -> str:
    """Render one validation summary into a short report string."""
    paths = ",".join(validation.required_paths)
    details = dict(validation.details)
    if extras:
        details.update(extras)
    rendered_details = ";".join(f"{key}={value}" for key, value in sorted(details.items()))
    return f"required_paths={paths};{rendered_details}"


def assess_destination(
    *,
    component: str,
    destination: Path,
    validator: Validator,
) -> tuple[str, ValidationResult | None, dict[str, object] | None]:
    """Classify one destination before preparation."""
    if not destination.exists():
        return "missing", None, None
    if not destination.is_dir():
        raise PrepareRuntimeDatabasesError(
            f"Destination must be a directory for {component}: {destination}"
        )
    try:
        validation = validator(destination)
    except PrepareRuntimeDatabasesError as error:
        if marker_exists(destination):
            raise PrepareRuntimeDatabasesError(
                f"Destination marker exists but validation failed for {component}: {error}"
            ) from error
        raise PrepareRuntimeDatabasesError(
            f"Destination is invalid or partially prepared for {component}: {error}"
        ) from error
    marker = load_marker(destination, component)
    if marker is None:
        raise PrepareRuntimeDatabasesError(
            f"Destination is valid but incomplete for {component}; "
            f"missing {MARKER_FILE_NAME}: {destination}"
        )
    return "ready", validation, marker


def load_remote_source_manifest(path: Path) -> dict[str, Any]:
    """Load the curated remote source manifest."""
    manifest_path = normalise_path(path)
    if not manifest_path.is_file():
        raise PrepareRuntimeDatabasesError(
            f"Missing remote source manifest: {manifest_path}"
        )
    try:
        data = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise PrepareRuntimeDatabasesError(
            f"Remote source manifest is not valid JSON: {manifest_path}"
        ) from error
    if not isinstance(data, dict):
        raise PrepareRuntimeDatabasesError(
            f"Remote source manifest must contain a JSON object: {manifest_path}"
        )
    return data


def resolve_remote_version(
    *,
    manifest: dict[str, Any],
    component: str,
    requested_version: str | None,
) -> tuple[str, dict[str, Any]]:
    """Resolve one remote component version from the manifest."""
    component_manifest = manifest.get(component)
    if not isinstance(component_manifest, dict):
        raise PrepareRuntimeDatabasesError(
            f"Remote source manifest is missing the {component} component."
        )
    versions = component_manifest.get("versions")
    if not isinstance(versions, dict) or not versions:
        raise PrepareRuntimeDatabasesError(
            f"Remote source manifest has no versions for {component}."
        )
    version_key = requested_version or component_manifest.get("default_version")
    if not isinstance(version_key, str) or not version_key:
        raise PrepareRuntimeDatabasesError(
            f"Remote source manifest has no default version for {component}."
        )
    version_config = versions.get(version_key)
    if not isinstance(version_config, dict):
        raise PrepareRuntimeDatabasesError(
            f"Remote source manifest does not define version {version_key!r} for {component}."
        )
    return version_key, version_config


def ensure_remote_download_support() -> None:
    """Require aria2 before attempting remote downloads."""
    if shutil.which("aria2c") is None:
        raise PrepareRuntimeDatabasesError(
            "aria2c is required for remote downloads but is not available on PATH."
        )


def infer_download_name(url: str, fallback_name: str) -> str:
    """Infer one local download file name from a URL."""
    parsed = urlparse(url)
    candidate = Path(parsed.path).name
    return candidate or fallback_name


def download_with_aria2(url: str, destination: Path) -> None:
    """Download one remote file with aria2 into a fixed destination path."""
    ensure_parent_directory(destination)
    command = [
        "aria2c",
        "--allow-overwrite=true",
        "--auto-file-renaming=false",
        "--console-log-level=warn",
        "--summary-interval=0",
        "--file-allocation=none",
        "--dir",
        str(destination.parent),
        "--out",
        destination.name,
        url,
    ]
    result = subprocess.run(command, check=False, capture_output=True, text=True)
    if result.returncode != 0:
        raise PrepareRuntimeDatabasesError(
            f"aria2 download failed for {url}: "
            f"{result.stderr.strip() or result.stdout.strip() or 'unknown error'}"
        )


def compute_checksum(path: Path, algorithm: str) -> str:
    """Compute one checksum for a file using the requested algorithm."""
    try:
        digest = hashlib.new(algorithm)
    except ValueError as error:
        raise PrepareRuntimeDatabasesError(
            f"Unsupported checksum algorithm: {algorithm}"
        ) from error
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_checksum_file(path: Path, expected_name: str | None = None) -> str:
    """Parse one checksum file into a plain checksum string."""
    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        parts = stripped.split()
        if expected_name and len(parts) >= 2:
            target_name = parts[-1].lstrip("*")
            if target_name != expected_name:
                continue
        return parts[0]
    raise PrepareRuntimeDatabasesError(f"Could not parse checksum from {path}")


def resolve_expected_checksum(
    checksum_config: dict[str, Any] | None,
    *,
    expected_name: str,
    scratch_dir: Path,
) -> tuple[str, str] | None:
    """Resolve one expected checksum from inline or remote metadata."""
    if checksum_config is None:
        return None
    checksum_type = checksum_config.get("type")
    if not isinstance(checksum_type, str) or not checksum_type:
        raise PrepareRuntimeDatabasesError("Remote checksum configuration is missing its type.")
    inline_value = checksum_config.get("value")
    if isinstance(inline_value, str) and inline_value:
        return checksum_type, inline_value
    checksum_url = checksum_config.get("url")
    if isinstance(checksum_url, str) and checksum_url:
        checksum_name = infer_download_name(checksum_url, f"{expected_name}.{checksum_type}")
        checksum_path = scratch_dir / checksum_name
        download_with_aria2(checksum_url, checksum_path)
        return checksum_type, parse_checksum_file(
            checksum_path,
            checksum_config.get("expected_name") or expected_name,
        )
    return None


def verify_download_checksum(
    path: Path,
    *,
    checksum_config: dict[str, Any] | None,
    scratch_dir: Path,
) -> str:
    """Verify one downloaded file checksum when the manifest provides one."""
    resolved = resolve_expected_checksum(
        checksum_config,
        expected_name=path.name,
        scratch_dir=scratch_dir,
    )
    if resolved is None:
        return "unchecked"
    checksum_type, expected = resolved
    observed = compute_checksum(path, checksum_type)
    if observed.lower() != expected.lower():
        raise ChecksumMismatchError(
            path=path,
            checksum_type=checksum_type,
            expected=expected,
            observed=observed,
        )
    return f"{checksum_type}:{expected.lower()}"


def download_with_checksum_retries(
    *,
    url: str,
    destination: Path,
    checksum_config: dict[str, Any] | None,
    scratch_dir: Path,
) -> str:
    """Download one file and retry transient MD5 checksum mismatches."""
    attempt = 0
    while True:
        attempt += 1
        destination.unlink(missing_ok=True)
        download_with_aria2(url, destination)
        try:
            return verify_download_checksum(
                destination,
                checksum_config=checksum_config,
                scratch_dir=scratch_dir,
            )
        except ChecksumMismatchError as error:
            if error.checksum_type.lower() != "md5" or attempt >= MAX_MD5_DOWNLOAD_ATTEMPTS:
                raise
            LOGGER.warning(
                "Checksum mismatch for %s on attempt %s/%s; retrying download.",
                destination.name,
                attempt,
                MAX_MD5_DOWNLOAD_ATTEMPTS,
            )


def build_remote_scratch_dir(destination: Path, scratch_root: Path | None) -> Path:
    """Create one scratch directory for remote downloads."""
    scratch_base = normalise_path(scratch_root) if scratch_root else destination.parent
    scratch_base.mkdir(parents=True, exist_ok=True)
    scratch_dir = scratch_base / f".download-{destination.name}-{timestamp_token()}"
    if scratch_dir.exists():
        raise PrepareRuntimeDatabasesError(f"Download scratch directory already exists: {scratch_dir}")
    scratch_dir.mkdir(parents=True, exist_ok=False)
    return scratch_dir


def prepare_remote_archive_component(
    *,
    component: str,
    destination: Path,
    validator: Validator,
    scratch_root: Path | None,
    url: str,
    archive_name: str,
    checksum_config: dict[str, Any] | None,
) -> tuple[ValidationResult, str, dict[str, str]]:
    """Download and prepare one remote archive-backed component."""
    scratch_dir = build_remote_scratch_dir(destination, scratch_root)
    try:
        archive_path = scratch_dir / archive_name
        checksum_status = download_with_checksum_retries(
            url=url,
            destination=archive_path,
            checksum_config=checksum_config,
            scratch_dir=scratch_dir,
        )
        validation = prepare_from_archive(
            component=component,
            source=archive_path,
            destination=destination,
            validator=validator,
            scratch_root=scratch_root,
        )
    finally:
        shutil.rmtree(scratch_dir, ignore_errors=True)
    metadata = {
        "source_mode": "remote",
        "transport": "aria2",
        "url": url,
        "checksum": checksum_status,
    }
    return validation, url, metadata


def prepare_remote_component(
    *,
    manifest: dict[str, Any],
    remote_component: str,
    component_label: str,
    version: str | None,
    destination: Path,
    validator: Validator,
    scratch_root: Path | None,
    lineage: str | None = None,
) -> tuple[ValidationResult, str, dict[str, str]]:
    """Prepare one component from the curated remote manifest."""
    ensure_remote_download_support()
    version_key, version_config = resolve_remote_version(
        manifest=manifest,
        component=remote_component,
        requested_version=version,
    )
    kind = version_config.get("kind")
    if not isinstance(kind, str) or not kind:
        raise PrepareRuntimeDatabasesError(
            f"Remote version {version_key!r} for {remote_component} is missing its kind."
        )

    metadata_prefix = {
        "component": component_label,
        "remote_component": remote_component,
        "version": version_key,
    }
    if kind == "archive":
        url = version_config.get("url")
        if not isinstance(url, str) or not url:
            raise PrepareRuntimeDatabasesError(
                f"Remote archive version {version_key!r} for {remote_component} is missing its URL."
            )
        archive_name = version_config.get("archive_name")
        if not isinstance(archive_name, str) or not archive_name:
            archive_name = infer_download_name(url, f"{remote_component}-{version_key}.archive")
        validation, source, metadata = prepare_remote_archive_component(
            component=component_label,
            destination=destination,
            validator=validator,
            scratch_root=scratch_root,
            url=url,
            archive_name=archive_name,
            checksum_config=version_config.get("checksum"),
        )
        metadata.update(metadata_prefix)
        return validation, source, metadata

    if kind == "lineage_archives":
        if not lineage:
            raise PrepareRuntimeDatabasesError(
                f"Remote BUSCO version {version_key!r} requires a lineage name."
            )
        url_template = version_config.get("lineage_url_template")
        if not isinstance(url_template, str) or not url_template:
            raise PrepareRuntimeDatabasesError(
                f"Remote BUSCO version {version_key!r} is missing its lineage URL template."
            )
        archive_name_template = version_config.get("archive_name_template")
        if not isinstance(archive_name_template, str) or not archive_name_template:
            archive_name_template = "{lineage}.tar.gz"
        url = url_template.format(lineage=lineage)
        archive_name = archive_name_template.format(lineage=lineage)
        checksum_config = None
        if "checksums" in version_config:
            checksum_config = version_config["checksums"].get(lineage)
            if not isinstance(checksum_config, dict) or not checksum_config:
                raise PrepareRuntimeDatabasesError(
                    f"Remote BUSCO version {version_key!r} has no checksum for {lineage}."
                )
        validation, source, metadata = prepare_remote_archive_component(
            component=component_label,
            destination=destination,
            validator=validator,
            scratch_root=scratch_root,
            url=url,
            archive_name=archive_name,
            checksum_config=checksum_config,
        )
        metadata.update(metadata_prefix)
        metadata["lineage"] = lineage
        return validation, source, metadata

    raise PrepareRuntimeDatabasesError(
        f"Unsupported remote source kind {kind!r} for {remote_component}."
    )


def resolve_existing_source(path: Path | None) -> Path | None:
    """Resolve one source path only when it exists."""
    if path is None:
        return None
    source = normalise_path(path)
    if source.exists():
        return source
    return None


def prepare_component(
    *,
    component: str,
    remote_component: str,
    source: Path | None,
    destination: Path,
    validator: Validator,
    link_mode: str,
    scratch_root: Path | None,
    force: bool,
    allow_file_source: bool = False,
    download: bool = False,
    version: str | None = None,
    manifest: dict[str, Any] | None = None,
    lineage: str | None = None,
) -> PreparationRecord:
    """Prepare one non-BUSCO database destination."""
    if component in RESOURCE_TOOLS:
        return prepare_annotation_component(
            component=component,
            source=source,
            destination=destination,
            download=download,
            version=version,
            manifest=manifest,
        )
    destination = normalise_path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    original_source = normalise_path(source) if source is not None else None
    local_source = resolve_existing_source(source)

    try:
        _state, existing_validation, marker = assess_destination(
            component=component,
            destination=destination,
            validator=validator,
        )
        if existing_validation is not None and marker is not None:
            LOGGER.info("Reusing ready %s destination: %s", component, destination)
            return PreparationRecord(
                component=component,
                status="present",
                source=str(marker.get("source", "NA")),
                destination=destination,
                details=describe_validation(existing_validation),
            )
    except PrepareRuntimeDatabasesError:
        if not destination.exists() or not force:
            raise
        LOGGER.warning("Replacing existing %s destination at %s", component, destination)
        remove_path(destination)

    if local_source is None and original_source is not None and not download:
        raise PrepareRuntimeDatabasesError(f"Missing {component} source: {original_source}")
    if local_source is None and original_source is not None and download:
        LOGGER.warning(
            "Local %s source does not exist and will be replaced by a remote download: %s",
            component,
            original_source,
        )

    with preparation_lock(destination, force=force):
        if local_source is not None:
            if local_source.is_dir():
                resolved_source = resolve_directory_source(local_source, validator, component)
                materialise_directory(resolved_source, destination, link_mode)
                validation = validator(destination)
            elif allow_file_source and local_source.is_file() and local_source.suffix == ".dmnd":
                materialise_checkm2_file(local_source, destination, link_mode)
                validation = validator(destination)
            else:
                validation = prepare_from_archive(
                    component=component,
                    source=local_source,
                    destination=destination,
                    validator=validator,
                    scratch_root=scratch_root,
                )
            source_description = str(local_source)
            source_metadata = {"source_mode": "local"}
        else:
            if not download:
                raise PrepareRuntimeDatabasesError(
                    f"No local source was supplied for {component}, and remote download is disabled."
                )
            if manifest is None:
                raise PrepareRuntimeDatabasesError(
                    f"Remote download was requested for {component}, but no manifest was loaded."
                )
            validation, source_description, source_metadata = prepare_remote_component(
                manifest=manifest,
                remote_component=remote_component,
                component_label=component,
                version=version,
                destination=destination,
                validator=validator,
                scratch_root=scratch_root,
                lineage=lineage,
            )
        details = describe_validation(validation, extras=source_metadata)
        write_marker(
            component=component,
            source=source_description,
            destination=destination,
            validation=validation,
            source_metadata=source_metadata,
        )
    return PreparationRecord(
        component=component,
        status="prepared",
        source=source_description,
        destination=destination,
        details=details,
    )


def prepare_annotation_component(
    *,
    component: str,
    source: Path | None,
    destination: Path,
    download: bool,
    version: str | None,
    manifest: dict[str, Any] | None,
) -> PreparationRecord:
    """Prepare immutable annotation content through the existing resource API."""
    destination = normalise_path(destination)
    if manifest is None:
        raise PrepareRuntimeDatabasesError("Annotation preparation requires the canonical source manifest")
    version_key, config = resolve_remote_version(
        manifest=manifest, component=component, requested_version=version
    )
    if config.get("kind") != "annotation_bundle":
        raise PrepareRuntimeDatabasesError(f"Unsupported annotation resource source kind for {component}")
    was_present = destination.exists()
    if source is None:
        if was_present:
            source = destination
        elif download:
            source = destination.with_name(destination.name + ".acquired")
            try:
                acquire_component(component, version_key, config, source)
            except (ResourceDownloadError, OSError, ValueError) as error:
                raise PrepareRuntimeDatabasesError(str(error)) from error
        else:
            raise PrepareRuntimeDatabasesError(f"Missing {component} source and download is disabled")
    source = normalise_path(source)
    try:
        resource = prepare_resource(component, version_key, source, destination)
    except (AnnotationResourceError, OSError, ValueError) as error:
        raise PrepareRuntimeDatabasesError(str(error)) from error
    validation = ValidationResult(
        required_paths=tuple(entry["path"] for entry in resource["contract"]["files"]),
        details={"resource_id": resource["resource_id"], "version": version_key},
    )
    write_marker(component=component, source=str(source), destination=destination, validation=validation)
    return PreparationRecord(component, "present" if was_present else "prepared", str(source), destination, describe_validation(validation))


def determine_busco_lineages(mapping: dict[str, Path]) -> tuple[str, ...]:
    """Determine the BUSCO lineages that should exist after preparation."""
    if mapping:
        return tuple(sorted(mapping))
    return DEFAULT_BUSCO_LINEAGES


def prepare_busco_lineages(
    *,
    source_mapping: dict[str, Path],
    destination_root: Path,
    link_mode: str,
    scratch_root: Path | None,
    force: bool,
    download: bool,
    version: str | None,
    manifest: dict[str, Any] | None,
) -> list[PreparationRecord]:
    """Prepare the requested BUSCO lineage directories beneath one root."""
    destination_root = normalise_path(destination_root)
    destination_root.mkdir(parents=True, exist_ok=True)
    lineages = determine_busco_lineages(source_mapping)
    records: list[PreparationRecord] = []
    for lineage in lineages:
        source = source_mapping.get(lineage)
        destination = destination_root / lineage
        validator = build_busco_lineage_validator(lineage)
        component = f"busco:{lineage}"
        records.append(
            prepare_component(
                component=component,
                remote_component="busco",
                source=source,
                destination=destination,
                validator=validator,
                link_mode=link_mode,
                scratch_root=scratch_root,
                force=force,
                download=download,
                version=version,
                manifest=manifest,
                lineage=lineage,
            )
        )

    root_validation = build_busco_validator(lineages)(destination_root)
    records.append(
        PreparationRecord(
            component="busco_root",
            status="prepared",
            source="derived",
            destination=destination_root,
            details=describe_validation(root_validation),
        )
    )
    return records


def build_component_records(args: argparse.Namespace) -> list[PreparationRecord]:
    """Prepare the requested database destinations and return report records."""
    annotation_requested = any(getattr(args, component + "_dest") for component in RESOURCE_TOOLS)
    manifest = load_remote_source_manifest(args.remote_source_manifest) if args.download or annotation_requested else None

    if args.taxdump_source and not args.taxdump_dest:
        raise PrepareRuntimeDatabasesError("Taxdump source requires --taxdump-dest.")
    if args.taxdump_version and not args.taxdump_dest:
        raise PrepareRuntimeDatabasesError("Taxdump version requires --taxdump-dest.")
    if args.checkm2_source and not args.checkm2_dest:
        raise PrepareRuntimeDatabasesError("CheckM2 source requires --checkm2-dest.")
    if args.checkm2_version and not args.checkm2_dest:
        raise PrepareRuntimeDatabasesError("CheckM2 version requires --checkm2-dest.")
    if args.codetta_source and not args.codetta_dest:
        raise PrepareRuntimeDatabasesError("Codetta source requires --codetta-dest.")
    if args.codetta_version and not args.codetta_dest:
        raise PrepareRuntimeDatabasesError("Codetta version requires --codetta-dest.")
    if args.eggnog_source and not args.eggnog_dest:
        raise PrepareRuntimeDatabasesError("eggNOG source requires --eggnog-dest.")
    if args.eggnog_version and not args.eggnog_dest:
        raise PrepareRuntimeDatabasesError("eggNOG version requires --eggnog-dest.")
    for component in ("cogclassifier", "pfam", "kofam"):
        if (getattr(args, component + "_source") or getattr(args, component + "_version")) and not getattr(args, component + "_dest"):
            raise PrepareRuntimeDatabasesError(f"{component} source/version requires --{component}-dest.")
    if args.padloc_source and not args.padloc_dest:
        raise PrepareRuntimeDatabasesError("PADLOC source requires --padloc-dest.")
    if args.padloc_version and not args.padloc_dest:
        raise PrepareRuntimeDatabasesError("PADLOC version requires --padloc-dest.")
    if args.busco_dest_root is None and (args.busco_lineage_source or args.busco_version):
        raise PrepareRuntimeDatabasesError(
            "BUSCO lineage sources or a BUSCO version require --busco-dest-root."
        )

    busco_sources = {
        lineage: path
        for lineage, path in (
            parse_name_path(token, "--busco-lineage-source")
            for token in args.busco_lineage_source
        )
    }

    records: list[PreparationRecord] = []
    if args.taxdump_dest:
        records.append(
            prepare_component(
                component="taxdump",
                remote_component="taxdump",
                source=args.taxdump_source,
                destination=args.taxdump_dest,
                validator=validate_taxdump,
                link_mode=args.link_mode,
                scratch_root=args.scratch_root,
                force=args.force,
                download=args.download,
                version=args.taxdump_version,
                manifest=manifest,
            )
        )
    if args.checkm2_dest:
        records.append(
            prepare_component(
                component="checkm2",
                remote_component="checkm2",
                source=args.checkm2_source,
                destination=args.checkm2_dest,
                validator=validate_checkm2,
                link_mode=args.link_mode,
                scratch_root=args.scratch_root,
                force=args.force,
                allow_file_source=True,
                download=args.download,
                version=args.checkm2_version,
                manifest=manifest,
            )
        )
    if args.busco_dest_root:
        records.extend(
            prepare_busco_lineages(
                source_mapping=busco_sources,
                destination_root=args.busco_dest_root,
                link_mode=args.link_mode,
                scratch_root=args.scratch_root,
                force=args.force,
                download=args.download,
                version=args.busco_version,
                manifest=manifest,
            )
        )
    if args.codetta_dest:
        records.append(
            prepare_component(
                component="codetta",
                remote_component="codetta",
                source=args.codetta_source,
                destination=args.codetta_dest,
                validator=validate_codetta,
                link_mode=args.link_mode,
                scratch_root=args.scratch_root,
                force=args.force,
                download=args.download,
                version=args.codetta_version,
                manifest=manifest,
            )
        )
    if args.eggnog_dest:
        records.append(
            prepare_component(
                component="eggnog",
                remote_component="eggnog",
                source=args.eggnog_source,
                destination=args.eggnog_dest,
                validator=validate_eggnog,
                link_mode=args.link_mode,
                scratch_root=args.scratch_root,
                force=args.force,
                download=args.download,
                version=args.eggnog_version,
                manifest=manifest,
            )
        )
    for component in ("cogclassifier", "pfam", "kofam"):
        if getattr(args, component + "_dest"):
            records.append(prepare_annotation_component(
                component=component,
                source=getattr(args, component + "_source"),
                destination=getattr(args, component + "_dest"),
                download=args.download,
                version=getattr(args, component + "_version"),
                manifest=manifest,
            ))
    if args.padloc_dest:
        records.append(
            prepare_component(
                component="padloc",
                remote_component="padloc",
                source=args.padloc_source,
                destination=args.padloc_dest,
                validator=validate_padloc,
                link_mode=args.link_mode,
                scratch_root=args.scratch_root,
                force=args.force,
                download=args.download,
                version=args.padloc_version,
                manifest=manifest,
            )
        )
    if not records:
        raise PrepareRuntimeDatabasesError(
            "No database targets were requested. Supply at least one destination."
        )
    return records


def write_report(path: Path, records: Sequence[PreparationRecord]) -> None:
    """Write one TSV report describing the prepared destinations."""
    report_path = normalise_path(path)
    ensure_parent_directory(report_path)
    with report_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for record in records:
            writer.writerow(report_row(record))


def build_nextflow_arguments(records: Sequence[PreparationRecord]) -> list[tuple[str, str]]:
    """Build one copy-pastable Nextflow argument block."""
    mapping: dict[str, str] = {}
    for record in records:
        if record.component == "taxdump":
            mapping["--taxdump"] = str(record.destination)
        elif record.component == "checkm2":
            mapping["--checkm2_db"] = str(record.destination)
        elif record.component == "codetta":
            mapping["--codetta_db"] = str(record.destination)
        elif record.component == "busco_root":
            mapping["--busco_db"] = str(record.destination)
        elif record.component == "eggnog":
            mapping["--eggnog_db"] = str(record.destination)
        elif record.component in ("cogclassifier", "pfam", "kofam"):
            mapping[f"--{record.component}_db"] = str(record.destination)
    return sorted(mapping.items())


def print_nextflow_arguments(records: Sequence[PreparationRecord]) -> None:
    """Print one ready-to-paste Nextflow command fragment."""
    arguments = build_nextflow_arguments(records)
    if not arguments:
        return
    print("Nextflow arguments:")
    print("nextflow run . \\")
    for index, (flag, value) in enumerate(arguments):
        suffix = " \\" if index < len(arguments) - 1 else ""
        print(f"  {flag} {shlex.quote(value)}{suffix}")


def main(argv: Sequence[str] | None = None) -> int:
    """Prepare requested runtime databases and print reusable paths."""
    args = parse_args(argv)
    configure_logging()
    try:
        records = build_component_records(args)
        if args.report:
            write_report(args.report, records)
        print_nextflow_arguments(records)
    except PrepareRuntimeDatabasesError as error:
        LOGGER.error("%s", error)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
