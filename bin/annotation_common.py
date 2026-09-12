"""Shared strict I/O and identities for the v0.4 annotation contract."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections.abc import Iterable
from pathlib import Path
from typing import Any

SCHEMA_VERSION = 1
TOOLS = ("eggnog", "cogclassifier", "pfam", "kofam", "padloc")
STATES = frozenset(
    {
        "success",
        "failed",
        "skipped_disabled",
        "skipped_not_selected",
        "upstream_failed",
        "incompatible_input",
    }
)
STATUS_COLUMNS = (
    "accession",
    "tool",
    "status",
    "reason",
    "exit_code",
    "input_proteins",
    "reported_hits",
    "accepted_proteins",
    "search_fingerprint",
    "normalization_fingerprint",
    "method_id",
    "schema_version",
    "field_errors",
)
PROTEIN_COLUMNS = (
    "accession",
    "gene_id",
    "protein_id",
    "tool_id",
    "description",
    "length",
    "sha256",
    "genetic_code",
    "contig_id",
    "source_contig_id",
    "contig_length",
    "strand",
    "segment_count",
    "topology",
    "genbank_location",
    "genbank_qualifiers",
)
COORDINATE_COLUMNS = (
    "accession",
    "gene_id",
    "feature_id",
    "contig_id",
    "source_contig_id",
    "segment",
    "start",
    "end",
    "strand",
    "phase",
    "attributes",
)


class AnnotationError(ValueError):
    """An input, resource or tool output violates an annotation invariant."""


def validate_accession(value: Any) -> None:
    """Require an accession that can identify one published sample directory."""
    if (
        not isinstance(value, str)
        or not value
        or value in (".", "..")
        or any(char in value for char in ("/", "\\"))
        or any(ord(char) < 32 for char in value)
    ):
        raise AnnotationError("Invalid published accession")


def digest(path: Path) -> str:
    """Hash file bytes with SHA-256."""
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def identity(value: Any) -> str:
    """Hash a canonical JSON representation of an analysis identity."""
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def write_json(path: Path, value: Any) -> None:
    """Write deterministic ASCII JSON with a final newline."""
    path.write_text(
        json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True) + "\n"
    )


def read_json(path: Path) -> Any:
    """Read JSON and report malformed content as an annotation error."""
    try:
        return json.loads(path.read_text())
    except (json.JSONDecodeError, UnicodeError) as exc:
        raise AnnotationError(f"Invalid JSON in {path}: {exc}") from exc


def write_tsv(
    path: Path, columns: Iterable[str], rows: Iterable[dict[str, Any]]
) -> None:
    """Write the declared TSV columns, representing unavailable values as NA."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {key: "NA" if value is None else value for key, value in row.items()}
            )


def read_tsv(path: Path, required: Iterable[str] = ()) -> list[dict[str, str]]:
    """Read a TSV with unique column names and complete rows."""
    # Detailed per-protein validation flags can exceed csv's default field
    # limit. A field cannot exceed its source file's byte length. Restore the
    # process-wide parser setting after this synchronous read, including errors.
    previous_limit = csv.field_size_limit()
    csv.field_size_limit(max(previous_limit, path.stat().st_size))
    try:
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = reader.fieldnames
            if (
                not header
                or any(not name for name in header)
                or len(header) != len(set(header))
                or not set(required) <= set(header)
            ):
                raise AnnotationError(
                    f"Invalid table header in {path}; required: {list(required)}"
                )
            rows = []
            for line, row in enumerate(reader, 2):
                if None in row or None in row.values():
                    raise AnnotationError(f"Malformed table row in {path}:{line}")
                rows.append(row)
            return rows
    finally:
        csv.field_size_limit(previous_limit)


def number(value: str, field: str, *, minimum: float | None = None) -> float:
    """Parse a finite numeric value and enforce its lower bound when specified."""
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise AnnotationError(f"Invalid {field}: {value!r}") from exc
    if not math.isfinite(result) or (minimum is not None and result < minimum):
        raise AnnotationError(f"Invalid {field}: {value!r}")
    return result


def integer(value: str, field: str, *, minimum: int = 0) -> int:
    """Parse an integer without rounding and enforce its lower bound."""
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise AnnotationError(f"Invalid {field}: {value!r}") from exc
    if result < minimum:
        raise AnnotationError(f"Invalid {field}: {value!r}")
    return result


def selected_tools(value: str) -> tuple[str, ...]:
    """Validate the explicit tool list without silently discarding unknown names."""
    selected = tuple(value.split(",")) if value else ()
    if len(selected) != len(set(selected)) or any(
        tool not in TOOLS for tool in selected
    ):
        raise AnnotationError(
            f"annotation_tools must be a unique comma-separated subset of {','.join(TOOLS)}"
        )
    return selected


def bundle_proteins(bundle: Path) -> tuple[dict[str, Any], list[dict[str, str]]]:
    """Verify all bundle file checksums and identities before exposing proteins."""
    manifest = read_json(bundle / "bundle.json")
    if (
        not isinstance(manifest, dict)
        or manifest.get("schema_version") != SCHEMA_VERSION
        or manifest.get("status") != "success"
    ):
        raise AnnotationError(
            f"Unsupported or unsuccessful annotation bundle: {bundle}"
        )
    required = {
        "proteins.faa",
        "protein_manifest.tsv",
        "gene_coordinates.tsv",
        "source.faa",
        "source.gff",
        "source.gbk",
        "genome.fasta",
    }
    files = manifest.get("files")
    if not isinstance(files, dict) or set(files) != required:
        raise AnnotationError(
            "Bundle must record all required source and derived files"
        )
    for name, expected in files.items():
        path = bundle / name
        if Path(name).name != name or not path.is_file() or digest(path) != expected:
            raise AnnotationError(f"Bundle content mismatch: {name}")
    rows = read_tsv(bundle / "protein_manifest.tsv", PROTEIN_COLUMNS)
    if not rows or len(rows) != manifest.get("protein_count"):
        raise AnnotationError("Bundle protein count mismatch")
    for key in ("gene_id", "protein_id", "tool_id"):
        if len({row[key] for row in rows}) != len(rows):
            raise AnnotationError(f"Duplicate bundle {key}")
    if any(
        row["accession"] != manifest.get("accession")
        or row["genetic_code"] != str(manifest.get("genetic_code"))
        for row in rows
    ):
        raise AnnotationError("Bundle accession or genetic-code mismatch")
    expected_input = identity(
        {
            "accession": manifest["accession"],
            "genetic_code": manifest["genetic_code"],
            "proteins": [
                {key: row[key] for key in ("gene_id", "tool_id", "sha256")}
                for row in rows
            ],
        }
    )
    if manifest.get("input_id") != expected_input:
        raise AnnotationError("Bundle protein identity mismatch")
    if manifest.get("coordinate_id") != files["gene_coordinates.tsv"]:
        raise AnnotationError("Bundle coordinate identity mismatch")
    if manifest.get("source_genome_sha256") != files["genome.fasta"]:
        raise AnnotationError("Bundle source genome identity mismatch")
    return manifest, rows


def protein_lookup(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    """Index original and compact IDs while rejecting cross-namespace ambiguity."""
    result = {}
    for row in rows:
        for key in (row["tool_id"], row["protein_id"]):
            if key in result and result[key]["gene_id"] != row["gene_id"]:
                raise AnnotationError(f"Ambiguous protein identifier {key!r}")
            result[key] = row
    return result


def query(lookup: dict[str, dict[str, str]], value: str) -> dict[str, str]:
    """Resolve a tool query ID to exactly one declared input protein."""
    if value not in lookup:
        raise AnnotationError(f"Tool returned unknown protein ID {value!r}")
    return lookup[value]
