"""Read the checked native Prodigal outcome without treating absence as zero."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from annotation_common import AnnotationError

PREFIX = "NF_PROKKA_PREDICTION_V1\t"


@dataclass(frozen=True)
class Prediction:
    mode: str
    cds_count: int


def parse_prediction(mode: str, count: str) -> Prediction:
    if mode not in ("single", "meta") or not re.fullmatch(r"0|[1-9][0-9]*", count):
        raise AnnotationError("Invalid native Prodigal prediction record")
    return Prediction(mode, int(count))


def read_prediction(path: Path) -> Prediction:
    lines = path.read_text().splitlines()
    exits = [
        line.removeprefix("exit_code=")
        for line in lines
        if line.startswith("exit_code=")
    ]
    if exits != ["0"]:
        raise AnnotationError("upstream_failed: Prokka did not complete successfully")
    # Earlier failed Prokka attempts may have completed their Prodigal child.
    # Only the final attempt can establish the published prediction outcome.
    starts = [index for index, line in enumerate(lines) if line.startswith("attempt=")]
    if starts:
        lines = lines[starts[-1] :]
    records = [
        line.removeprefix(PREFIX).split("\t")
        for line in lines
        if line.startswith(PREFIX)
    ]
    if len(records) != 1 or len(records[0]) != 2:
        raise AnnotationError(
            "upstream_failed: Missing or ambiguous native Prodigal prediction record"
        )
    try:
        return parse_prediction(*records[0])
    except AnnotationError as error:
        raise AnnotationError(f"upstream_failed: {error}") from error


def status_from_manifest(row: dict[str, str]) -> tuple[str, list[str]]:
    if row.get("exit_code") != "0":
        return "failed", ["prokka_failed"]
    try:
        prediction = parse_prediction(
            row.get("prodigal_mode", "NA"), row.get("cds_count", "NA")
        )
    except AnnotationError:
        return "failed", ["invalid_prokka_prediction"]
    try:
        gff_size, faa_size = int(row["gff_size"]), int(row["faa_size"])
    except (KeyError, ValueError):
        return "failed", ["missing_prokka_outputs"]
    if gff_size <= 0 or faa_size < 0:
        return "failed", ["missing_prokka_outputs"]
    if prediction.cds_count == 0 and faa_size == 0:
        return "no_cds", []
    if prediction.cds_count > 0 and faa_size > 0:
        return "done", []
    return "failed", ["inconsistent_prokka_outputs"]
