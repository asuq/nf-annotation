"""Shared evidence bookkeeping for the five pinned annotation adapters."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from annotation_common import AnnotationError, write_tsv


@dataclass
class Normalized:
    """Validated evidence, accepted gene/feature sets, and optional-field errors."""

    tables: dict[str, tuple[tuple[str, ...], list[dict[str, Any]]]] = field(
        default_factory=dict
    )
    features: dict[str, dict[str, set[str]]] = field(default_factory=dict)
    definitions: dict[str, dict[str, str]] = field(default_factory=dict)
    mapped: set[str] = field(default_factory=set)
    accepted: set[str] = field(default_factory=set)
    errors: list[dict[str, str]] = field(default_factory=list)
    reported_hits: int = 0

    def add(self, name: str, gene: str, values: list[str] | set[str]) -> None:
        """Accumulate distinct features without discarding repeated evidence rows."""
        if values:
            self.features.setdefault(name, {}).setdefault(gene, set()).update(values)
            self.accepted.add(gene)

    def error(self, gene: str, name: str, raw: str, reason: str) -> None:
        """Flag an optional field; aggregation must invalidate its whole genome cell."""
        self.errors.append(dict(gene_id=gene, field=name, raw_value=raw, reason=reason))

    def publish(self, root: Path) -> dict[str, Any]:
        """Write declared evidence tables and return JSON-compatible summary data."""
        root.mkdir(exist_ok=True)
        for name, (columns, rows) in self.tables.items():
            write_tsv(root / name, columns, rows)
        return {
            "features": {
                name: {gene: sorted(values) for gene, values in sorted(genes.items())}
                for name, genes in sorted(self.features.items())
            },
            "definitions": self.definitions,
            "mapped_genes": sorted(self.mapped),
            "accepted_genes": sorted(self.accepted),
            "field_errors": self.errors,
            "reported_hits": self.reported_hits,
        }


def json_cell(value: Any) -> str:
    """Encode structured evidence into one lossless TSV cell."""
    return json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(",", ":"))


def native_rows(
    path: Path, columns: tuple[str, ...], *, delimiter: str = "\t"
) -> list[dict[str, str]]:
    """Read an exact native header and reject truncated or over-wide rows."""
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        if next(reader, None) != list(columns):
            raise AnnotationError(f"Unexpected native schema: {path}")
        result = []
        for line, row in enumerate(reader, 2):
            if len(row) != len(columns):
                raise AnnotationError(f"Malformed native row: {path}:{line}")
            result.append(dict(zip(columns, row, strict=True)))
        return result


def ordered_categories(value: str, vocabulary: set[str]) -> list[str]:
    """Validate native ordered category letters without trimming malformed values."""
    if not value or any(letter not in vocabulary for letter in value):
        raise AnnotationError("category contains an empty or unknown value")
    if len(value) != len(set(value)):
        raise AnnotationError("category contains duplicate letters")
    return list(value)
