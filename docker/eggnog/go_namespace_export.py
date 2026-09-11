"""Export native eggNOG namespace winners without changing native annotations."""

from collections.abc import Iterable, Mapping
from typing import Any, TextIO

KEY = "_nf_go_namespaces"
FIELDS = ("gos_mf", "gos_bp", "gos_cc")
HEADER = [
    "query",
    *[item for field in FIELDS for item in (field, f"{field}_confidence")],
]


def record_namespaces(
    annotations: dict[str, Any],
    values: Mapping[str, Iterable[str]],
    tiers: Mapping[str, int],
    confidence_labels: Mapping[int, str],
) -> None:
    """Copy the engine's existing selections; leave its merged fields untouched."""
    if values:
        annotations[KEY] = {
            field: {
                "terms": sorted(values[field]),
                "confidence": confidence_labels[tiers[field]],
            }
            for field in FIELDS
            if field in values
        }


def write_row(out: TextIO, annotation: tuple) -> None:
    """Write one current-format annotation, checking namespace/merged agreement."""
    if len(annotation) != 14:
        raise ValueError("GO export requires the pinned eggNOG v3 annotation schema")
    query, annotations = annotation[0], annotation[4]
    namespaces = annotations.get(KEY, {})
    merged = set()
    row = [query]
    for field in FIELDS:
        namespace = namespaces.get(field)
        if namespace is None:
            row.extend(("-", "-"))
        else:
            terms, confidence = namespace["terms"], namespace["confidence"]
            if not terms or confidence not in {"high", "medium", "low"}:
                raise ValueError(
                    f"Invalid native GO namespace evidence for {query}/{field}"
                )
            merged.update(terms)
            row.extend((",".join(terms), confidence))
    if merged != set(annotations.get("GOs", ())):
        raise ValueError(f"GO namespace evidence disagrees with native GOs for {query}")
    if any(
        not isinstance(value, str) or any(c in value for c in "\t\r\n") for value in row
    ):
        raise ValueError("GO namespace export contains an invalid TSV value")
    print("\t".join(row), file=out)
