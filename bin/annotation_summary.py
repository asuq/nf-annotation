"""Shared count definitions for annotation master fields and feature matrices."""

from __future__ import annotations

from collections import Counter
from decimal import Decimal, localcontext
from fractions import Fraction
from typing import Any

from annotation_common import TOOLS

MATRICES = {
    "eggnog_ko_counts": ("eggnog", "ko"),
    "eggnog_go_counts": ("eggnog", "go"),
    "eggnog_ec_counts": ("eggnog", "ec"),
    "kofam_ko_counts": ("kofam", "ko"),
    "cog_counts": ("cogclassifier", "cog"),
    "pfam_gene_counts": ("pfam", "pfam"),
}
FIELDS = {
    "eggnog": (
        "Preferred_name",
        "ko",
        "go",
        "ec",
        "KEGG_Pathway",
        "KEGG_Module",
        "KEGG_Reaction",
        "KEGG_rclass",
        "BRITE",
        "KEGG_TC",
        "CAZy",
        "BiGG_Reaction",
        "pfam",
    ),
    "cogclassifier": ("cog", "categories"),
    "pfam": ("pfam",),
    "kofam": ("ko",),
    "padloc": ("systems",),
}
COG_CATEGORIES = tuple(sorted("JAKLBDYVTMNZWUOCEFGHIPQRSX"))
TOOL_METRICS = (
    "status",
    "analysed_proteins",
    "mapped_proteins",
    "accepted_proteins",
    "mapped_fraction",
    "accepted_fraction",
)
ANNOTATION_COLUMNS = (
    "annotation_input_proteins",
    "annotation_complete",
    *[f"{tool}_{metric}" for tool in TOOLS for metric in TOOL_METRICS],
    *[
        f"{tool}_{field}_{metric}"
        for tool, field in MATRICES.values()
        for metric in ("proteins", "assignments")
    ],
    *[
        f"cog_category_{letter}_{metric}"
        for letter in COG_CATEGORIES
        for metric in ("gene_count", "fractional_gene_count")
    ],
    "padloc_systems",
    "annotation_field_errors",
)


def decimal_fraction(value: Fraction) -> str:
    """Serialize exact internal counts to 12 decimal places, without binary drift."""
    with localcontext() as context:
        context.prec = 40
        number = (Decimal(value.numerator) / Decimal(value.denominator)).quantize(
            Decimal("0.000000000001")
        )
    return format(number, "f").rstrip("0").rstrip(".") or "0"


def bad_fields(result: dict[str, Any]) -> set[str]:
    """Return every field whose whole-genome summary must remain unavailable."""
    return {
        error["field"]
        for error in (result.get("evidence") or {}).get("field_errors", [])
    }


def accepted_count(result: dict[str, Any]) -> int | None:
    """Do not turn a malformed functional field into a complete coverage count."""
    if result["status"] != "success" or bad_fields(result) & set(
        FIELDS[result["tool"]]
    ) - {"categories"}:
        return None
    return len(result["evidence"]["accepted_genes"])


def feature_counts(result: dict[str, Any], field: str) -> Counter[str] | None:
    """Count distinct genes per feature, or return unknown for an invalid field."""
    if result["status"] != "success" or field in bad_fields(result):
        return None
    genes = result["evidence"]["features"].get(field, {})
    return Counter(feature for values in genes.values() for feature in set(values))


def genome_summary(
    results: dict[str, dict[str, Any]], input_count: int | None, complete: bool
) -> dict[str, Any]:
    """Derive the master row using the same gene sets consumed by matrices."""
    row: dict[str, Any] = dict.fromkeys(ANNOTATION_COLUMNS)
    row.update(
        annotation_input_proteins=input_count, annotation_complete=str(complete).lower()
    )
    errors = []
    for tool, result in results.items():
        row[f"{tool}_status"] = result["status"]
        if result["status"] != "success":
            continue
        evidence = result["evidence"]
        mapped, accepted = len(evidence["mapped_genes"]), accepted_count(result)
        row.update(
            {
                f"{tool}_analysed_proteins": input_count,
                f"{tool}_mapped_proteins": mapped,
                f"{tool}_accepted_proteins": accepted,
                f"{tool}_mapped_fraction": decimal_fraction(
                    Fraction(mapped, input_count)
                )
                if input_count
                else None,
                f"{tool}_accepted_fraction": decimal_fraction(
                    Fraction(accepted, input_count)
                )
                if input_count and accepted is not None
                else None,
            }
        )
        errors.extend(f"{tool}:{field}" for field in sorted(bad_fields(result)))
    for tool, field in MATRICES.values():
        result = results[tool]
        counts = feature_counts(result, field)
        if counts is not None:
            row[f"{tool}_{field}_assignments"] = sum(counts.values())
            row[f"{tool}_{field}_proteins"] = len(
                result["evidence"]["features"].get(field, {})
            )
    cog = results["cogclassifier"]
    counts = feature_counts(cog, "categories")
    if counts is not None:
        fractional: dict[str, Fraction] = {
            letter: Fraction(0) for letter in COG_CATEGORIES
        }
        for categories in cog["evidence"]["features"].get("categories", {}).values():
            for letter in set(categories):
                fractional[letter] += Fraction(1, len(set(categories)))
        for letter in COG_CATEGORIES:
            row[f"cog_category_{letter}_gene_count"] = counts[letter]
            row[f"cog_category_{letter}_fractional_gene_count"] = decimal_fraction(
                fractional[letter]
            )
    systems = feature_counts(results["padloc"], "systems")
    if systems is not None:
        row["padloc_systems"] = len(systems)
    row["annotation_field_errors"] = ";".join(errors) or "none"
    return row
