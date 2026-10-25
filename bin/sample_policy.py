"""Explicit domain and QC publication contracts; no taxonomic inference."""

from __future__ import annotations

import re
from collections.abc import Mapping
from typing import Any

DOMAINS = {"Bacteria": "bac", "Archaea": "arc"}
DOMAIN_FIELDS = (
    "annotation_domain",
    "domain_source",
    "prokka_kingdom",
    "barrnap_kingdom",
    "busco_lineages",
    "primary_busco_lineage",
)
PROVENANCE_COLUMNS = (
    "Annotation_Domain",
    "Annotation_Domain_Source",
    "Primary_BUSCO_Lineage",
    "Primary_BUSCO_Score",
    "QC_Scope",
)


def domain(value: Any) -> str:
    """Accept only explicit bacterial or archaeal domain names."""
    if isinstance(value, str):
        for name in DOMAINS:
            if value.strip().casefold() == name.casefold():
                return name
    raise ValueError(f"Domain must be Bacteria or Archaea; received {value!r}")


def profiles(value: Any) -> dict[str, list[str]]:
    """Validate operator-selected lineage profiles without inferring subclades."""
    if not isinstance(value, dict):
        raise ValueError(  # noqa: TRY004 - invalid configuration is a ValueError contract
            "busco_domain_lineages must be a map of domains to lineage lists"
        )
    result: dict[str, list[str]] = {}
    used: set[str] = set()
    for raw, lineages in value.items():
        name = domain(raw)
        if name in result or not isinstance(lineages, list) or not lineages:
            raise ValueError(
                "Domain profiles require unique domains and nonempty lineage lists"
            )
        checked = []
        for lineage in lineages:
            if not isinstance(lineage, str) or not re.fullmatch(
                r"[a-z][a-z0-9_]*(?:\.[0-9]+)?", lineage
            ):
                raise ValueError(f"Invalid BUSCO lineage: {lineage!r}")
            if lineage in used:
                raise ValueError(
                    f"BUSCO lineage occurs in multiple profile entries: {lineage}"
                )
            opposite = "archaea" if name == "Bacteria" else "bacteria"
            if lineage.startswith(opposite + "_odb"):
                raise ValueError(
                    f"BUSCO lineage {lineage} conflicts with domain {name}"
                )
            used.add(lineage)
            checked.append(lineage)
        result[name] = checked
    return result


def resolve_sample(
    metadata: Mapping[str, str], settings: Mapping[str, Any]
) -> dict[str, str]:
    """Resolve one declared domain and derive all tool-specific settings."""
    removed = {"domain_column", "prokka_kingdom", "barrnap_kingdom"} & settings.keys()
    if removed:
        raise ValueError(
            f"Unsupported options: {', '.join(sorted(removed))}; use --domain"
        )
    mode = settings.get("domain")
    if mode not in ("bacteria", "archaea", "metadata"):
        raise ValueError(
            "--domain is required and must be bacteria, archaea or metadata"
        )
    configured = profiles(settings.get("busco_domain_lineages", {}))
    declared = metadata.get("Domain")
    name = domain(declared) if mode == "metadata" else domain(mode)
    if (
        mode != "metadata"
        and declared is not None
        and declared.strip()
        and domain(declared) != name
    ):
        raise ValueError("Domain metadata conflicts with --domain")
    # Explicit bacteria without profiles keeps the documented v0.4.3 contract.
    if mode == "bacteria" and not configured:
        return {}
    if name not in configured:
        raise ValueError(f"Missing explicit BUSCO lineage profile for {name}")
    if settings.get("busco_primary_column") is not None:
        raise ValueError(
            "A cohort-wide busco_primary_column cannot be used with domain profiles"
        )
    lineages = configured[name]
    return {
        "annotation_domain": name,
        "domain_source": "metadata:Domain" if mode == "metadata" else "domain",
        "prokka_kingdom": name,
        "barrnap_kingdom": DOMAINS[name],
        "busco_lineages": ",".join(lineages),
        "primary_busco_lineage": lineages[0],
    }


def sample_lineages(
    sample: Mapping[str, str], legacy: list[str] | tuple[str, ...]
) -> list[str]:
    """Read and validate the complete per-sample domain record, or the legacy contract."""
    present = set(DOMAIN_FIELDS) & sample.keys()
    if not present:
        return list(legacy)
    if any(
        not isinstance(sample[field], str) or not sample[field] for field in present
    ):
        raise ValueError("Invalid per-sample domain policy field")
    if present != set(DOMAIN_FIELDS):
        raise ValueError("Incomplete per-sample domain policy")
    if sample["domain_source"] not in ("domain", "metadata:Domain"):
        raise ValueError("Missing or unsupported domain declaration provenance")
    name = domain(sample["annotation_domain"])
    if sample["prokka_kingdom"] != name or sample["barrnap_kingdom"] != DOMAINS[name]:
        raise ValueError("Conflicting per-sample domain and kingdom settings")
    lineages = profiles({name: sample["busco_lineages"].split(",")})[name]
    if sample["primary_busco_lineage"] != lineages[0]:
        raise ValueError("Primary BUSCO lineage differs from the declared profile")
    return lineages


def primary_column(sample: Mapping[str, str], legacy: str | None) -> str | None:
    """Select the actually run lineage, preserving the original ANI score formula."""
    lineages = sample_lineages(sample, [])
    return "BUSCO_" + lineages[0] if lineages else legacy


def qc_scope(sample: Mapping[str, str]) -> str:
    """Absent scope denotes the documented complete v0.4 publication contract."""
    value = sample.get("qc_scope", "full")
    if value not in ("full", "qc"):
        raise ValueError(f"Unknown published QC scope: {value!r}")
    return value
