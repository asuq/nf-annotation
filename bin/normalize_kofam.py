"""Normalize native KOfamScan candidates and its adaptive-threshold pass marker."""

from __future__ import annotations

import csv
import re
from pathlib import Path

from annotation_common import AnnotationError, number, protein_lookup, query, read_tsv
from annotation_normalization import Normalized

COLUMNS = (
    "accession",
    "gene_id",
    "protein_id",
    "ko",
    "native_marker",
    "threshold",
    "printed_threshold",
    "score_type",
    "score",
    "evalue",
    "definition",
    "accepted",
    "state",
)


def normalize(raw: Path, proteins: list[dict[str, str]], resource: Path) -> Normalized:
    """Preserve all reported candidates and trust native unrounded threshold decisions."""
    result, lookup = Normalized(), protein_lookup(proteins)
    metadata = read_tsv(
        resource / "ko_list", ("knum", "threshold", "score_type", "definition")
    )
    kos = {row["knum"]: row for row in metadata}
    if len(kos) != len(metadata):
        raise AnnotationError("Duplicate KOfam KO metadata")
    profiles = read_tsv(
        resource / "prokaryote_profiles.tsv",
        ("ko", "profile", "threshold", "score_type"),
    )
    selected = {row["ko"] for row in profiles}
    evidence, seen = [], set()
    with (raw / "kofam.tsv").open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if next(reader, None) != [
            "#",
            "gene name",
            "KO",
            "thrshld",
            "score",
            "E-value",
            "KO definition",
        ]:
            raise AnnotationError("Unexpected KOfamScan detail-tsv header")
        if next(reader, None) != [
            "#",
            "---------",
            "------",
            "-------",
            "------",
            "---------",
            "-------------",
        ]:
            raise AnnotationError("Missing KOfamScan detail-tsv header delimiter")
        for line, row in enumerate(reader, 3):
            if len(row) != 7:
                raise AnnotationError(f"Malformed KOfamScan row {line}")
            mark, identifier, ko, printed_threshold, score, evalue, definition = row
            protein = query(lookup, identifier)
            gene = protein["gene_id"]
            if (
                mark not in ("", "*")
                or not re.fullmatch(r"K[0-9]{5}", ko)
                or ko not in selected
                or ko not in kos
                or (gene, ko) in seen
            ):
                raise AnnotationError(
                    "Invalid, unselected or duplicate KOfam candidate"
                )
            seen.add((gene, ko))
            data = kos[ko]
            number(score, "KOfam score")
            number(evalue, "KOfam E-value", minimum=0)
            available = data["threshold"] != "-"
            if available:
                threshold = number(data["threshold"], "adaptive threshold")
                if printed_threshold != f"{threshold:.2f}" or data[
                    "score_type"
                ] not in ("full", "domain"):
                    raise AnnotationError(
                        "KOfam output and pinned threshold metadata disagree"
                    )
            elif printed_threshold or mark:
                raise AnnotationError("KOfam threshold-unavailable hit cannot pass")
            if definition != data["definition"]:
                raise AnnotationError(
                    "KOfam definition differs from pinned KO metadata"
                )
            accepted = mark == "*"
            state = (
                "threshold_pass"
                if accepted
                else ("threshold_fail" if available else "threshold_unavailable")
            )
            result.mapped.add(gene)
            if accepted:
                result.add("ko", gene, [ko])
                result.definitions.setdefault("ko", {})[ko] = definition
            evidence.append(
                dict(
                    accession=protein["accession"],
                    gene_id=gene,
                    protein_id=protein["protein_id"],
                    ko=ko,
                    native_marker=mark,
                    threshold=data["threshold"],
                    printed_threshold=printed_threshold,
                    score_type=data["score_type"],
                    score=score,
                    evalue=evalue,
                    definition=definition,
                    accepted=str(accepted).lower(),
                    state=state,
                )
            )
    result.reported_hits = len(evidence)
    columns = COLUMNS
    result.tables["kofam_hits.tsv"] = columns, evidence
    return result


TABLE_COLUMNS = {"kofam_hits.tsv": COLUMNS}
