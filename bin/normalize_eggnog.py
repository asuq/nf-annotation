"""Normalize the pinned eggNOG v3 schema and native GO namespace evidence."""

from __future__ import annotations

import csv
import re
import sqlite3
from pathlib import Path
from typing import Any

from annotation_common import (
    AnnotationError,
    integer,
    number,
    query,
)
from annotation_normalization import (
    Normalized,
    json_cell,
    native_rows,
    ordered_categories,
)
from eggnog_native import (
    SEED_NATIVE_COLUMNS,
    member_name,
    query_lookup,
    read_native_seeds,
    validate_execution,
)

FIELDS = (
    "Preferred_name",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "KEGG_rclass",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "BiGG_Reaction",
    "PFAMs",
)
HEADER = (
    "query",
    "seed_ortholog",
    "evalue",
    "score",
    "eggNOG_OGs",
    "tax_ceiling",
    "farthest_donor_lineage",
    "COG_category",
    *FIELDS,
    "annotation_confidence",
)
GO_FIELDS = ("gos_mf", "gos_bp", "gos_cc")
GO_HEADER = (
    "query",
    *[key for name in GO_FIELDS for key in (name, name + "_confidence")],
)
PATTERNS = {
    "KEGG_ko": r"K[0-9]{5}",
    "GOs": r"GO:[0-9]{7}",
    "EC": r"[1-7]\.(?:[0-9]+|-)\.(?:[0-9]+|-)\.(?:[0-9]+|-)",
}
FEATURE_NAMES = {"KEGG_ko": "ko", "GOs": "go", "EC": "ec", "PFAMs": "pfam"}
# COG's official vocabulary; X (mobilome) is included, and R/S remain valid.
COG_CATEGORIES = set("JAKLBDYVTMNZWUOCEFGHIPQRSX")
GO_NAMESPACES = dict(
    gos_mf="molecular_function",
    gos_bp="biological_process",
    gos_cc="cellular_component",
)


def ontology_terms(path: Path) -> dict[str, tuple[str, str]]:
    """Read live GO IDs and aliases without creating any ancestor assignments."""
    terms: dict[str, tuple[str, str]] = {}
    stanza: dict[str, list[str]] = {}
    in_term = False

    def finish() -> None:
        if not in_term or stanza.get("is_obsolete") == ["true"]:
            return
        if any(len(stanza.get(key, [])) != 1 for key in ("id", "name", "namespace")):
            raise AnnotationError("Incomplete or duplicate GO term metadata")
        namespace, name = stanza["namespace"][0], stanza["name"][0]
        if namespace not in GO_NAMESPACES.values():
            raise AnnotationError("Unknown GO namespace in prepared ontology")
        for identifier in stanza["id"] + stanza.get("alt_id", []):
            if not re.fullmatch(r"GO:[0-9]{7}", identifier) or identifier in terms:
                raise AnnotationError("Invalid or duplicate GO term identifier")
            terms[identifier] = namespace, name

    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith("["):
                finish()
                stanza = {}
                in_term = line == "[Term]"
            elif in_term and ": " in line:
                key, value = line.split(": ", 1)
                if key in ("id", "alt_id", "name", "namespace", "is_obsolete"):
                    stanza.setdefault(key, []).append(value)
        finish()
    if not terms:
        raise AnnotationError("Prepared GO ontology contains no live terms")
    return terms


def values(raw: str, name: str) -> list[str]:
    """Decode a native comma list, preserving missing and malformed distinctions."""
    if raw == "-":
        return []
    terms = raw.split(",")
    if any(not value or value.strip() != value for value in terms):
        raise AnnotationError("empty term or surrounding whitespace")
    if name == "EC":
        if any(not term.startswith("ec:") for term in terms):
            raise AnnotationError("missing native eggNOG v7 EC namespace")
        terms = [term.removeprefix("ec:") for term in terms]
    if name in PATTERNS and any(
        not re.fullmatch(PATTERNS[name], term) for term in terms
    ):
        raise AnnotationError(f"invalid {name} identifier")
    return terms


def _seed_evidence(
    rows: list[dict[str, str]],
    lookup: dict[str, dict[str, str]],
    resource: Path,
    result: Normalized,
) -> dict[str, dict[str, str]]:
    """Resolve validated native integer seeds without rewriting their input rows."""
    seeds, evidence = {}, []
    with sqlite3.connect(
        f"file:{resource / 'eggnog.db'}?mode=ro", uri=True
    ) as database:
        for native in rows:
            row = dict(native)
            protein = query(lookup, row["qseqid"])
            seed_id = integer(row["sseqid"], "eggNOG integer seed ID")
            display = database.execute(
                "SELECT name FROM protein_names WHERE id = ?", (seed_id,)
            ).fetchone()
            if display is None or not display[0]:
                raise AnnotationError(
                    "eggNOG search seed is absent from the annotation database"
                )
            row["seed_display"] = display[0]
            seeds[row["qseqid"]] = row
            result.mapped.add(protein["gene_id"])
            evidence.append(
                dict(accession=protein["accession"], gene_id=protein["gene_id"], **row)
            )
    result.tables["eggnog_seed_hits.tsv"] = (
        ("accession", "gene_id", *SEED_NATIVE_COLUMNS, "seed_display"),
        evidence,
    )
    result.reported_hits = len(rows)
    return seeds


def seed_hits(
    raw: Path,
    lookup: dict[str, dict[str, str]],
    resource: Path,
    result: Normalized,
) -> dict[str, dict[str, str]]:
    """Separate complete native seed evidence from functional annotation."""
    rows = read_native_seeds(raw / "eggnog.emapper.seed_orthologs", lookup)
    return _seed_evidence(rows, lookup, resource, result)


COLUMNS = (
    "accession",
    "gene_id",
    *HEADER,
    "confidence_decoded",
    "accepted_fields",
    "excluded_fields",
    "go_namespaces",
    "categories",
    "field_errors",
)


def normalize(raw: Path, proteins: list[dict[str, str]], resource: Path) -> Normalized:
    """Filter each native field independently; never infer GO namespace confidence."""
    result = Normalized()
    lookup = query_lookup(proteins)
    seeds = seed_hits(raw, lookup, resource, result)
    ontology = ontology_terms(resource / "go-basic.obo")
    return _normalize_annotations(raw, lookup, seeds, result, ontology)


def _normalize_annotations(
    raw: Path,
    lookup: dict[str, dict[str, str]],
    seeds: dict[str, dict[str, str]],
    result: Normalized,
    ontology: dict[str, tuple[str, str]],
) -> Normalized:
    """Interpret one proteome's native annotation and its own validated seed rows."""
    path = raw / "eggnog.emapper.annotations"
    comments, rows, header = [], [], None
    with path.open(newline="") as handle:
        for line in handle:
            if line.startswith("##"):
                comments.append(line.rstrip("\n"))
            elif line.startswith("#"):
                if header is not None or line.rstrip("\r\n").removeprefix("#").split(
                    "\t"
                ) != list(HEADER):
                    raise AnnotationError("Unexpected eggNOG v3 annotation header")
                header = HEADER
            elif line.strip():
                row = next(csv.reader([line], delimiter="\t"))
                if header is None or len(row) != len(HEADER):
                    raise AnnotationError(
                        "Missing header or truncated eggNOG annotation row"
                    )
                rows.append(dict(zip(HEADER, row, strict=True)))
    if header is None:
        raise AnnotationError("Missing eggNOG v3 annotation schema")
    footer_counts = [
        int(match.group(1))
        for line in comments
        if (match := re.fullmatch(r"## ([0-9]+) queries scanned", line))
    ]
    if len(footer_counts) != 1 or not len(rows) <= footer_counts[0] <= len(seeds):
        raise AnnotationError(
            "Missing or inconsistent eggNOG annotation completion count"
        )
    legends = [line for line in comments if line.startswith("## confidence codes:")]
    orders = [
        line for line in comments if line.startswith("## confidence field order:")
    ]
    if legends != [
        "## confidence codes: h=high m=medium l=low -=not annotated"
    ] or orders != ["## confidence field order: " + " ".join(FIELDS)]:
        raise AnnotationError("Missing or unsupported eggNOG confidence legend")
    go_rows = native_rows(
        raw / "eggnog.emapper.annotations.go_namespaces.tsv", GO_HEADER
    )
    go = {row["query"]: row for row in go_rows}
    if len(go) != len(go_rows):
        raise AnnotationError("Duplicate eggNOG GO sidecar query")
    seen, evidence = set(), []
    for row in rows:
        protein = query(lookup, row["query"])
        gene = protein["gene_id"]
        if gene in seen or row["query"] not in go:
            raise AnnotationError(
                "Duplicate eggNOG query or missing GO sidecar evidence"
            )
        seen.add(gene)
        if row["seed_ortholog"] in ("", "-", "ERROR"):
            raise AnnotationError("eggNOG annotation lacks a valid seed ortholog")
        seed = seeds.get(row["query"])
        if seed is None or row["seed_ortholog"] != seed["seed_display"]:
            raise AnnotationError(
                "eggNOG annotation seed differs from its retained search evidence"
            )
        if any(
            number(row[name], name) != number(seed[native_name], name)
            for name, native_name in (("evalue", "evalue"), ("score", "bitscore"))
        ):
            raise AnnotationError(
                "eggNOG annotation scores differ from its retained search evidence"
            )
        number(row["evalue"], "seed E-value", minimum=0)
        number(row["score"], "seed score")
        result.mapped.add(gene)
        before = len(result.errors)
        confidence = row["annotation_confidence"]
        accepted, decoded, excluded = {}, {}, {}
        for index, name in enumerate(FIELDS):
            field_name = FEATURE_NAMES.get(name, name)
            code = confidence[index] if len(confidence) == len(FIELDS) else "?"
            decoded[name] = {
                "h": "high",
                "m": "medium",
                "l": "low",
                "-": "not annotated",
            }.get(code, "invalid")
            try:
                terms = values(row[name], name)
                if code not in "hml-" or (bool(terms) != (code != "-")):
                    raise AnnotationError("confidence code and native value disagree")
                if name != "GOs":
                    if code in "hm":
                        accepted[name] = terms
                        result.add(field_name, gene, terms)
                    elif terms:
                        excluded[name] = "low confidence"
            except AnnotationError as error:
                result.error(
                    gene,
                    field_name,
                    json_cell({"value": row[name], "confidence": confidence}),
                    str(error),
                )
        # GO terms are accepted using the three already-computed namespace winners.
        merged, go_accepted, namespace_evidence = set(), set(), {}
        for name in GO_FIELDS:
            value, conf = go[row["query"]][name], go[row["query"]][name + "_confidence"]
            namespace_evidence[name] = dict(terms=value, confidence=conf)
            try:
                terms = values(value, "GOs")
                if any(
                    term not in ontology or ontology[term][0] != GO_NAMESPACES[name]
                    for term in terms
                ):
                    raise AnnotationError(
                        "GO term is absent, obsolete or assigned to the wrong namespace"
                    )
                if conf not in {"high", "medium", "low", "-"} or bool(terms) != (
                    conf != "-"
                ):
                    raise AnnotationError("namespace confidence and terms disagree")
                merged.update(terms)
                if conf in {"high", "medium"}:
                    go_accepted.update(terms)
            except AnnotationError as error:
                result.error(
                    gene, "go", json_cell(namespace_evidence[name]), str(error)
                )
        try:
            if merged != set(values(row["GOs"], "GOs")):
                raise AnnotationError("namespace terms disagree with native merged GOs")
        except AnnotationError as error:
            result.error(gene, "go", row["GOs"], str(error))
        if not any(error["field"] == "go" for error in result.errors[before:]):
            result.add("go", gene, go_accepted)
            result.definitions.setdefault("go", {}).update(
                {term: ontology[term][1] for term in go_accepted}
            )
            accepted["GOs"] = sorted(go_accepted)
        try:
            categories = (
                []
                if row["COG_category"] == "-"
                else ordered_categories(row["COG_category"], COG_CATEGORIES)
            )
        except AnnotationError as error:
            categories = None
            result.error(gene, "categories", row["COG_category"], str(error))
        evidence.append(
            dict(
                accession=protein["accession"],
                gene_id=gene,
                **row,
                confidence_decoded=json_cell(decoded),
                accepted_fields=json_cell(accepted),
                excluded_fields=json_cell(excluded),
                go_namespaces=json_cell(namespace_evidence),
                categories=json_cell(categories),
                field_errors=json_cell(result.errors[before:]),
            )
        )
    if set(go) != {row["query"] for row in rows}:
        raise AnnotationError("eggNOG GO sidecar and annotation query sets differ")
    columns = COLUMNS
    result.tables["eggnog_annotations.tsv"] = columns, evidence
    return result


def normalize_batch(
    raw: Path,
    proteins: list[dict[str, str]],
    resource: Path,
    *,
    batch: dict[str, Any],
) -> dict[str, Normalized]:
    """Normalize a shared search and genuinely separate per-proteome annotations."""
    partitions = validate_execution(raw, batch, proteins)
    lookup = query_lookup(proteins)
    search = Normalized()
    seeds = _seed_evidence(
        [row for rows in partitions.values() for row in rows],
        lookup,
        resource,
        search,
    )
    ontology = ontology_terms(resource / "go-basic.obo")
    columns, evidence = search.tables["eggnog_seed_hits.tsv"]
    grouped_proteins = {member["accession"]: [] for member in batch["members"]}
    grouped_evidence = {accession: [] for accession in grouped_proteins}
    for protein in proteins:
        if protein["accession"] not in grouped_proteins:
            raise AnnotationError("Protein belongs to an undeclared annotation member")
        grouped_proteins[protein["accession"]].append(protein)
    for row in evidence:
        grouped_evidence[row["accession"]].append(row)
    results = {}
    for index, member in enumerate(batch["members"]):
        accession = member["accession"]
        selected = grouped_evidence[accession]
        result = Normalized()
        result.tables["eggnog_seed_hits.tsv"] = columns, selected
        result.reported_hits = len(selected)
        result.mapped = {row["gene_id"] for row in selected}
        results[accession] = _normalize_annotations(
            raw / "annotations" / member_name(index),
            query_lookup(grouped_proteins[accession]),
            {row["qseqid"]: seeds[row["qseqid"]] for row in selected},
            result,
            ontology,
        )
    return results


TABLE_COLUMNS = {"eggnog_annotations.tsv": COLUMNS}
TABLE_COLUMNS["eggnog_seed_hits.tsv"] = (
    "accession",
    "gene_id",
    *SEED_NATIVE_COLUMNS,
    "seed_display",
)
