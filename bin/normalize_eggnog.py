"""Normalize the pinned eggNOG v3 schema and native GO namespace evidence."""

from __future__ import annotations

import csv
import json
import re
import sqlite3
from pathlib import Path

from annotation_common import (
    AnnotationError,
    integer,
    number,
    query,
    validate_accession,
)
from annotation_normalization import (
    Normalized,
    json_cell,
    native_rows,
    ordered_categories,
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


SEED_NATIVE_COLUMNS = (
    "qseqid",
    "sseqid",
    "evalue",
    "bitscore",
    "qstart",
    "qend",
    "sstart",
    "send",
    "pident",
    "qcov",
    "scov",
)


def query_lookup(proteins: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    """Index only canonical native tool IDs; original protein IDs are sample-local."""
    if not proteins:
        raise AnnotationError("eggNOG input contains no declared proteins")
    lookup, genes = {}, set()
    for protein in proteins:
        identifier, gene = protein.get("tool_id"), protein.get("gene_id")
        validate_accession(protein.get("accession"))
        if (
            not isinstance(identifier, str)
            or not identifier
            or any(char.isspace() for char in identifier)
            or not isinstance(gene, str)
            or not gene
        ):
            raise AnnotationError("Invalid canonical eggNOG query or gene ID")
        if identifier in lookup or gene in genes:
            raise AnnotationError("Duplicate canonical eggNOG query or gene ID")
        lookup[identifier] = protein
        genes.add(gene)
    return lookup


def seed_hits(
    raw: Path,
    lookup: dict[str, dict[str, str]],
    resource: Path,
    result: Normalized,
) -> dict[str, dict[str, str]]:
    """Separate mapped seed hits from proteins receiving functional annotations."""
    columns = SEED_NATIVE_COLUMNS
    text = (raw / "eggnog.emapper.seed_orthologs").read_text()
    lines = text.splitlines()
    if not text.endswith("\n") or [
        line for line in lines if line.startswith("#") and not line.startswith("##")
    ] != ["#" + "\t".join(columns)]:
        raise AnnotationError(
            "Missing native eggNOG seed header or complete final line"
        )
    rows = [line.split("\t") for line in lines if line and not line.startswith("#")]
    if [line for line in lines if re.fullmatch(r"## [0-9]+ queries scanned", line)] != [
        f"## {len(rows)} queries scanned"
    ]:
        raise AnnotationError(
            "Missing or inconsistent native eggNOG seed completion count"
        )
    seeds, evidence = {}, []
    with sqlite3.connect(
        f"file:{resource / 'eggnog.db'}?mode=ro", uri=True
    ) as database:
        for values in rows:
            if len(values) != len(columns):
                raise AnnotationError("Truncated eggNOG seed ortholog row")
            row = dict(zip(columns, values, strict=True))
            protein = query(lookup, row["qseqid"])
            if row["qseqid"] in seeds:
                raise AnnotationError("Duplicate eggNOG seed query")
            seed_id = integer(row["sseqid"], "eggNOG integer seed ID")
            for key in ("evalue", "bitscore"):
                number(row[key], key, minimum=0)
            for key in ("pident", "qcov", "scov"):
                if number(row[key], key, minimum=0) > 100:
                    raise AnnotationError(
                        "eggNOG seed identity/coverage exceeds 100 percent"
                    )
            for key in ("qstart", "qend", "sstart", "send"):
                integer(row[key], key, minimum=1)
            if (
                not 1
                <= int(row["qstart"])
                <= int(row["qend"])
                <= int(protein["length"])
            ):
                raise AnnotationError(
                    "eggNOG seed coordinates exceed the input protein"
                )
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
        ("accession", "gene_id", *columns, "seed_display"),
        evidence,
    )
    result.reported_hits = len(rows)
    return seeds


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
    raw: Path, proteins: list[dict[str, str]], resource: Path
) -> dict[str, Normalized]:
    """Parse one complete native batch once and project its evidence by accession.

    Native files and their completion markers remain untouched. Each member
    receives normalized tables in native row order, including header-only tables
    for members without hits. Used definitions and validated empty GO evidence
    remain specific to each member.
    """
    complete = normalize(raw, proteins, resource)
    accessions = {protein["gene_id"]: protein["accession"] for protein in proteins}
    members = {
        accession: Normalized() for accession in dict.fromkeys(accessions.values())
    }

    def member_for_gene(gene: str) -> Normalized:
        if gene not in accessions:
            raise AnnotationError("Normalized eggNOG evidence contains an unknown gene")
        return members[accessions[gene]]

    for name, (columns, rows) in complete.tables.items():
        for member in members.values():
            member.tables[name] = columns, []
        for row in rows:
            member = member_for_gene(row["gene_id"])
            if row["accession"] != accessions[row["gene_id"]]:
                raise AnnotationError(
                    "Normalized eggNOG row has a mismatched accession"
                )
            member.tables[name][1].append(row)
            if name == "eggnog_annotations.tsv" and "GOs" in json.loads(
                row["accepted_fields"]
            ):
                # The single-proteome parser records valid empty GO evidence as
                # an empty definition dictionary; absent/invalid GO has no key.
                member.definitions.setdefault("go", {})
    for name, genes in complete.features.items():
        for gene, values in genes.items():
            member_for_gene(gene).features.setdefault(name, {})[gene] = set(values)
    for gene in complete.mapped:
        member_for_gene(gene).mapped.add(gene)
    for gene in complete.accepted:
        member_for_gene(gene).accepted.add(gene)
    for error in complete.errors:
        member_for_gene(error["gene_id"]).errors.append(dict(error))
    for member in members.values():
        member.reported_hits = len(member.tables["eggnog_seed_hits.tsv"][1])
        for name, definitions in complete.definitions.items():
            used = {
                term
                for terms in member.features.get(name, {}).values()
                for term in terms
            }
            selected = {
                term: definitions[term] for term in sorted(used) if term in definitions
            }
            if selected:
                member.definitions[name] = selected
    if (
        sum(member.reported_hits for member in members.values())
        != complete.reported_hits
    ):
        raise AnnotationError("Batch projection does not conserve reported seed hits")
    return members


TABLE_COLUMNS = {"eggnog_annotations.tsv": COLUMNS}
TABLE_COLUMNS["eggnog_seed_hits.tsv"] = (
    "accession",
    "gene_id",
    *SEED_NATIVE_COLUMNS,
    "seed_display",
)
