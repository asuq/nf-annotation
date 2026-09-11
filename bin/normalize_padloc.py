"""Normalize native PADLOC systems and member genes with exact coordinate joins."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote

from annotation_common import (
    AnnotationError,
    integer,
    number,
    protein_lookup,
    query,
    read_tsv,
)
from annotation_normalization import Normalized, json_cell, native_rows

HEADER = (
    "system.number",
    "seqid",
    "system",
    "target.name",
    "hmm.accession",
    "hmm.name",
    "protein.name",
    "full.seq.E.value",
    "domain.iE.value",
    "target.coverage",
    "hmm.coverage",
    "start",
    "end",
    "strand",
    "target.description",
    "relative.position",
    "contig.end",
    "all.domains",
    "best.hits",
)


TABLE_COLUMNS = {
    "defence_systems.tsv": (
        "accession",
        "system_id",
        "native_system_number",
        "system",
        "contig_id",
        "member_genes",
        "gene_ids",
    ),
    "defence_genes.tsv": (
        "accession",
        "gene_id",
        "protein_id",
        "system_id",
        "segment",
        *HEADER,
    ),
}


def normalize(raw: Path, proteins: list[dict[str, str]], bundle: Path) -> Normalized:
    """Count native systems once and retain every member/model and genomic segment."""
    result, lookup = Normalized(), protein_lookup(proteins)
    lines = (raw / "input.domtblout").read_text().splitlines()
    if not any(line.strip() == "# [ok]" for line in lines) or not any(
        line.startswith("# Program:") and line.split()[-1] == "hmmsearch"
        for line in lines
    ):
        raise AnnotationError("PADLOC lacks a complete HMMER domain report")
    domain_identities = set()
    for line in lines:
        if line.strip() and not line.startswith("#"):
            values = line.split(maxsplit=22)
            if len(values) != 23:
                raise AnnotationError("Malformed PADLOC HMMER domain row")
            protein = query(lookup, values[0])
            if integer(values[2], "PADLOC query length", minimum=1) != int(
                protein["length"]
            ):
                raise AnnotationError("PADLOC HMMER query length differs from bundle")
            domain_identities.add((protein["gene_id"], values[4], values[3]))
            result.mapped.add(protein["gene_id"])
    path = raw / "input_padloc.csv"
    if path.is_file():
        native = native_rows(path, HEADER, delimiter=",")
        if not native:
            raise AnnotationError(
                "Unexpected empty PADLOC CSV; native zero results omit it"
            )
    else:
        log = (raw / "tool.log").read_text()
        if not re.search(r"\bNothing found for input(?:\s|$)", log):
            raise AnnotationError(
                "PADLOC CSV missing without its native zero-system declaration"
            )
        native = []
    coordinates = read_tsv(
        bundle / "gene_coordinates.tsv",
        ("gene_id", "contig_id", "start", "end", "strand", "segment"),
    )
    by_gene = defaultdict(list)
    for row in coordinates:
        by_gene[row["gene_id"]].append(row)
    systems, members, seen = {}, [], set()
    for row in native:
        protein = query(lookup, row["target.name"])
        gene, accession = protein["gene_id"], protein["accession"]
        number_id = integer(row["system.number"], "PADLOC system number", minimum=1)
        system_id = quote(accession, safe="") + "::padloc::" + str(number_id)
        if (
            not row["system"]
            or not row["hmm.accession"]
            or not row["hmm.name"]
            or not row["protein.name"]
        ):
            raise AnnotationError("PADLOC system/member identity is missing")
        if (gene, row["hmm.accession"], row["hmm.name"]) not in domain_identities:
            raise AnnotationError("PADLOC member lacks retained native HMM evidence")
        for name in ("full.seq.E.value", "domain.iE.value"):
            number(row[name], name, minimum=0)
        for name in ("target.coverage", "hmm.coverage"):
            if number(row[name], name, minimum=0) > 1:
                raise AnnotationError("PADLOC coverage exceeds one")
        for name in ("start", "end", "relative.position", "contig.end"):
            integer(row[name], name, minimum=1)
        if int(row["relative.position"]) > int(row["contig.end"]):
            raise AnnotationError(
                "PADLOC relative position exceeds its contig gene count"
            )
        segments = [
            segment
            for segment in by_gene[gene]
            if all(
                segment[key] == row[native_key]
                for key, native_key in (
                    ("contig_id", "seqid"),
                    ("start", "start"),
                    ("end", "end"),
                    ("strand", "strand"),
                )
            )
        ]
        if len(segments) != 1:
            raise AnnotationError(
                "PADLOC member does not join exactly one canonical coordinate segment"
            )
        key = (
            system_id,
            gene,
            segments[0]["segment"],
            row["hmm.accession"],
            row["hmm.name"],
        )
        if key in seen:
            raise AnnotationError("Duplicate PADLOC member/model/segment relationship")
        seen.add(key)
        system = systems.setdefault(
            system_id,
            dict(
                accession=accession,
                system_id=system_id,
                native_system_number=number_id,
                system=row["system"],
                contig_id=row["seqid"],
                genes=set(),
            ),
        )
        if system["system"] != row["system"] or system["contig_id"] != row["seqid"]:
            raise AnnotationError(
                "PADLOC system number has conflicting system or contig identities"
            )
        system["genes"].add(gene)
        result.add("systems", gene, [system_id])
        members.append(
            dict(
                accession=accession,
                gene_id=gene,
                protein_id=protein["protein_id"],
                system_id=system_id,
                segment=segments[0]["segment"],
                **row,
            )
        )
    system_rows = []
    for system in systems.values():
        genes = sorted(system.pop("genes"))
        system_rows.append(
            dict(**system, member_genes=len(genes), gene_ids=json_cell(genes))
        )
    result.reported_hits = len(native)
    result.tables["defence_systems.tsv"] = (
        TABLE_COLUMNS["defence_systems.tsv"],
        system_rows,
    )
    result.tables["defence_genes.tsv"] = TABLE_COLUMNS["defence_genes.tsv"], members
    return result
