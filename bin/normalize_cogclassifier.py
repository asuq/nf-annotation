"""Normalize native COGclassifier assignments without changing best-hit selection."""

from __future__ import annotations

import csv
from pathlib import Path

from annotation_common import (
    AnnotationError,
    integer,
    number,
    protein_lookup,
    query,
)
from annotation_normalization import (
    Normalized,
    json_cell,
    native_rows,
    ordered_categories,
)

NATIVE_COLUMNS = (
    "QUERY_ID",
    "COG_ID",
    "CDD_ID",
    "EVALUE",
    "IDENTITY",
    "GENE_NAME",
    "COG_NAME",
    "COG_LETTER",
    "COG_DESCRIPTION",
)

BLAST_COLUMNS = (
    "query",
    "subject",
    "identity",
    "alignment_length",
    "mismatches",
    "gapopens",
    "query_start",
    "query_end",
    "subject_start",
    "subject_end",
    "evalue",
    "bitscore",
)


COLUMNS = (
    "accession",
    "gene_id",
    "protein_id",
    "cog_id",
    "cdd_id",
    "native_primary_category",
    "categories_raw",
    "categories",
    "definition",
    "accepted",
    "reason",
    "evalue",
    "bitscore",
    "identity",
    "query_start",
    "query_end",
    "subject_start",
    "subject_end",
    "field_errors",
)


def normalize(raw: Path, proteins: list[dict[str, str]], resource: Path) -> Normalized:
    """Validate native assignments and restore the full pinned category strings."""
    result, lookup = Normalized(), protein_lookup(proteins)
    path = raw / "rpsblast.tsv"
    # Validate every reported ID/row, including alternatives ignored by the native classifier.
    with path.open(newline="") as handle:
        for line, values in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(values) != len(BLAST_COLUMNS):
                raise AnnotationError(f"Malformed RPS-BLAST row {line}")
            row = dict(zip(BLAST_COLUMNS, values, strict=True))
            protein = query(lookup, row["query"])
            if (
                not row["subject"].startswith("CDD:")
                or not row["subject"][4:].isdigit()
            ):
                raise AnnotationError("Invalid RPS-BLAST subject identifier")
            identity = number(row["identity"], "percent identity", minimum=0)
            if identity > 100:
                raise AnnotationError("RPS-BLAST identity exceeds 100 percent")
            for name in ("evalue", "bitscore"):
                number(row[name], name, minimum=0)
            for name in (
                "alignment_length",
                "query_start",
                "query_end",
                "subject_start",
                "subject_end",
            ):
                integer(row[name], name, minimum=1)
            for name in ("mismatches", "gapopens"):
                integer(row[name], name)
            if (
                not 1
                <= int(row["query_start"])
                <= int(row["query_end"])
                <= int(protein["length"])
            ):
                raise AnnotationError("RPS-BLAST protein coordinates exceed input")
            result.reported_hits += 1
            result.mapped.add(protein["gene_id"])
    assignments = native_rows(raw / "cogclassifier.native.tsv", NATIVE_COLUMNS)
    native_by_query = {row["QUERY_ID"]: row for row in assignments}
    if len(native_by_query) != len(assignments):
        raise AnnotationError("Duplicate native COGclassifier assignment")
    with (resource / "cog_func_category.tsv").open(newline="") as handle:
        vocabulary = {row[0] for row in csv.reader(handle, delimiter="\t")}
    with (resource / "cog_definition.tsv").open(newline="") as handle:
        definitions = {row[0]: row for row in csv.reader(handle, delimiter="\t")}
    with (resource / "cddid.tbl").open(newline="") as handle:
        mapping = {
            row[0]: row[1]
            for row in csv.reader(handle, delimiter="\t")
            if row[1].startswith("COG")
        }
    top_hits = {}
    with path.open(newline="") as handle:
        for values in csv.reader(handle, delimiter="\t"):
            top_hits.setdefault(
                values[0], dict(zip(BLAST_COLUMNS, values, strict=True))
            )
    evidence = []
    for identifier, aln in top_hits.items():
        protein = query(lookup, identifier)
        gene = protein["gene_id"]
        cdd_id = aln["subject"].removeprefix("CDD:")
        if cdd_id not in mapping:
            raise AnnotationError("COG hit is absent from the pinned CDD mapping")
        cog_id = mapping[cdd_id]
        definition = definitions.get(cog_id)
        fields, errors_before = None, len(result.errors)
        accepted = identifier in native_by_query
        if accepted != (definition is not None):
            raise AnnotationError(
                "Native COGclassifier assignment set disagrees with its definitions"
            )
        if accepted and (
            native_by_query[identifier]["COG_ID"] != cog_id
            or native_by_query[identifier]["CDD_ID"] != cdd_id
            or native_by_query[identifier]["COG_LETTER"] != definition[1][0]
        ):
            raise AnnotationError(
                "Native COGclassifier assignment differs from its first reported hit"
            )
        if definition is not None:
            try:
                fields = ordered_categories(definition[1], vocabulary)
                result.add("categories", gene, fields)
            except AnnotationError as error:
                result.error(gene, "categories", definition[1], str(error))
        if accepted:
            result.add("cog", gene, [cog_id])
            result.definitions.setdefault("cog", {})[cog_id] = definition[2]
        row = dict(
            accession=protein["accession"],
            gene_id=gene,
            protein_id=protein["protein_id"],
            cog_id=cog_id,
            cdd_id=cdd_id,
            native_primary_category=native_by_query[identifier]["COG_LETTER"]
            if accepted
            else None,
            categories_raw=definition[1] if definition else None,
            categories=json_cell(fields),
            definition=definition[2] if definition else None,
            accepted=str(accepted).lower(),
            reason="native_assignment" if accepted else "missing_native_definition",
            evalue=aln["evalue"],
            bitscore=aln["bitscore"],
            identity=aln["identity"],
            query_start=aln["query_start"],
            query_end=aln["query_end"],
            subject_start=aln["subject_start"],
            subject_end=aln["subject_end"],
            field_errors=json_cell(result.errors[errors_before:]),
        )
        evidence.append(row)
    if not set(native_by_query) <= set(top_hits):
        raise AnnotationError(
            "Native COGclassifier assignment lacks retained search evidence"
        )
    columns = COLUMNS
    result.tables["cog_assignments.tsv"] = columns, evidence
    return result


TABLE_COLUMNS = {"cog_assignments.tsv": COLUMNS}
