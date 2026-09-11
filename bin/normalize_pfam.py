"""Join retained HMMER domain evidence to native PfamScan clan resolution."""

from __future__ import annotations

import re
from pathlib import Path

from annotation_common import AnnotationError, integer, number, protein_lookup, query
from annotation_normalization import Normalized, json_cell


def scan_rows(path: Path, resolved: bool) -> dict[tuple[str, ...], list[str]]:
    """Require the qualified PfamScan format, including its declared overlap mode."""
    lines = path.read_text().splitlines()
    mode = "on" if resolved else "off"
    if (
        not any(line.startswith("# <seq id> <alignment start>") for line in lines)
        or f"#    resolve clan overlaps: {mode}" not in lines
    ):
        raise AnnotationError(f"Invalid PfamScan metadata: {path}")
    rows = {}
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        values = line.split()
        if len(values) != 15 or values[13] not in {"0", "1"}:
            raise AnnotationError(f"Malformed PfamScan row: {path}")
        key = tuple(values[index] for index in (0, 5, 1, 2, 3, 4, 8, 9))
        if key in rows:
            raise AnnotationError("Duplicate PfamScan domain identity")
        rows[key] = values
    return rows


COLUMNS = (
    "accession",
    "gene_id",
    "protein_id",
    "pfam_accession",
    "family",
    "model",
    "clan",
    "sequence_evalue",
    "sequence_score",
    "conditional_evalue",
    "independent_evalue",
    "domain_score",
    "domain_number",
    "domain_total",
    "protein_start",
    "protein_end",
    "envelope_start",
    "envelope_end",
    "hmm_start",
    "hmm_end",
    "hmm_length",
    "accepted",
    "state",
    "native_scan",
)


def normalize(raw: Path, proteins: list[dict[str, str]], resource: Path) -> Normalized:
    """Retain raw alternatives; accept only significant, natively resolved domains."""
    result, lookup = Normalized(), protein_lookup(proteins)
    original = scan_rows(raw / "pfam.raw.tsv", False)
    resolved = scan_rows(raw / "pfam.resolved.tsv", True)
    if any(
        key not in original or original[key] != row for key, row in resolved.items()
    ):
        raise AnnotationError(
            "Resolved PfamScan results are not a subset of the raw results"
        )
    metadata, current = {}, {}
    with (resource / "Pfam-A.hmm.dat").open() as handle:
        for line in handle:
            if line.startswith(("#=GF ID ", "#=GF AC ", "#=GF DE ", "#=GF CL ")):
                _, name, value = line.split(maxsplit=2)
                current[name] = value.strip()
            elif line.strip() == "//":
                if (
                    "ID" not in current
                    or "AC" not in current
                    or current["AC"] in metadata
                ):
                    raise AnnotationError("Invalid Pfam model metadata")
                metadata[current["AC"]] = current
                current = {}
    path = raw / "hmmscan.domtblout"
    lines = path.read_text().splitlines()
    if not any(
        line.startswith("# Program:") and line.split()[-1] == "hmmscan"
        for line in lines
    ) or not any(line.strip() == "# [ok]" for line in lines):
        raise AnnotationError("Missing complete HMMER hmmscan output footer")
    evidence, seen = [], set()
    for line in lines:
        if not line.strip() or line.startswith("#"):
            continue
        columns = line.split(maxsplit=22)
        if len(columns) != 23:
            raise AnnotationError("Truncated HMMER domain row")
        (
            model,
            accession,
            model_length,
            identifier,
            _,
            protein_length,
            seq_evalue,
            seq_score,
            seq_bias,
            ordinal,
            total,
            conditional_evalue,
            independent_evalue,
            domain_score,
            domain_bias,
            hmm_start,
            hmm_end,
            ali_start,
            ali_end,
            env_start,
            env_end,
            accuracy,
            description,
        ) = columns
        protein = query(lookup, identifier)
        gene = protein["gene_id"]
        if (
            not re.fullmatch(r"PF[0-9]{5}\.[0-9]+", accession)
            or accession not in metadata
            or metadata[accession]["ID"] != model
        ):
            raise AnnotationError("HMMER output contains an unknown Pfam model")
        lengths = [
            integer(value, "HMMER length", minimum=1)
            for value in (model_length, protein_length)
        ]
        if lengths[1] != int(protein["length"]):
            raise AnnotationError("HMMER query length differs from the input protein")
        for begin, end, limit in (
            (ali_start, ali_end, lengths[1]),
            (env_start, env_end, lengths[1]),
            (hmm_start, hmm_end, lengths[0]),
        ):
            if (
                not 1
                <= integer(begin, "domain start", minimum=1)
                <= integer(end, "domain end", minimum=1)
                <= limit
            ):
                raise AnnotationError(
                    "HMMER domain coordinates exceed the declared sequence"
                )
        if not int(env_start) <= int(ali_start) <= int(ali_end) <= int(env_end):
            raise AnnotationError("HMMER alignment is outside its envelope")
        if integer(ordinal, "domain number", minimum=1) > integer(
            total, "domain total", minimum=1
        ):
            raise AnnotationError("HMMER domain number exceeds the total")
        for value in (
            seq_evalue,
            seq_bias,
            conditional_evalue,
            independent_evalue,
            domain_bias,
        ):
            number(value, "HMMER probability or bias", minimum=0)
        for value in (seq_score, domain_score):
            number(value, "HMMER score")
        if number(accuracy, "alignment accuracy", minimum=0) > 1:
            raise AnnotationError("HMMER alignment accuracy exceeds one")
        key = (
            identifier,
            accession,
            ali_start,
            ali_end,
            env_start,
            env_end,
            hmm_start,
            hmm_end,
        )
        if key in seen:
            raise AnnotationError("Duplicate HMMER domain identity")
        seen.add(key)
        native = original.get(key)
        clan = metadata[accession].get("CL", "No_clan")
        if native and (
            native[6] != model or native[10] != model_length or native[14] != clan
        ):
            raise AnnotationError("PfamScan model metadata differs from HMMER/resource")
        accepted = key in resolved and resolved[key][13] == "1"
        state = (
            "accepted"
            if accepted
            else "clan_overlap"
            if key in original and key not in resolved
            else "below_gathering_threshold"
            if native
            else "not_reported_by_pfam_scan"
        )
        result.mapped.add(gene)
        if accepted:
            family = accession.split(".")[0]
            result.add("pfam", gene, [family])
            result.definitions.setdefault("pfam", {})[family] = metadata[accession].get(
                "DE", model
            )
        evidence.append(
            dict(
                accession=protein["accession"],
                gene_id=gene,
                protein_id=protein["protein_id"],
                pfam_accession=accession,
                family=accession.split(".")[0],
                model=model,
                clan=clan,
                sequence_evalue=seq_evalue,
                sequence_score=seq_score,
                conditional_evalue=conditional_evalue,
                independent_evalue=independent_evalue,
                domain_score=domain_score,
                domain_number=ordinal,
                domain_total=total,
                protein_start=ali_start,
                protein_end=ali_end,
                envelope_start=env_start,
                envelope_end=env_end,
                hmm_start=hmm_start,
                hmm_end=hmm_end,
                hmm_length=model_length,
                accepted=str(accepted).lower(),
                state=state,
                native_scan=json_cell(native),
            )
        )
    if not set(original) <= seen:
        raise AnnotationError(
            "PfamScan result lacks matching retained HMMER domain evidence"
        )
    result.reported_hits = len(evidence)
    columns = COLUMNS
    result.tables["pfam_domains.tsv"] = columns, evidence
    return result


TABLE_COLUMNS = {"pfam_domains.tsv": COLUMNS}
