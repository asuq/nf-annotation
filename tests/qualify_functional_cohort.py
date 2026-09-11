#!/usr/bin/env python3
"""Reconcile native v0.4 cohort outputs against the recorded biological controls."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from decimal import Decimal
from pathlib import Path

from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
from annotation_common import (
    TOOLS,
    AnnotationError,
    bundle_proteins,
    read_json,
    read_tsv,
)
from annotation_source import validate_source
from annotation_summary import COG_CATEGORIES, MATRICES
from validate_inputs import detect_metadata_key_column


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AnnotationError(message)


def native_integer(value: object, name: str, minimum: int) -> int:
    require(
        type(value) is int and value >= minimum, f"Invalid native {name}: {value!r}"
    )
    return value


def crispr_coordinates(payload: dict, lengths: dict[str, int]) -> list[dict]:
    """Join the wrapper's explicit version-stripped IDs to source coordinates."""
    aliases = {}
    for source_id, length in lengths.items():
        native_id = re.sub(r"\.[0-9]+$", "", source_id)
        require(native_id not in aliases, f"Ambiguous CRISPR contig alias: {native_id}")
        aliases[native_id] = (source_id, length)
    sequences = payload.get("Sequences")
    require(isinstance(sequences, list), "Missing native CRISPR sequence list")
    seen, names, arrays = set(), set(), []
    for sequence in sequences:
        require(isinstance(sequence, dict), "Malformed native CRISPR sequence")
        native_id = sequence.get("Id")
        require(
            native_id in aliases and native_id not in seen,
            f"Unknown/duplicate CRISPR contig: {native_id}",
        )
        seen.add(native_id)
        source_id, length = aliases[native_id]
        require(
            native_integer(sequence.get("Length"), "contig length", 1) == length,
            f"CRISPR contig length differs: {source_id}",
        )
        candidates = sequence.get("Crisprs")
        require(isinstance(candidates, list), f"Missing CRISPR candidates: {source_id}")
        for candidate in candidates:
            require(isinstance(candidate, dict), "Malformed native CRISPR candidate")
            level = native_integer(candidate.get("Evidence_Level"), "evidence level", 1)
            require(level <= 4, "Unknown native CRISPR evidence level")
            # Preserve the established summary policy: only level 1 is excluded.
            if level == 1:
                continue
            start = native_integer(candidate.get("Start"), "array start", 1)
            end = native_integer(candidate.get("End"), "array end", 1)
            spacers = native_integer(candidate.get("Spacers"), "spacer count", 0)
            require(start <= end <= length, f"CRISPR array lies outside {source_id}")
            name = candidate.get("Name")
            require(
                isinstance(name, str) and bool(name) and (native_id, name) not in names,
                "Missing/duplicate native CRISPR array identity",
            )
            names.add((native_id, name))
            arrays.append(
                dict(
                    contig_id=native_id,
                    source_contig_id=source_id,
                    crispr_id=name,
                    start=start,
                    end=end,
                    evidence_level=level,
                    spacer_count=spacers,
                )
            )
    require(seen == set(aliases), "Native CRISPR results omit source contigs")
    return arrays


def matrix_counts(proteins: list[dict[str, str]], column: str) -> tuple[Counter, bool]:
    """Count genes independently of the aggregation implementation's counters."""
    counts = Counter()
    unavailable = False
    for row in proteins:
        value = row[column]
        if value == "NA":
            unavailable = True
            continue
        features = json.loads(value)
        require(
            isinstance(features, list)
            and all(isinstance(feature, str) and feature for feature in features),
            f"Invalid protein feature cell: {column}",
        )
        require(
            len(features) == len(set(features)),
            f"Duplicate protein feature assignment: {column}",
        )
        counts.update(features)
    return counts, unavailable


def category_counts(proteins: list[dict[str, str]]) -> tuple[Counter, dict, bool]:
    """Split each gene's unit weight equally among its declared COG categories."""
    counts, unavailable = matrix_counts(proteins, "cogclassifier_categories")
    require(
        set(counts) <= set(COG_CATEGORIES), "Unknown COG category in protein summary"
    )
    weights = dict.fromkeys(COG_CATEGORIES, Decimal(0))
    for row in proteins:
        if row["cogclassifier_categories"] == "NA":
            continue
        categories = json.loads(row["cogclassifier_categories"])
        for category in categories:
            weights[category] += Decimal(1) / len(categories)
    return counts, weights, unavailable


def qualify(results: Path, controls: list[dict[str, str]]) -> dict:
    manifest = validate_source(results)
    accessions = [row["accession"] for row in controls]
    require(manifest["complete"] is True, "Annotation acceptance is incomplete")
    require(
        manifest["accessions"] == accessions,
        "Qualified cohort differs from requested controls",
    )
    require(
        set(manifest["enabled_tools"]) == set(TOOLS),
        "Qualification requires all five core callers",
    )
    master = read_tsv(results / "tables/master_table.tsv")
    master_key = detect_metadata_key_column(list(master[0]))
    by_accession = {row[master_key]: row for row in master}
    summaries = read_tsv(results / "tables/protein_function_summary.tsv")
    proteins_by_accession = defaultdict(list)
    for row in summaries:
        proteins_by_accession[row["accession"]].append(row)
    require(
        set(proteins_by_accession) == set(accessions), "Protein summaries omit a genome"
    )
    reports = {}
    for control in controls:
        accession = control["accession"]
        row = by_accession[accession]
        sample = results / "samples" / accession
        metadata, proteins = bundle_proteins(sample / "annotation/bundle")
        expected_genes = {protein["gene_id"] for protein in proteins}
        observed = [protein["gene_id"] for protein in proteins_by_accession[accession]]
        require(
            len(observed) == len(set(observed)) and set(observed) == expected_genes,
            f"Protein summary identity mismatch: {accession}",
        )
        require(
            row["Gcode"] == control["expected_gcode"],
            f"Genetic-code regression: {accession}",
        )
        four, eleven = (
            Decimal(row["Average_Gene_Length_gcode4"]),
            Decimal(row["Average_Gene_Length_gcode11"]),
        )
        require(
            four.is_finite() and eleven.is_finite() and four > 0 and eleven > 0,
            "Invalid paired gene lengths",
        )
        require(
            row["Gcode"] == ("4" if four > Decimal("1.5") * eleven else "11"),
            "Genetic-code decision does not match the paired native values",
        )
        lengths = {}
        for record in SeqIO.parse(sample / "annotation/bundle/genome.fasta", "fasta"):
            require(record.id not in lengths, "Duplicate source FASTA contig")
            lengths[record.id] = len(record.seq)
        arrays = crispr_coordinates(read_json(sample / "ccfinder/result.json"), lengths)
        retained = read_tsv(sample / "ccfinder/ccfinder_crisprs.tsv")
        keys = (
            "contig_id",
            "crispr_id",
            "start",
            "end",
            "evidence_level",
            "spacer_count",
        )
        expected_rows = {tuple(str(array[key]) for key in keys) for array in arrays}
        observed_rows = [tuple(array[key] for key in keys) for array in retained]
        require(
            len(observed_rows) == len(set(observed_rows))
            and set(observed_rows) == expected_rows,
            "Published CRISPR coordinates differ from native retained arrays",
        )
        total_spacers = sum(array["spacer_count"] for array in arrays)
        require(
            len(arrays)
            == int(control["reported_crispr_arrays"])
            == int(row["CRISPRS"]),
            f"CRISPR-array regression: {accession}",
        )
        require(
            total_spacers
            == int(control["reported_spacers"])
            == int(row["SPACERS_SUM"]),
            f"CRISPR-spacer regression: {accession}",
        )
        coordinates = read_tsv(sample / "annotation/bundle/gene_coordinates.tsv")
        for coordinate in coordinates:
            source_id = coordinate["source_contig_id"]
            require(
                source_id in lengths
                and 1
                <= int(coordinate["start"])
                <= int(coordinate["end"])
                <= lengths[source_id],
                "CDS cannot join to original genome coordinates",
            )
        overlaps = [
            dict(
                contig_id=array["source_contig_id"],
                crispr_id=array["crispr_id"],
                gene_ids=sorted(
                    {
                        coordinate["gene_id"]
                        for coordinate in coordinates
                        if coordinate["source_contig_id"] == array["source_contig_id"]
                        and int(coordinate["start"]) <= array["end"]
                        and int(coordinate["end"]) >= array["start"]
                    }
                ),
            )
            for array in arrays
        ]
        for tool in TOOLS:
            require(
                row[f"{tool}_status"] == "success",
                f"Native {tool} is incomplete: {accession}",
            )
            require(
                int(row[f"{tool}_analysed_proteins"]) == len(proteins),
                f"Native {tool} input count differs: {accession}",
            )
        counts, weights, unavailable = category_counts(proteins_by_accession[accession])
        for category in COG_CATEGORIES:
            count = row[f"cog_category_{category}_gene_count"]
            weight = row[f"cog_category_{category}_fractional_gene_count"]
            if unavailable:
                require(count == weight == "NA", "Invalid COG categories became counts")
            else:
                require(
                    int(count) == counts[category], "COG category gene count differs"
                )
                # Published fractional counts round to 12 decimal places.
                require(
                    abs(Decimal(weight) - weights[category])
                    <= Decimal("0.0000000000005"),
                    "Fractional COG category count differs",
                )
        reports[accession] = dict(
            genetic_code=row["Gcode"],
            protein_count=len(proteins),
            input_id=metadata["input_id"],
            crispr_arrays=len(arrays),
            spacers=total_spacers,
            source_contigs=len(lengths),
            cds_segments=len(coordinates),
            crispr_overlapping_genes=overlaps,
            sixteen_s=row["16S"],
            codetta=row["Codetta_Genetic_Code"],
            field_errors=row["annotation_field_errors"],
        )
    for matrix, (tool, field) in MATRICES.items():
        path = results / f"tables/functional_matrices/{matrix}.tsv"
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            features = reader.fieldnames[1:]
            rows = list(reader)
        require(
            [row["accession"] for row in rows] == accessions,
            f"Matrix row order differs: {matrix}",
        )
        require(
            features == sorted(set(features)), f"Matrix feature order differs: {matrix}"
        )
        observed_features = set()
        for row in rows:
            protein_rows = proteins_by_accession[row["accession"]]
            counts, unavailable = matrix_counts(protein_rows, f"{tool}_{field}")
            observed_features.update(counts)
            master_row = by_accession[row["accession"]]
            assignments = "NA" if unavailable else str(sum(counts.values()))
            assigned_proteins = (
                "NA"
                if unavailable
                else str(
                    sum(
                        bool(json.loads(protein[f"{tool}_{field}"]))
                        for protein in protein_rows
                    )
                )
            )
            require(
                master_row[f"{tool}_{field}_assignments"] == assignments
                and master_row[f"{tool}_{field}_proteins"] == assigned_proteins,
                f"Protein/master feature counts differ: {matrix}, {row['accession']}",
            )
            for feature in features:
                expected = "NA" if unavailable else str(counts[feature])
                require(
                    row[feature] == expected,
                    f"Protein/matrix mismatch: {matrix}, {row['accession']}, {feature}",
                )
        require(
            set(features) == observed_features,
            f"Matrix omits or invents accepted features: {matrix}",
        )
    systems = read_tsv(results / "tables/defence_systems.tsv")
    system_ids = [(row["accession"], row["system_id"]) for row in systems]
    require(len(system_ids) == len(set(system_ids)), "Duplicate PADLOC system identity")
    for accession in accessions:
        require(
            int(by_accession[accession]["padloc_systems"])
            == sum(row["accession"] == accession for row in systems),
            "PADLOC system/master counts differ",
        )
    return dict(
        complete=True,
        accessions=accessions,
        genomes=reports,
        matrices=list(MATRICES),
        control_interpretation="Historical genetic-code and CRISPR regression controls; not experimental biological truth",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument(
        "--controls", type=Path, default=ROOT / "assets/tables/functional/cohort.tsv"
    )
    parser.add_argument(
        "--accession",
        action="append",
        help="Qualify a declared subset, in control-table order.",
    )
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    controls = read_tsv(args.controls)
    if args.accession:
        require(
            set(args.accession) <= {row["accession"] for row in controls},
            "Unknown requested control accession",
        )
        controls = [row for row in controls if row["accession"] in args.accession]
    try:
        report = qualify(args.results, controls)
    except (AnnotationError, OSError, ValueError) as error:
        report = dict(complete=False, error=str(error))
        args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(str(error), file=sys.stderr)
        return 1
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"Qualified {len(controls)} genomes and six source-specific matrices")
    return 0


if __name__ == "__main__":
    sys.exit(main())
