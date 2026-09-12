#!/usr/bin/env python3
"""Build a validated, reusable protein/coordinate bundle from Prokka outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import shutil
import sys
from collections import defaultdict
from io import StringIO
from pathlib import Path
from urllib.parse import quote, unquote

from annotation_common import (
    COORDINATE_COLUMNS,
    PROTEIN_COLUMNS,
    SCHEMA_VERSION,
    AnnotationError,
    digest,
    identity,
    integer,
    write_json,
    write_tsv,
)
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

LOGGER = logging.getLogger(__name__)


def fasta_records(text: str, label: str) -> list[SeqRecord]:
    """Read unique non-empty records without silently removing sequence whitespace."""
    if not text.strip().startswith(">"):
        reason = (
            "empty_proteome"
            if label == "protein input" and not text.strip()
            else "invalid_fasta"
        )
        raise AnnotationError(f"{reason}: empty or invalid FASTA: {label}")
    for line in text.splitlines():
        if line and not line.startswith(">") and any(char.isspace() for char in line):
            raise AnnotationError(f"Whitespace within FASTA sequence: {label}")
    try:
        records = list(SeqIO.parse(StringIO(text), "fasta"))
    except ValueError as exc:
        raise AnnotationError(f"Invalid FASTA: {label}: {exc}") from exc
    ids = [record.id for record in records]
    if (
        not records
        or len(ids) != len(set(ids))
        or any(not record.seq for record in records)
    ):
        raise AnnotationError(f"Empty sequences or duplicate identifiers: {label}")
    return records


def attributes(text: str) -> dict[str, str]:
    """Decode GFF3 attributes while rejecting duplicate keys."""
    result = {}
    for field in text.split(";"):
        if not field:
            continue
        if "=" not in field:
            raise AnnotationError(f"Invalid GFF attribute: {field!r}")
        key, value = field.split("=", 1)
        key, value = unquote(key), unquote(value)
        if key in result:
            raise AnnotationError(f"Duplicate GFF attribute: {key}")
        result[key] = value
    return result


def validate_genbank(
    gbk: Path,
    contigs: list[SeqRecord],
    proteins: list[SeqRecord],
    cds_by_protein: dict,
    gcode: int,
) -> dict:
    """Validate native translations/locations and retain biological segment order."""
    try:
        records = list(SeqIO.parse(gbk, "genbank"))
    except ValueError as exc:
        raise AnnotationError(f"Invalid GenBank: {exc}") from exc
    gff_sequences = {record.id: str(record.seq).upper() for record in contigs}
    protein_sequences = {record.id: str(record.seq) for record in proteins}
    by_sequence = defaultdict(list)
    for contig, sequence in gff_sequences.items():
        by_sequence[sequence].append(contig)
    seen_contigs, result = set(), {}
    for record in records:
        sequence = str(record.seq).upper()
        matches = {
            key
            for key in (record.id, record.name)
            if gff_sequences.get(key) == sequence
        }
        if not matches:
            matches = set(by_sequence.get(sequence, []))
        if len(matches) != 1:
            raise AnnotationError(
                f"Cannot uniquely map GenBank record {record.id} to GFF"
            )
        contig = matches.pop()
        if contig in seen_contigs:
            raise AnnotationError(f"Repeated GenBank contig: {contig}")
        seen_contigs.add(contig)
        topology = record.annotations.get("topology")
        if topology is not None and topology not in ("linear", "circular"):
            raise AnnotationError(f"Unknown GenBank topology: {topology}")
        for feature in record.features:
            if feature.type != "CDS":
                continue
            translations = feature.qualifiers.get("translation", [])
            loci = feature.qualifiers.get("locus_tag", [])
            if not translations and not set(loci) & protein_sequences.keys():
                continue  # Non-translated features remain in the original GenBank/GFF.
            if (
                len(loci) != 1
                or loci[0] not in protein_sequences
                or len(translations) != 1
            ):
                raise AnnotationError(
                    "GenBank translated CDS lacks a unique FAA locus/translation"
                )
            locus = loci[0]
            if locus in result:
                raise AnnotationError(f"Repeated translated CDS in GenBank: {locus}")
            if translations[0] != protein_sequences[locus]:
                raise AnnotationError(
                    f"GenBank translation disagrees with FAA: {locus}"
                )
            if feature.qualifiers.get("transl_table") != [str(gcode)]:
                raise AnnotationError(
                    f"GenBank translation table disagrees with selected code: {locus}"
                )
            if feature.location is None:
                raise AnnotationError(f"GenBank CDS has no location: {locus}")
            segments = cds_by_protein[locus]
            by_position = {}
            for segment in segments:
                position = (
                    segment["contig_id"],
                    segment["start"],
                    segment["end"],
                    segment["strand"],
                )
                if position in by_position:
                    raise AnnotationError(f"Duplicate CDS segment: {locus}")
                by_position[position] = segment
            ordered = []
            for part in feature.location.parts:
                if part.ref or part.ref_db or part.strand not in (-1, 1):
                    raise AnnotationError(
                        f"Unsupported remote or unstranded CDS location: {locus}"
                    )
                try:
                    position = (
                        contig,
                        int(part.start) + 1,
                        int(part.end),
                        "+" if part.strand == 1 else "-",
                    )
                except TypeError as exc:
                    raise AnnotationError(f"Unknown CDS position: {locus}") from exc
                if position not in by_position:
                    raise AnnotationError(
                        f"GenBank CDS coordinates disagree with GFF: {locus}"
                    )
                ordered.append(by_position.pop(position))
            if by_position:
                raise AnnotationError(f"GenBank is missing GFF CDS segments: {locus}")
            codon_start = feature.qualifiers.get("codon_start")
            if codon_start not in (["1"], ["2"], ["3"]):
                raise AnnotationError(f"Invalid GenBank codon_start: {locus}")
            if int(ordered[0]["phase"]) != int(codon_start[0]) - 1:
                raise AnnotationError(
                    f"GenBank codon_start disagrees with GFF phase: {locus}"
                )
            cds_by_protein[locus] = ordered
            result[locus] = {
                "topology": topology,
                "genbank_location": str(feature.location),
                "genbank_qualifiers": json.dumps(
                    {
                        key: value
                        for key, value in feature.qualifiers.items()
                        if key != "translation"
                    },
                    sort_keys=True,
                    separators=(",", ":"),
                    ensure_ascii=True,
                ),
            }
    if seen_contigs != set(gff_sequences):
        raise AnnotationError("GenBank contigs differ from GFF contigs")
    if set(result) != set(protein_sequences):
        raise AnnotationError(
            "GenBank translated CDS identifiers differ from FAA identifiers"
        )
    return result


def build_bundle(
    accession: str,
    gcode: int,
    genome: Path,
    faa: Path,
    gff: Path,
    gbk: Path,
    output: Path,
) -> dict:
    """Create the canonical protein bundle without altering native translations."""
    if output.exists() and any(output.iterdir()):
        raise AnnotationError("Bundle output must be a new empty directory")
    if (
        gcode not in (4, 11)
        or not accession.strip()
        or any(ord(c) < 32 for c in accession)
    ):
        raise AnnotationError(
            "Bundle requires a valid accession and genetic code 4 or 11"
        )
    proteins = fasta_records(faa.read_text(), "protein input")
    for protein in proteins:
        if set(str(protein.seq)) - set("ACDEFGHIKLMNPQRSTVWYBXZJUO"):
            raise AnnotationError(f"Invalid amino-acid symbols in {protein.id}")
    source_contigs = fasta_records(genome.read_text(), "source genome")
    for contig in source_contigs:
        if set(str(contig.seq).upper()) - set("ACGTRYSWKMBDHVN"):
            raise AnnotationError(f"Invalid DNA symbols in {contig.id}")
    text = gff.read_text()
    if not text.startswith("##gff-version 3") or "\n##FASTA\n" not in text:
        raise AnnotationError(
            "Prokka GFF3 must contain its authoritative FASTA sequences"
        )
    features_text, contigs_text = text.split("\n##FASTA\n", 1)
    contigs = fasta_records(contigs_text, "GFF genome")
    original_sequences = {
        record.id: str(record.seq).upper() for record in source_contigs
    }
    # Prokka 1.15.6 uppercases DNA and replaces non-ACGT IUPAC bases with N.
    # Source gaps/pads are rejected above: their removal would shift coordinates.
    masking = str.maketrans({base: "N" for base in "RYSWKMBDHV"})
    source_ids = {
        name: sequence.translate(masking)
        for name, sequence in original_sequences.items()
    }
    source_by_sequence = defaultdict(list)
    for name, sequence in source_ids.items():
        source_by_sequence[sequence].append(name)
    contig_sources = {}
    sequence_provenance = []
    lengths = {}
    for record in contigs:
        sequence = str(record.seq).upper()
        if source_ids.get(record.id) == sequence:
            source_id = record.id
        else:
            matches = source_by_sequence.get(sequence, [])
            if len(matches) != 1:
                raise AnnotationError(
                    f"Cannot uniquely map GFF contig {record.id} to the source genome"
                )
            source_id = matches[0]
        if source_id in contig_sources.values():
            raise AnnotationError(
                f"GFF contigs map repeatedly to source contig {source_id}"
            )
        contig_sources[record.id] = source_id
        lengths[record.id] = len(record.seq)
        original = original_sequences[source_id]
        sequence_provenance.append(
            {
                "contig_id": record.id,
                "source_contig_id": source_id,
                "length": len(sequence),
                "source_sequence_sha256": hashlib.sha256(original.encode()).hexdigest(),
                "native_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                "masked_iupac_bases": sum(
                    source != native
                    for source, native in zip(original, sequence, strict=True)
                ),
            }
        )

    protein_ids = {record.id for record in proteins}
    cds_by_protein = defaultdict(list)
    for line in features_text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) != 9:
            raise AnnotationError("GFF feature must contain nine columns")
        if fields[2] != "CDS":
            continue
        seqid, _, _, start, end, _, strand, phase, encoded = fields
        attr = attributes(encoded)
        matches = {
            attr.get(key) for key in ("ID", "locus_tag", "protein_id")
        } & protein_ids
        if len(matches) != 1:
            raise AnnotationError(
                f"CDS does not identify exactly one FAA protein: {encoded}"
            )
        protein_id = matches.pop()
        start, end = (
            integer(start, "CDS start", minimum=1),
            integer(end, "CDS end", minimum=1),
        )
        if seqid not in lengths or start > end or end > lengths[seqid]:
            raise AnnotationError(
                f"CDS coordinates outside source contig: {protein_id}"
            )
        if (
            strand not in ("+", "-")
            or phase not in ("0", "1", "2")
            or not attr.get("ID")
        ):
            raise AnnotationError(f"Invalid CDS strand, phase or ID: {protein_id}")
        cds_by_protein[protein_id].append(
            {
                "feature_id": attr["ID"],
                "contig_id": seqid,
                "source_contig_id": contig_sources[seqid],
                "start": start,
                "end": end,
                "strand": strand,
                "phase": phase,
                "attributes": encoded,
            }
        )
    if set(cds_by_protein) != protein_ids:
        raise AnnotationError(
            f"FAA proteins without CDS features: {sorted(protein_ids - set(cds_by_protein))}"
        )

    genbank_metadata = validate_genbank(gbk, contigs, proteins, cds_by_protein, gcode)

    output.mkdir(parents=True, exist_ok=True)
    rows, coordinate_rows, tool_ids = [], [], set()
    with (output / "proteins.faa").open("w") as handle:
        for protein in proteins:
            gene_id = f"{quote(accession, safe='')}::{quote(protein.id, safe='')}"
            tool_id = "p" + hashlib.sha256(gene_id.encode()).hexdigest()[:24]
            if tool_id in tool_ids:
                raise AnnotationError("Canonical tool-ID collision")
            tool_ids.add(tool_id)
            segments = cds_by_protein[protein.id]
            if len({(s["contig_id"], s["strand"]) for s in segments}) != 1:
                raise AnnotationError(
                    f"CDS segments disagree on contig/strand: {protein.id}"
                )
            signatures = [(s["start"], s["end"], s["phase"]) for s in segments]
            if len(signatures) != len(set(signatures)):
                raise AnnotationError(f"Duplicate CDS segment: {protein.id}")
            first = segments[0]
            rows.append(
                {
                    "accession": accession,
                    "gene_id": gene_id,
                    "protein_id": protein.id,
                    "tool_id": tool_id,
                    "description": protein.description,
                    "length": len(protein.seq),
                    "sha256": hashlib.sha256(str(protein.seq).encode()).hexdigest(),
                    "genetic_code": gcode,
                    "contig_id": first["contig_id"],
                    "source_contig_id": first["source_contig_id"],
                    "contig_length": lengths[first["contig_id"]],
                    "strand": first["strand"],
                    "segment_count": len(segments),
                    **genbank_metadata[protein.id],
                }
            )
            for index, segment in enumerate(segments, 1):
                coordinate_rows.append(
                    {
                        "accession": accession,
                        "gene_id": gene_id,
                        "segment": index,
                        **segment,
                    }
                )
            handle.write(f">{tool_id}\n{protein.seq}\n")
    write_tsv(output / "protein_manifest.tsv", PROTEIN_COLUMNS, rows)
    write_tsv(output / "gene_coordinates.tsv", COORDINATE_COLUMNS, coordinate_rows)
    for source, name in (
        (faa, "source.faa"),
        (gff, "source.gff"),
        (gbk, "source.gbk"),
        (genome, "genome.fasta"),
    ):
        shutil.copyfile(source, output / name)
    files = {
        path.name: digest(path)
        for path in sorted(output.iterdir())
        if path.is_file() and path.name != "bundle.json"
    }
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "accession": accession,
        "status": "success",
        "genetic_code": gcode,
        "protein_count": len(rows),
        "files": files,
        "input_id": identity(
            {
                "accession": accession,
                "genetic_code": gcode,
                "proteins": [
                    {key: row[key] for key in ("gene_id", "tool_id", "sha256")}
                    for row in rows
                ],
            }
        ),
        "coordinate_id": files["gene_coordinates.tsv"],
        "source_genome_sha256": digest(genome),
        "source_sequence_policy": "prokka_uppercase_iupac_to_n",
        "contig_sequence_provenance": sequence_provenance,
        "coordinate_convention": "1-based-closed; segments in GenBank biological order",
        "producer_sha256": digest(Path(__file__)),
        "common_contract_sha256": digest(
            Path(__file__).with_name("annotation_common.py")
        ),
    }
    write_json(output / "bundle.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--gcode", required=True, type=int)
    for option in ("genome", "faa", "gff", "gbk", "outdir"):
        parser.add_argument(f"--{option}", required=True, type=Path)
    parser.add_argument("--prokka-log", required=True, type=Path)
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args.outdir.mkdir(parents=True, exist_ok=True)
    if any(args.outdir.iterdir()):
        parser.error("Bundle output must be a new empty directory")
    try:
        codes = [
            line.split("=", 1)[1]
            for line in args.prokka_log.read_text().splitlines()
            if line.startswith("exit_code=")
        ]
        if len(codes) != 1 or codes[0] != "0":
            raise AnnotationError(
                "upstream_failed: Prokka did not complete successfully"
            )
        manifest = build_bundle(
            args.accession,
            args.gcode,
            args.genome,
            args.faa,
            args.gff,
            args.gbk,
            args.outdir,
        )
        manifest["prokka_log_sha256"] = digest(args.prokka_log)
        write_json(args.outdir / "bundle.json", manifest)
    except (AnnotationError, OSError, UnicodeError) as exc:
        LOGGER.error("%s", exc)
        write_json(
            args.outdir / "bundle.json",
            {
                "schema_version": SCHEMA_VERSION,
                "accession": args.accession,
                "status": "upstream_failed"
                if str(exc).startswith("upstream_failed:")
                else "incompatible_input",
                "reason": str(exc),
            },
        )
        # Emit the failed bundle to allow independent samples to finish; final
        # scientific acceptance checks every requested sample/tool combination.
    return 0


if __name__ == "__main__":
    sys.exit(main())
