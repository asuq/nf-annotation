"""Explicit reviewed genetic-code decisions, keyed by accession and FASTA identity."""

from __future__ import annotations

import base64
import gzip
import hashlib
import json
import re
from pathlib import Path

FIELDS = {"accession", "fasta_identity_md5", "genetic_code", "review_id", "reason"}


def encode_candidates(entries: list[dict]) -> str:
    """Keep internal TSV fields safe for Nextflow's lightweight CSV parser."""
    return base64.b64encode(
        json.dumps(entries, sort_keys=True, separators=(",", ":")).encode()
    ).decode("ascii")


def recorded_candidates(sample: dict) -> list[dict]:
    return json.loads(
        base64.b64decode(sample.get("review_candidates_base64", "W10="), validate=True)
    )


def validate_entry(entry: dict) -> dict:
    if not isinstance(entry, dict) or set(entry) != FIELDS:
        raise ValueError(
            "Reviewed code entries require accession, fasta_identity_md5, genetic_code, review_id and reason"
        )
    if any(
        not isinstance(entry[k], str)
        or not entry[k].strip()
        or any(ord(c) < 32 for c in entry[k])
        for k in ("accession", "review_id", "reason")
    ):
        raise ValueError(
            "Reviewed code identity and review provenance must be nonempty text"
        )
    if not isinstance(entry["fasta_identity_md5"], str) or not re.fullmatch(
        "[0-9a-f]{32}", entry["fasta_identity_md5"]
    ):
        raise ValueError("Reviewed FASTA identity must be a lowercase MD5 digest")
    if type(entry["genetic_code"]) is not int or entry["genetic_code"] not in (4, 11):
        raise ValueError("Reviewed genetic code must be 4 or 11")
    return entry


def load_manifest(path: Path | None) -> dict[str, list[dict]]:
    if path is None:
        return {}
    record = json.loads(path.read_text())
    if (
        not isinstance(record, dict)
        or set(record) != {"schema", "reviews"}
        or record["schema"] != "reviewed-genetic-codes-v1"
        or not isinstance(record["reviews"], list)
    ):
        raise ValueError("Unsupported reviewed genetic-code manifest")
    result, seen = {}, set()
    for entry in record["reviews"]:
        validate_entry(entry)
        key = entry["accession"], entry["fasta_identity_md5"]
        if key in seen:
            raise ValueError(
                "Duplicate accession/FASTA identity in reviewed code manifest"
            )
        seen.add(key)
        result.setdefault(entry["accession"], []).append(entry)
    return {
        key: sorted(values, key=lambda x: x["fasta_identity_md5"])
        for key, values in result.items()
    }


def fasta_identity(path: Path) -> str:
    """MD5 of ordered IDs and case-sensitive DNA, independent of FASTA wrapping."""
    with path.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    checksum = hashlib.md5(usedforsecurity=False)
    seen, length = set(), 0
    with gzip.open(path, "rb") if compressed else path.open("rb") as stream:
        for raw in stream:
            line = raw.rstrip(b"\r\n")
            if not line:
                continue
            if line.startswith(b">"):
                fields = line[1:].split()
                if not fields or fields[0] in seen or (seen and not length):
                    raise ValueError(
                        "Empty or duplicate FASTA record in reviewed genome"
                    )
                seen.add(fields[0])
                length = 0
                checksum.update(b"\0record\0" + fields[0] + b"\0")
            else:
                if not seen or not re.fullmatch(
                    rb"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", line
                ):
                    raise ValueError(
                        "Invalid DNA or sequence whitespace in reviewed genome"
                    )
                checksum.update(line)
                length += len(line)
    if not seen or not length:
        raise ValueError("Empty reviewed genome or final contig")
    return checksum.hexdigest()


def select_review(
    accession: str, genome: Path | None, candidates: list[dict]
) -> dict | None:
    if not isinstance(candidates, list):
        raise ValueError("Reviewed code candidates must be a list")  # noqa: TRY004
    if not candidates:
        return None
    for entry in candidates:
        validate_entry(entry)
        if entry["accession"] != accession:
            raise ValueError("Foreign reviewed code accession")
    if genome is None:
        raise ValueError("A reviewed code requires the current genome FASTA")
    fingerprint = fasta_identity(genome)
    matches = [
        entry for entry in candidates if entry["fasta_identity_md5"] == fingerprint
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Reviewed FASTA identity is missing or ambiguous for {accession}"
        )
    return matches[0]


def validate_published_decision(sample: dict, native: Path, master: dict) -> None:
    """Replay an explicit review from its retained paired reports and genome."""
    candidates = recorded_candidates(sample)
    if not candidates and master.get("Gcode_Review", "NA") == "NA":
        return
    import summarise_checkm2
    from master_table_contract import GCODE_PROVENANCE_COLUMNS

    expected = summarise_checkm2.build_output_row(
        sample["accession"],
        summarise_checkm2.parse_report(native / "checkm2_gcode4/quality_report.tsv"),
        summarise_checkm2.parse_report(native / "checkm2_gcode11/quality_report.tsv"),
        warnings=[],
        gcode_rule=master["Gcode_Rule"],
        genome=native / f"staged/{sample['internal_id']}.fasta",
        genome_id=sample["internal_id"],
        review_candidates=candidates,
    )
    for key in ("Gcode", "Low_quality", *GCODE_PROVENANCE_COLUMNS):
        if master.get(key) != expected[key]:
            raise ValueError(
                f"Published reviewed genetic-code decision disagrees on {key}"
            )
