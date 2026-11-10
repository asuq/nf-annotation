"""Verify published CRISPRCasFinder header restoration from retained native evidence.

Sequence identity is checked by the native adapter at execution and by the
existing publication fingerprints at reuse; this check reads only compact JSON.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

IDENTIFIER_FIELDS = {
    "Id",
    "Name",
    "Version",
    "Gene",
    "Sequence",
    "Seqid",
    "seq_id",
    "contig",
    "Summary_CRISPR-Cas",
}
ALIAS = re.compile(r"(?<![A-Za-z0-9])\bctg[0-9]{8}(?=_|\b)")


def restore(value, tool_id: str, source_id: str, field: str = ""):
    if isinstance(value, dict):
        return {
            key: restore(item, tool_id, source_id, key) for key, item in value.items()
        }
    if isinstance(value, list):
        return [restore(item, tool_id, source_id, field) for item in value]
    if isinstance(value, str) and field in IDENTIFIER_FIELDS:

        def replace(match):
            if match.group() != tool_id:
                raise ValueError("Foreign native contig reference")
            return source_id

        return ALIAS.sub(replace, value)
    return value


def validate_publication(root: Path) -> None:
    mapping = json.loads((root / "contig_mapping.json").read_text())
    native = json.loads((root / "native_result.json").read_text())
    actual = json.loads((root / "result.json").read_text())
    if mapping.get("schema") != "ccfinder-contig-map-v1":
        raise ValueError("Unsupported CRISPRCasFinder contig mapping")
    entries = mapping["contigs"]
    if not entries or len({entry["source_id"] for entry in entries}) != len(entries):
        raise ValueError("Empty or duplicate source contig mapping")
    if [entry["tool_id"] for entry in entries] != [
        f"ctg{i:08d}" for i in range(1, len(entries) + 1)
    ]:
        raise ValueError("Invalid ordered tool contig mapping")
    rows = native["Sequences"]
    by_id = {row["Id"]: row for row in rows}
    if len(by_id) != len(rows) or set(by_id) != {entry["tool_id"] for entry in entries}:
        raise ValueError("Native contig membership differs from mapping")
    expected = dict(native, Sequences=[])
    for entry in entries:
        row = by_id[entry["tool_id"]]
        if row["Length"] != entry["length"] or not all(
            isinstance(row[key], list) for key in ("Crisprs", "Cas")
        ):
            raise ValueError("Native contig length or CRISPR/Cas structure differs")
        converted = restore(row, entry["tool_id"], entry["source_id"])
        converted.update(
            Id=entry["source_id"],
            Version=entry["source_id"],
            Description=entry["description"],
        )
        expected["Sequences"].append(converted)
    if actual != expected:
        raise ValueError("Published CRISPRCasFinder result differs from native mapping")
