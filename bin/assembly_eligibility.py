#!/usr/bin/env python3
"""Attach the downstream size decision to ordered sample records."""

import argparse
from pathlib import Path

from annotation_common import AnnotationError, read_tsv, write_tsv

SHORT = "assembly_below_20000_bp"


def exclusion(bp: int) -> str:
    if bp < 1:
        raise AnnotationError("Assembled length must be positive")
    return SHORT if bp < 20000 else "NA"


def master_fields(sample: dict[str, str], checkm2_low_quality: str) -> dict[str, str]:
    reason = sample["Assembly_Exclusion_Reason"]
    return {
        "Assembly_Bases": sample["Assembly_Bases"],
        "Assembly_Exclusion_Reason": reason,
        "CheckM2_Low_quality": checkm2_low_quality,
        "Low_quality": "true" if reason == SHORT else checkm2_low_quality,
    }


def classify(
    samples: list[dict[str, str]], stats: list[dict[str, str]]
) -> list[dict[str, str]]:
    sizes = {row["accession"]: int(row["genome_size"]) for row in stats}
    if len(sizes) != len(stats) or set(sizes) != {row["accession"] for row in samples}:
        raise AnnotationError("Assembly statistics must match sample membership")
    return [
        dict(
            row,
            Assembly_Bases=str(sizes[row["accession"]]),
            Assembly_Exclusion_Reason=exclusion(sizes[row["accession"]]),
        )
        for row in samples
    ]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("samples", "stats", "output"):
        parser.add_argument(key, type=Path)
    args = parser.parse_args()
    rows = classify(read_tsv(args.samples), read_tsv(args.stats))
    write_tsv(args.output, tuple(rows[0]), rows)
