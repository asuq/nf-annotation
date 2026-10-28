#!/usr/bin/env python3
"""Build the four existing 16S views in one pass using unchanged selection rules."""

from __future__ import annotations

import argparse
from contextlib import ExitStack
from pathlib import Path

from concat_best_16s import (
    Cohort16SError,
    append_fasta_content,
    build_manifest_row,
    load_metadata_atypical_index,
    load_qc_row,
    load_status_row,
    should_include_in_cohort,
    write_manifest,
)

OUTPUTS = {
    "intact": "all_best_16S",
    "partial": "all_partial_16S",
    "low_quality_intact": "low_quality_best_16S",
    "low_quality_partial": "low_quality_partial_16S",
}


def build(
    status_paths: list[Path], qc_directory: Path, metadata: Path | None, output: Path
) -> None:
    """Retain supplied shell-glob order, sample IDs, inclusion rules and output bytes."""
    if len(status_paths) != len(set(status_paths)):
        raise Cohort16SError("Duplicate 16S status input")
    atypical = load_metadata_atypical_index(metadata)
    output.mkdir(parents=True, exist_ok=True)
    manifests = {kind: [] for kind in OUTPUTS}
    observed = set()
    with ExitStack() as stack:
        handles = {
            kind: stack.enter_context(
                (output / (name + ".fna")).open("w", encoding="utf-8")
            )
            for kind, name in OUTPUTS.items()
        }
        for status_path in status_paths:
            if not status_path.name.endswith("_16S_status.tsv"):
                raise Cohort16SError(f"Unexpected 16S status filename: {status_path}")
            internal_id = status_path.name.removesuffix("_16S_status.tsv")
            if not internal_id.strip():
                raise Cohort16SError(f"Missing internal ID for {status_path}")
            fasta = status_path.with_name(internal_id + "_best_16S.fna")
            qc_path = qc_directory / (internal_id + "_checkm2_summary.tsv")
            if not fasta.is_file() or not qc_path.is_file():
                raise Cohort16SError(
                    f"Missing matching FASTA or QC summary for {status_path}"
                )
            status = load_status_row(status_path)
            qc = load_qc_row(qc_path)
            accession = status["accession"]
            if accession in observed or qc["accession"] != accession:
                raise Cohort16SError(
                    f"Duplicate or mismatched 16S/QC accession: {accession}"
                )
            observed.add(accession)
            for kind, handle in handles.items():
                if should_include_in_cohort(status, qc, kind, atypical):
                    append_fasta_content(handle, fasta, internal_id=internal_id.strip())
                    manifests[kind].append(
                        build_manifest_row(status, qc, cohort_kind=kind)
                    )
    for kind, name in OUTPUTS.items():
        write_manifest(output / (name + "_manifest.tsv"), manifests[kind])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--statuses",
        type=Path,
        required=True,
        help="NUL-delimited paths in the original glob order",
    )
    parser.add_argument("--qc-dir", type=Path, required=True)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        paths = [
            Path(value.decode("utf-8"))
            for value in args.statuses.read_bytes().split(b"\0")
            if value
        ]
        build(paths, args.qc_dir, args.metadata, args.output)
    except (Cohort16SError, OSError, UnicodeError) as error:
        parser.exit(1, str(error) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
