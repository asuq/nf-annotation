#!/usr/bin/env python3
"""Seal explicit dependency bindings after all reference-import batches finish."""

import argparse
from pathlib import Path

from annotation_common import AnnotationError, identity, read_json, read_tsv, write_json
from annotation_path_lists import read_path_list
from publication_references import REFERENCE_FILE, Publication


def complete(audits: list[Path], samples: Path, destination: Path) -> dict:
    accessions = [row["accession"] for row in read_tsv(samples, ["accession"])]
    if not accessions or len(accessions) != len(set(accessions)):
        raise AnnotationError("Reference publication requires unique requested samples")
    expected = set(accessions)
    observed = {}
    snapshots = {}
    for path in audits:
        audit = read_json(path)
        if audit.get("complete") is not True:
            raise AnnotationError("Cannot publish an incomplete reference import")
        for row in audit["samples"]:
            accession = row["accession"]
            if (
                accession in observed
                or accession not in expected
                or row["status"] != "done"
            ):
                raise AnnotationError("Duplicate, foreign or failed reference sample")
            reference = row["reference"]
            observed[accession] = reference
            for source in (reference["native"], reference["bundle"]):
                if source is None or source == {"local": True}:
                    continue
                root = source["root"]
                if root in snapshots and snapshots[root] != source:
                    raise AnnotationError("Conflicting reference source snapshots")
                snapshots[root] = source
    if set(observed) != expected:
        raise AnnotationError("Reference import omitted requested samples")
    roots = sorted(snapshots)
    indexes = {root: index for index, root in enumerate(roots)}
    records = {}
    for accession in accessions:
        reference = observed[accession]
        bundle = reference["bundle"]
        records[accession] = {
            "native_source": indexes[reference["native"]["root"]],
            "entries": reference["entries"],
            "bundle": bundle
            if bundle is None or bundle == {"local": True}
            else {"source": indexes[bundle["root"]]},
        }
    record = {
        "schema_version": 1,
        "mode": "reference",
        "sources": [snapshots[root] for root in roots],
        "samples": records,
    }
    record["reference_id"] = identity(record)
    target = destination / REFERENCE_FILE
    if target.exists() and read_json(target) != record:
        raise AnnotationError(
            "Existing reference publication has different dependencies; use a fresh --outdir"
        )
    Publication(destination, pending=True, document=record)
    temporary = destination / (REFERENCE_FILE + ".tmp")
    if temporary.exists() or temporary.is_symlink():
        raise AnnotationError("Unexpected reference publication temporary file")
    write_json(temporary, record)
    temporary.replace(target)
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audits", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    args = parser.parse_args()
    try:
        record = complete(
            read_path_list(args.audits),
            args.samples,
            args.destination.resolve(strict=True),
        )
        write_json(args.receipt, {"reference_id": record["reference_id"]})
    except (AnnotationError, OSError, ValueError, KeyError) as error:
        parser.exit(1, str(error) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
