#!/usr/bin/env python3
"""Content-bound, code-bound source validation for cumulative ANI cohorts."""

from __future__ import annotations

import argparse
import hashlib
import time
from pathlib import Path

from annotation_common import AnnotationError, identity, read_json, write_json
from annotation_source import validate_source
from fingerprinting.source_fingerprint import fingerprint


def code_id() -> str:
    """Identify every Python helper available to the validation contract."""
    result = hashlib.sha256()
    root = Path(__file__).resolve().parent
    for path in sorted(
        root.rglob("*.py"), key=lambda item: item.relative_to(root).as_posix()
    ):
        result.update(path.relative_to(root).as_posix().encode() + b"\0")
        result.update(path.read_bytes())
        result.update(b"\0")
    return result.hexdigest()


def validate(roots: list[Path], fingerprints: list[dict]) -> list[dict]:
    """Validate one bounded source group after fresh content fingerprinting."""
    if len(roots) != len(fingerprints) or not roots:
        raise AnnotationError("Source validation roots and fingerprints must match")
    validator = code_id()
    receipts = []
    validated = {}
    for root, record in zip(roots, fingerprints, strict=True):
        if str(root.resolve(strict=True)) != record["root"]:
            raise AnnotationError("Source fingerprint belongs to a different root")
        validate_source(root, _content_checked=True, _validated=validated)
        receipt = {"source": record, "validator": validator}
        receipt["receipt_id"] = identity(receipt)
        receipts.append(receipt)
    return receipts


def checked_roots(receipts: list[dict], expected: list[dict]) -> set[Path]:
    """Accept only internally generated receipts for the current bytes and code."""
    validator = code_id()
    wanted = {record["root"]: record for record in expected}
    if len(wanted) != len(expected) or len(receipts) != len(expected):
        raise AnnotationError("Duplicate or incomplete validated source receipts")
    observed = set()
    for receipt in receipts:
        record = receipt.get("source", {})
        root = record.get("root")
        if (
            root in observed
            or record != wanted.get(root)
            or receipt.get("validator") != validator
            or receipt.get("receipt_id")
            != identity(
                {key: value for key, value in receipt.items() if key != "receipt_id"}
            )
        ):
            raise AnnotationError(
                "Source validation receipt does not match current content and code"
            )
        observed.add(root)
    if observed != set(wanted):
        raise AnnotationError("Missing validated source receipt")
    return {Path(root) for root in observed}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation", choices=("fingerprint", "validate"))
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--fingerprints", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metrics", type=Path)
    args = parser.parse_args()
    start = time.monotonic()
    try:
        roots = [Path(value) for value in read_json(args.sources)]
        cache = {}
        records = (
            [fingerprint(root, cache=cache) for root in roots]
            if args.operation == "fingerprint"
            else validate(roots, read_json(args.fingerprints))
        )
        write_json(args.output, records)
        if args.metrics:
            write_json(
                args.metrics,
                {"operation": args.operation, "seconds": time.monotonic() - start},
            )
    except (AnnotationError, OSError, ValueError, KeyError) as error:
        parser.exit(1, str(error) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
