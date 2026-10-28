#!/usr/bin/env python3
"""Content-bound, code-bound source validation for cumulative ANI cohorts."""

from __future__ import annotations

import argparse
import hashlib
import os
import stat
import time
from pathlib import Path

from annotation_common import AnnotationError, digest, identity, read_json, write_json
from annotation_source import validate_source
from publication_references import REFERENCE_CONTRACT, REFERENCE_FILE, Publication


def code_id() -> str:
    """Identify every Python helper available to the validation contract."""
    result = hashlib.sha256()
    for path in sorted(Path(__file__).resolve().parent.glob("*.py")):
        result.update(path.name.encode() + b"\0")
        result.update(path.read_bytes())
        result.update(b"\0")
    return result.hexdigest()


def fingerprint(
    root: Path, *, cache: dict | None = None, active: tuple[Path, ...] = ()
) -> dict:
    """Hash scientific publication contents, excluding mutable execution logs."""
    root = root.resolve(strict=True)
    cache = {} if cache is None else cache
    if root in active:
        raise AnnotationError(f"Circular publication reference: {root}")
    if root in cache:
        return cache[root]
    manifest = root / "annotation_results.json"
    if manifest.is_symlink() or not manifest.is_file():
        raise AnnotationError(f"Invalid publication manifest: {manifest}")
    # DirEntry retains directory type/stat information. Avoid repeated is_file,
    # is_dir and stat round trips for every small file on a shared filesystem.
    paths = [(manifest, manifest.stat())]
    reference = read_json(manifest).get("pipeline_contract") == REFERENCE_CONTRACT
    publication = Publication(root) if reference else None
    dependencies = []
    if publication:
        paths.append((root / REFERENCE_FILE, (root / REFERENCE_FILE).stat()))
        for source in publication.document["sources"]:
            current = fingerprint(
                Path(source["root"]), cache=cache, active=(*active, root)
            )
            if current["sha256"] != source["sha256"]:
                raise AnnotationError(f"Referenced source changed: {source['root']}")
            dependencies.append({"root": current["root"], "sha256": current["sha256"]})
    elif (root / REFERENCE_FILE).exists():
        raise AnnotationError(
            "Reference manifest lacks the reference publication contract"
        )
    links = []
    pending = []
    for name in ("tables", "samples", "annotation_batches"):
        directory = root / name
        if name == "annotation_batches" and not directory.exists():
            continue
        if directory.is_symlink() or not directory.is_dir():
            raise AnnotationError(f"Invalid publication directory: {directory}")
        pending.append(directory)
    while pending:
        with os.scandir(pending.pop()) as entries:
            for entry in entries:
                if entry.is_symlink():
                    relative = Path(entry.path).relative_to(root).as_posix()
                    if publication is None or relative not in publication.links:
                        raise AnnotationError(
                            f"Publication fingerprint cannot follow a symlink: {entry.path}"
                        )
                    links.append([relative, "reference", os.readlink(entry.path)])
                    continue
                if entry.is_dir(follow_symlinks=False):
                    pending.append(Path(entry.path))
                else:
                    info = entry.stat(follow_symlinks=False)
                    if not stat.S_ISREG(info.st_mode):
                        raise AnnotationError(f"Special publication file: {entry.path}")
                    paths.append((Path(entry.path), info))
    tree = hashlib.sha256()
    files = total = 0
    for path, before in sorted(paths, key=lambda item: item[0]):
        checksum = digest(path)
        after = path.stat()
        if (before.st_ino, before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
            after.st_ctime_ns,
        ):
            raise AnnotationError(f"Publication changed during fingerprinting: {path}")
        tree.update(
            identity(
                [path.relative_to(root).as_posix(), before.st_size, checksum]
            ).encode()
        )
        files += 1
        total += before.st_size
    for link in sorted(links):
        tree.update(identity(link).encode())
    result = {
        "root": str(root),
        "sha256": tree.hexdigest(),
        "files": files,
        "bytes": total,
    }
    if dependencies:
        result["references"] = dependencies
    cache[root] = result
    return result


def validate(roots: list[Path], fingerprints: list[dict]) -> list[dict]:
    """Validate one bounded source group after fresh content fingerprinting."""
    if len(roots) != len(fingerprints) or not roots:
        raise AnnotationError("Source validation roots and fingerprints must match")
    validator = code_id()
    receipts = []
    for root, record in zip(roots, fingerprints, strict=True):
        if str(root.resolve(strict=True)) != record["root"]:
            raise AnnotationError("Source fingerprint belongs to a different root")
        if (
            read_json(root / "annotation_results.json").get("pipeline_contract")
            == REFERENCE_CONTRACT
        ):
            validate_source(root, _content_checked=True)
        else:
            validate_source(root)
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
