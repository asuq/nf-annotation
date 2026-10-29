#!/usr/bin/env python3
"""Recheck immutable reference dependencies before accepting a completed run."""

import argparse
from pathlib import Path

from annotation_common import AnnotationError, read_json
from fingerprinting.source_fingerprint import FingerprintScanner
from publication_references import Publication


def verify(
    publication: Path,
    acceptance: Path,
    snapshot: Path | None = None,
    *,
    workers: int = 1,
    progress: Path | None = None,
    source_validation: Path | None = None,
) -> None:
    status = read_json(acceptance)
    if status.get("complete") is not True:
        raise AnnotationError("Requested annotation analyses did not all succeed")
    with FingerprintScanner(workers, progress=progress) as scanner:
        if snapshot is not None:
            expected = read_json(snapshot)
            current = scanner.fingerprint(publication)
            if (
                current["root"] != expected["root"]
                or current["sha256"] != expected["sha256"]
            ):
                raise AnnotationError("Reference input changed during annotation")
            scanner.metrics["completed_roots"] = 1
        else:
            view = Publication(publication, pending=True)
            scanner.metrics["requested_roots"] = len(view.snapshots)
            for source in view.snapshots:
                if (
                    scanner.fingerprint(Path(source["root"]))["sha256"]
                    != source["sha256"]
                ):
                    raise AnnotationError(
                        f"Referenced source changed during analysis: {source['root']}"
                    )
                scanner.metrics["completed_roots"] += 1
                scanner.report(force=True)
            view.check_links()
            if source_validation is not None:
                expected_sources = read_json(source_validation)["fingerprints"]
                scanner.metrics["requested_roots"] += len(expected_sources)
                for expected in expected_sources:
                    current = scanner.fingerprint(Path(expected["root"]))
                    if current["sha256"] != expected["sha256"]:
                        raise AnnotationError(
                            f"Reference input changed during analysis: {expected['root']}"
                        )
                    scanner.metrics["completed_roots"] += 1
                    scanner.report(force=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--publication", type=Path, required=True)
    parser.add_argument("--acceptance", type=Path, required=True)
    parser.add_argument("--snapshot", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--progress", type=Path)
    parser.add_argument("--source-validation", type=Path)
    args = parser.parse_args()
    try:
        verify(
            args.publication.resolve(strict=True),
            args.acceptance,
            args.snapshot,
            workers=args.workers,
            progress=args.progress,
            source_validation=args.source_validation,
        )
        args.output.write_text(
            "Reference dependencies and requested analyses passed validation.\n"
        )
    except (AnnotationError, OSError, ValueError) as error:
        parser.exit(1, str(error) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
