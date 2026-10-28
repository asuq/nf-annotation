#!/usr/bin/env python3
"""Recheck immutable reference dependencies before accepting a completed run."""

import argparse
from pathlib import Path

from annotation_common import AnnotationError, read_json
from publication_references import Publication
from reuse_source_validation import fingerprint


def verify(publication: Path, acceptance: Path, snapshot: Path | None = None) -> None:
    status = read_json(acceptance)
    if status.get("complete") is not True:
        raise AnnotationError("Requested annotation analyses did not all succeed")
    if snapshot is not None:
        expected = read_json(snapshot)
        current = fingerprint(publication)
        if (
            current["root"] != expected["root"]
            or current["sha256"] != expected["sha256"]
        ):
            raise AnnotationError("Reference input changed during annotation")
    else:
        view = Publication(publication, pending=True)
        cache = {}
        for source in view.document["sources"]:
            if (
                fingerprint(Path(source["root"]), cache=cache)["sha256"]
                != source["sha256"]
            ):
                raise AnnotationError(
                    f"Referenced source changed during analysis: {source['root']}"
                )
        view.check_links()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--publication", type=Path, required=True)
    parser.add_argument("--acceptance", type=Path, required=True)
    parser.add_argument("--snapshot", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        verify(args.publication.resolve(strict=True), args.acceptance, args.snapshot)
        args.output.write_text(
            "Reference dependencies and requested analyses passed validation.\n"
        )
    except (AnnotationError, OSError, ValueError) as error:
        parser.exit(1, str(error) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
