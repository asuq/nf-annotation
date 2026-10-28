#!/usr/bin/env python3
"""Publish prepared cohort tables only after checking the existing run identity."""

import argparse
import json
import shutil
from pathlib import Path


def check(prepared: Path, previous: Path | None, output: Path) -> None:
    current = json.loads((prepared / "cohort_update_run.json").read_text())
    if previous is not None and json.loads(previous.read_text()) != current:
        raise ValueError(
            "The existing update output belongs to different inputs or settings; choose a fresh --outdir."
        )
    output.mkdir(parents=True, exist_ok=True)
    for path in prepared.iterdir():
        if not path.is_file() or path.is_symlink():
            raise ValueError("Prepared cohort tables must be ordinary files")
        shutil.copyfile(path, output / path.name)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--previous-update", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        check(args.prepared, args.previous_update, args.output)
    except (OSError, ValueError) as error:
        parser.exit(1, str(error) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
