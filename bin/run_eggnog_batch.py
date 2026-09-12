#!/usr/bin/env python3
"""Run one native search followed by isolated native annotation of each proteome."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

# Staged helper modules belong to the immutable Nextflow input directory.
sys.dont_write_bytecode = True

from annotation_common import AnnotationError, digest, read_json, write_json
from eggnog_batches import validate_batch
from eggnog_native import (
    EXECUTION_SCHEMA,
    NATIVE_CODE_FILES,
    stage_plan,
    validate_execution,
    write_partitions,
)


def run(batchdir: Path, commands: Path, raw: Path) -> int:
    """Retain each actual native exit status and stop the batch on its first failure."""
    batchdir = batchdir.resolve(strict=True)
    batch, proteins = validate_batch(batchdir)
    code = {name: digest(Path(__file__).with_name(name)) for name in NATIVE_CODE_FILES}
    if batch["search_method"]["native_code"] != code:
        raise AnnotationError(
            "Native eggNOG execution code differs from its planned identity"
        )
    if read_json(commands) != batch["search_method"]["command"]["native_stages"]:
        raise AnnotationError("Native eggNOG phase commands changed after planning")
    plans = stage_plan(batch)
    fasta = Path("task/input.faa")
    if fasta.is_symlink() or digest(fasta) != batch["files"]["input.faa"]:
        raise AnnotationError("Native search input differs from its validated batch")
    if raw != Path("raw") or not raw.is_dir() or (raw / "execution.json").exists():
        raise AnnotationError(
            "Native batch requires the fresh declared raw output directory"
        )
    receipt = {"schema": EXECUTION_SCHEMA, "batch_id": batch["batch_id"], "stages": []}
    write_json(raw / "execution.json", receipt)
    for stage in plans:
        destination = raw / stage["directory"]
        destination.mkdir(parents=True)
        temporary = stage["command"][stage["command"].index("--temp_dir") + 1]
        Path(temporary).mkdir(parents=True)
        started = time.monotonic()
        with (destination / "tool.log").open("w") as handle:
            completed = subprocess.run(
                stage["command"], stdout=handle, stderr=subprocess.STDOUT, check=False
            )
        exit_code = (
            completed.returncode
            if completed.returncode >= 0
            else 128 - completed.returncode
        )
        (destination / "exit_code.txt").write_text(f"{exit_code}\n")
        record = dict(
            stage, exit_code=exit_code, wall_seconds=time.monotonic() - started
        )
        receipt["stages"].append(record)
        write_json(raw / "execution.json", receipt)
        print(json.dumps(record, sort_keys=True), flush=True)
        if exit_code:
            return exit_code
        if stage["name"] == "search":
            write_partitions(batchdir, raw, batch, proteins)
    validate_execution(raw, batch, proteins)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    parser.add_argument("--batch", type=Path, required=True)
    parser.add_argument("--commands", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        return run(args.batch, args.commands, args.output)
    except (OSError, ValueError) as error:
        print(str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
