#!/usr/bin/env python3
"""Read-only fingerprint timings after creating synthetic local I/O fixtures."""

import argparse
import importlib.util
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
from fingerprinting.source_fingerprint import fingerprint_sources

BASELINE = "6e60ca05c963a1d5ab71042dd4680e277f0df245"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", required=True, type=Path)
    args = parser.parse_args()
    args.directory.mkdir()
    roots = []
    for index in range(4):
        root = args.directory / f"source-{index}"
        (root / "tables").mkdir(parents=True)
        (root / "samples/A").mkdir(parents=True)
        (root / "annotation_results.json").write_text("{}\n")
        block = bytes(range(256)) * 4096
        for number in range(8):
            with (root / f"samples/A/large-{number}.bin").open("wb") as handle:
                for _ in range(32):
                    handle.write(block)
        for number in range(1024):
            (root / f"tables/small-{number}.tsv").write_bytes(block[:512])
        roots.append(root)
    old_path = args.directory / "baseline.py"
    old_path.write_bytes(
        subprocess.check_output(
            ["git", "show", f"{BASELINE}:bin/reuse_source_validation.py"], cwd=ROOT
        )
    )
    spec = importlib.util.spec_from_file_location("baseline_fingerprint", old_path)
    baseline = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(baseline)
    reports, expected = [], None
    for workers in (0, 1, 4, 8, 8, 4, 1, 0):
        start = time.monotonic()
        if workers == 0:
            cache = {}
            records = [baseline.fingerprint(root, cache=cache) for root in roots]
            metrics = {"workers": 1, "implementation": "v0.4.5 serial"}
        else:
            records, metrics = fingerprint_sources(roots, workers=workers)
            metrics["implementation"] = "bounded readers"
        elapsed = time.monotonic() - start
        if expected is None:
            expected = records
        if records != expected:
            raise AssertionError("Fingerprint bytes differ from v0.4.5")
        reports.append({**metrics, "seconds": elapsed})
    result = {
        "baseline": BASELINE,
        "platform": platform.platform(),
        "python": platform.python_version(),
        "files": sum(row["files"] for row in expected),
        "bytes": sum(row["bytes"] for row in expected),
        "identical_records": True,
        "runs": reports,
        "limits": "Synthetic local filesystem; caches were not flushed; recently written inputs and repeated reads may be cached. No shared-storage or native scientific benchmark.",
    }
    (args.directory / "benchmark.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
