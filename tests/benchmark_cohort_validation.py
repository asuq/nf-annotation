#!/usr/bin/env python3
"""Controlled local validation benchmark; no native scientific jobs or production reads."""

import copy
import importlib.util
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import prepare_cohort_update
from test_prepare_cohort_update import PrepareCohortUpdateTestCase


def main():
    destination = Path(sys.argv[1])
    case = PrepareCohortUpdateTestCase()
    case.setUp()
    try:
        accessions = [f"SYN{i:04d}" for i in range(64)]
        case.published_cohort(accessions)
        sequence = ">contig1 synthetic storage fixture\n" + ("ACGT" * 20 + "\n") * 62500
        for accession in accessions:
            case.genomes[accession].write_text(sequence)
            (case.source / f"samples/{accession}/staged/{accession}.fasta").write_text(
                sequence
            )
        args = case.arguments(list(reversed(accessions)))
        baseline = case.root / "baseline.py"
        baseline.write_bytes(
            subprocess.check_output(
                [
                    "git",
                    "show",
                    "b1278586a745f20a4f754f3251eac4810dd38417:bin/prepare_cohort_update.py",
                ],
                cwd=ROOT,
            )
        )
        spec = importlib.util.spec_from_file_location("baseline_cohort", baseline)
        old = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(old)
        expected = None
        runs = []
        for index, workers in enumerate((0, 1, 4, 8, 4, 1, 0)):
            current = copy.deepcopy(args)
            current.outdir = case.root / f"output-{index}"
            current.workers = workers or 1
            start = time.monotonic()
            (old if workers == 0 else prepare_cohort_update).run_prepare(current)
            elapsed = time.monotonic() - start
            outputs = {p.name: p.read_bytes() for p in current.outdir.iterdir()}
            if expected is None:
                expected = outputs
            if outputs != expected:
                raise AssertionError("Scientific tables or run identity changed")
            runs.append(
                {
                    "implementation": "v0.4.6" if workers == 0 else "candidate",
                    "workers": workers or 1,
                    "seconds": elapsed,
                }
            )
        destination.write_text(
            json.dumps(
                {
                    "samples": len(accessions),
                    "bases_per_genome": 5_000_000,
                    "platform": platform.platform(),
                    "python": platform.python_version(),
                    "runs": runs,
                    "identical_outputs": True,
                    "limits": "Synthetic local fixtures; caches not flushed and host load uncontrolled. Includes source validation and sample checks; no shared-storage or native scientific qualification.",
                },
                indent=2,
            )
            + "\n"
        )
        print(destination.read_text())
    finally:
        case.doCleanups()


if __name__ == "__main__":
    main()
