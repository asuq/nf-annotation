#!/usr/bin/env python3
"""Native parity against unmodified pinned FastAAI on frozen public proteomes."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import prepare_aai
from annotation_common import digest, protein_query_id, write_json
from compare_aai import compare

EXPECTED = {
    "ecoli.faa.gz": "283c91ff66d0597f2488fd1c00ad40f8d0a0b37a92f981b6148d75fa4a934896",
    "mycoplasma.faa.gz": "7117e926bed866db7fab6d53567d6de4fe55845deb2dc995a86720850e808c58",
    "archaea.faa.gz": "cbd09de5edd1193e81c0e93f290e8a9b0d992055817512007e68b3781514d58a",
}


def run(fixtures, output):
    import fastaai.fastaai as native  # ty: ignore[unresolved-import]  # Provided by the pinned native image.

    output.mkdir(parents=True)
    for name, sha in EXPECTED.items():
        if digest(fixtures / name) != sha:
            raise AssertionError("Frozen proteome changed: " + name)
    upstream = Path(native.__file__).with_suffix(".upstream.py")
    original = upstream.read_text()
    adapted = Path(native.__file__).read_text()
    import ast

    class WithoutBinaryWriter(ast.NodeTransformer):
        def visit_If(self, node):
            if ast.unparse(node.test) == "style == 'binary'":
                return self.visit(node.orelse[0])
            return self.generic_visit(node)

    assert ast.dump(WithoutBinaryWriter().visit(ast.parse(adapted))) == ast.dump(
        ast.parse(original)
    )
    counts = {}
    times = {}
    profiles = []
    records: list[dict[str, Any]] = []
    method = json.loads(prepare_aai.RUNTIME.read_text())
    inputs = output / "inputs"
    inputs.mkdir()
    canonical = output / "canonical"
    canonical.mkdir()
    members = [
        ("ecoli.faa.gz", "Bacteria", 11),
        ("mycoplasma.faa.gz", "Bacteria", 4),
        ("archaea.faa.gz", "Archaea", 11),
        ("archaea.faa.gz", "Archaea", 11),
    ]

    def forbidden(*args, **kwargs):
        raise AssertionError("Gene prediction was invoked")

    native.input_file.genome_to_protein = forbidden
    prepare_aai.native_init()
    for i, (name, domain, code) in enumerate(members):
        alias = f"S{i:08d}"
        protein = inputs / (alias + ".faa.gz")
        protein.symlink_to((fixtures / name).resolve())
        readable = canonical / (alias + ".faa")
        accession = {
            "ecoli.faa.gz": "GCF_000005845.2",
            "mycoplasma.faa.gz": "GCF_000027325.1",
            "archaea.faa.gz": "GCF_000091665.1",
        }[name]
        with gzip.open(protein, "rt") as src, readable.open("w") as dst:
            ordinal = 0
            for line in src:
                if line.startswith(">"):
                    ordinal += 1
                    line = ">" + protein_query_id(accession, ordinal) + "\n"
                dst.write(line)
        work = output / ("profile-" + alias)
        work.mkdir()
        start = time.monotonic()
        count, crystal = prepare_aai.extract((readable, work / "native"))
        times[alias] = time.monotonic() - start
        assert count > 0
        for filename in ("markers.faa", "markers.hmm"):
            shutil.copyfile(crystal / filename, work / filename)
        write_json(
            work / "profile.json", {"profile_id": work.name, "marker_count": count}
        )
        profiles.append(work)
        counts[alias] = count
        records.append(
            {
                "accession": alias,
                "domain": domain,
                "status": "ready",
                "profile": str(work),
                "marker_count": count,
                "gcode": code,
            }
        )
    short = inputs / "empty.faa"
    short.write_text(">short\nMWA\n")
    count, crystal = prepare_aai.extract((short, output / "no-markers"))
    assert count == 0 and crystal is None
    short.unlink()
    baseline = output / "upstream"
    commands = [
        [
            sys.executable,
            str(upstream),
            "build_db",
            "--proteins",
            str(inputs),
            "--output",
            str(baseline),
            "--database",
            "baseline.db",
            "--threads",
            "2",
        ],
        [
            sys.executable,
            str(upstream),
            "db_query",
            "--query",
            str(baseline / "database/baseline.db"),
            "--target",
            str(baseline / "database/baseline.db"),
            "--output",
            str(baseline),
            "--threads",
            "2",
            "--output_style",
            "tsv",
        ],
    ]
    for command in commands:
        with (output / "upstream.log").open("a") as log:
            subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    expected = {}
    for path in (baseline / "results").glob("*_results.txt"):
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                expected[
                    row["query"].removesuffix(".faa.gz"),
                    row["target"].removesuffix(".faa.gz"),
                ] = row
    compared = 0
    for domain in ("Bacteria", "Archaea", "NumericalParityControl"):
        selected = (
            records
            if domain == "NumericalParityControl"
            else [r for r in records if r["domain"] == domain]
        )
        destination = output / domain
        start = time.monotonic()
        compare(selected, [Path(r["profile"]) for r in selected], destination, 2)
        times[domain + "_compare"] = time.monotonic() - start
        from aai_evidence import iter_pairs

        for row in iter_pairs(destination):
            ref = expected[row["query"], row["target"]]
            assert round(float(row["avg_jacc_sim"]), 4) == float(ref["avg_jacc_sim"]), (
                row,
                ref,
            )
            assert int(row["num_shared_SCPs"]) == int(ref["num_shared_SCPs"])
            assert int(row["poss_shared_SCPs"]) == int(ref["poss_shared_SCPs"])
            if ref["AAI_estimate"].startswith(("<", ">")):
                assert row["AAI_estimate"] == ref["AAI_estimate"]
            else:
                assert (
                    abs(float(row["AAI_estimate"]) - float(ref["AAI_estimate"]))
                    <= 0.0050001
                )
            compared += 1
    report = {
        "runtime": method,
        "proteome_sha256": EXPECTED,
        "genetic_codes": [4, 11],
        "domains": ["Bacteria", "Archaea"],
        "markers": counts,
        "native_pairs_compared": compared,
        "no_markers_passed": True,
        "gene_prediction_called": False,
        "readable_ids_match_original": True,
        "seconds": times,
    }
    write_json(output / "qualification.json", report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fixtures", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    run(a.fixtures, a.output)
