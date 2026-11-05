#!/usr/bin/env python3
"""Engineering timings for native profile reuse and condensed AAI clustering."""

from __future__ import annotations
import argparse
import gzip
import json
import resource
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from annotation_common import digest, read_json, write_json, write_tsv


def rss_bytes(usage):
    return usage.ru_maxrss if sys.platform == "darwin" else usage.ru_maxrss * 1024


def native_cache(fixtures, output):
    import prepare_aai
    from types import SimpleNamespace

    output.mkdir(parents=True)
    bundle = output / "bundle"
    bundle.mkdir()
    with gzip.open(fixtures / "mycoplasma.faa.gz", "rb") as src:
        (bundle / "proteins.faa").write_bytes(src.read())
    write_json(
        bundle / "bundle.json",
        dict(
            accession="GCF_000027325.1",
            status="success",
            files={"proteins.faa": digest(bundle / "proteins.faa")},
        ),
    )
    acc = "GCF_000027325.1"
    args = SimpleNamespace(
        **{
            k: output / (k + ".json")
            for k in (
                "records",
                "metadata",
                "checkm2",
                "busco",
                "stats",
                "bundles",
                "sources",
                "output",
            )
        },
        primary_busco_column="BUSCO_fixture",
        domain="bacteria",
        threads=2,
        destination=output / "publication",
    )
    args.destination.mkdir()
    write_json(
        args.records, [dict(accession=acc, is_new="false", assembly_level="Scaffold")]
    )
    write_tsv(
        args.metadata,
        ("accession", "Assembly_Level"),
        [dict(accession=acc, Assembly_Level="Scaffold")],
    )
    write_tsv(
        args.checkm2,
        (
            "accession",
            "Gcode",
            "Low_quality",
            "Completeness_gcode4",
            "Contamination_gcode4",
        ),
        [
            dict(
                accession=acc,
                Gcode="4",
                Low_quality="false",
                Completeness_gcode4="95",
                Contamination_gcode4="1",
            )
        ],
    )
    write_tsv(
        args.stats,
        ("accession", "n50", "scaffolds", "genome_size", "gc_content"),
        [
            dict(
                accession=acc,
                n50="1000",
                scaffolds="1",
                genome_size="1000",
                gc_content="30",
            )
        ],
    )
    busco = output / "busco.tsv"
    write_tsv(
        busco, ("accession", "BUSCO_fixture"), [dict(accession=acc, BUSCO_fixture="NA")]
    )
    write_json(args.busco, [str(busco)])
    write_json(args.bundles, [str(bundle)])
    write_json(args.sources, [])
    runs = []
    for name in ("first", "reuse"):
        prepare_aai.prepare(args)
        metrics = read_json(output / "preprocessing_metrics.json")
        metrics["run"] = name
        runs.append(metrics)
    assert (
        runs[0]["extracted"] == 1
        and runs[1]["extracted"] == 0
        and runs[1]["reused"] == 1
    )
    result = dict(
        native_protein_input=True,
        synthetic_QC_metadata=True,
        runs=runs,
        input_protein_bytes=(bundle / "proteins.faa").stat().st_size,
        profile_bytes=sum(
            p.stat().st_size
            for p in (args.destination / "cohort/fastaai/profiles").rglob("*")
            if p.is_file()
        ),
        peak_worker_rss_bytes=rss_bytes(resource.getrusage(resource.RUSAGE_CHILDREN)),
    )
    write_json(output / "benchmark.json", result)
    print(json.dumps(result, indent=2))


def synthetic_clustering(samples, output):
    from cluster_aai import cluster_domain
    from test_aai_clustering import record

    output.mkdir(parents=True)
    from aai_evidence import create

    pairs = output / "evidence"
    records = [record(f"SYN{i:06d}") for i in range(samples)]
    arrays = create(pairs, [r["accession"] for r in records])
    arrays["jaccard"][:] = 0.55
    arrays["shared"][:] = 60
    arrays["possible"][:] = 80
    for array in arrays.values():
        if hasattr(array, "flush"):
            array.flush()
            array._mmap.close()
    del array, arrays
    start = time.monotonic()
    rows, reps, members, _ = cluster_domain(records, pairs, 0.65, "default", "Bacteria")
    elapsed = time.monotonic() - start
    assert len(rows) == samples and len(reps) == 1 and len(members) == samples
    result = dict(
        synthetic=True,
        samples=samples,
        pairs=samples * (samples - 1) // 2,
        seconds=elapsed,
        pair_file_bytes=sum(p.stat().st_size for p in pairs.iterdir()),
        condensed_vector_bytes=samples * (samples - 1) // 2 * 8,
        process_peak_rss_bytes=rss_bytes(resource.getrusage(resource.RUSAGE_SELF)),
    )
    write_json(output / "benchmark.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fixtures", type=Path)
    p.add_argument("--samples", type=int, default=2500)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.fixtures:
        native_cache(a.fixtures, a.output)
    else:
        synthetic_clustering(a.samples, a.output)
