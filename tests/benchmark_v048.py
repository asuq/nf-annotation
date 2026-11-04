#!/usr/bin/env python3
"""Compare v0.4.7/v0.4.8 reporting with valid synthetic multi-gene bundles."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
selected = (
    Path(sys.argv[sys.argv.index("--source") + 1]) if "--worker" in sys.argv else ROOT
)
sys.path.insert(0, str(selected / "bin"))
from annotation_common import TOOLS, digest, identity, read_json, write_json, write_tsv


def generate(root: Path, samples: int, genes: int) -> None:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import SeqFeature, SimpleLocation
    from Bio.SeqRecord import SeqRecord
    from prepare_annotation_bundle import build_bundle

    root.mkdir(parents=True)
    protein = "M" + "A" * 99
    coding = "ATG" + "GCT" * 99 + "TAA"
    unit = coding + "A" * 30
    sequence = unit * genes
    (root / "genome.fna").write_text(">contig\n" + sequence + "\n")
    (root / "source.faa").write_text(
        "".join(f">gene_{i:06}\n{protein}\n" for i in range(genes))
    )
    gff = ["##gff-version 3"]
    record = SeqRecord(
        Seq(sequence),
        id="contig",
        name="contig",
        description="synthetic engineering benchmark",
    )
    record.annotations["molecule_type"] = "DNA"
    record.features = []
    for i in range(genes):
        start = i * len(unit)
        gene = f"gene_{i:06}"
        record.features.append(
            SeqFeature(
                SimpleLocation(start, start + len(coding), strand=1),
                type="CDS",
                qualifiers={
                    "locus_tag": [gene],
                    "translation": [protein],
                    "transl_table": ["11"],
                    "codon_start": ["1"],
                },
            )
        )
        gff.append(
            f"contig\tProkka\tCDS\t{start + 1}\t{start + len(coding)}\t.\t+\t0\tID={gene};locus_tag={gene}"
        )
    (root / "source.gff").write_text(
        "\n".join(gff) + "\n##FASTA\n>contig\n" + sequence + "\n"
    )
    SeqIO.write(record, root / "source.gbk", "genbank")
    accessions = [f"SYN{i:06}" for i in range(samples)]
    bundles = []
    for accession in accessions:
        bundle = root / accession
        build_bundle(
            accession,
            11,
            root / "genome.fna",
            root / "source.faa",
            root / "source.gff",
            root / "source.gbk",
            bundle,
        )
        bundles.append(str(bundle))
    write_json(root / "bundles.json", bundles)
    write_tsv(
        root / "master.tsv",
        ("Accession", "Gcode"),
        [{"Accession": a, "Gcode": "11"} for a in accessions],
    )
    write_tsv(
        root / "status.tsv", ("accession",), [{"accession": a} for a in accessions]
    )
    plan = dict(
        schema_version=1,
        accessions=accessions,
        enabled_tools=[],
        tasks=[
            dict(
                accession=a,
                tool=t,
                action="skip",
                status="skipped_disabled",
                reason="tool_disabled",
                input_proteins=genes,
            )
            for a in accessions
            for t in TOOLS
        ],
    )
    plan["plan_id"] = identity(plan)
    write_json(root / "plan.json", plan)


def worker(args):
    import resource

    # Import the selected implementation before any candidate reporting modules.
    sys.path.insert(0, str(args.source / "bin"))
    import aggregate_annotations

    root = args.fixture
    bundles = [Path(p) for p in read_json(root / "bundles.json")]
    start = time.monotonic()
    if args.mode == "ani":
        from report_genomes import report

        report(
            root / "plan.json",
            root / "master.tsv",
            root / "status.tsv",
            [p / "bundle.json" for p in bundles],
            args.output,
        )
    else:
        aggregate_annotations.aggregate(
            root / "plan.json",
            root / "master.tsv",
            root / "status.tsv",
            bundles,
            [],
            args.output,
        )
    usage = resource.getrusage(resource.RUSAGE_SELF)
    print(
        json.dumps(
            dict(
                seconds=time.monotonic() - start,
                maxrss=usage.ru_maxrss,
                user_seconds=usage.ru_utime,
                system_seconds=usage.ru_stime,
                inblock=usage.ru_inblock,
                oublock=usage.ru_oublock,
            )
        )
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=64)
    parser.add_argument("--genes", type=int, default=3000)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--mode", choices=("full", "ani"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worker", action="store_true")
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    generate(args.fixture, args.samples, args.genes)
    runs = []
    for name, source, mode in [
        ("v047", args.baseline, "full"),
        ("v048_full", ROOT, "full"),
        ("v048_ani", ROOT, "ani"),
    ]:
        output = args.fixture / name
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--worker",
            "--fixture",
            str(args.fixture),
            "--source",
            str(source),
            "--mode",
            mode,
            "--output",
            str(output),
        ]
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        peak_spool = 0
        while process.poll() is None:
            peak_spool = max(
                peak_spool,
                sum(
                    p.stat().st_size
                    for d in args.fixture.glob(".annotation-aggregation-*")
                    for p in d.glob("*.sqlite")
                    if p.exists()
                ),
            )
            time.sleep(0.1)
        stdout, stderr = process.communicate()
        if process.returncode:
            raise RuntimeError(stderr)
        measured = json.loads(stdout.splitlines()[-1])
        measured.update(
            name=name,
            peak_sqlite_bytes=peak_spool,
            published_bytes=sum(
                p.stat().st_size for p in output.rglob("*") if p.is_file()
            ),
        )
        runs.append(measured)
    baseline = read_json(args.fixture / "v047/annotation_results.json")
    for name in ("v048_full", "v048_ani"):
        current = read_json(args.fixture / name / "annotation_results.json")
        if name == "v048_full" and set(current["tables"]) != set(baseline["tables"]):
            raise AssertionError("Full output table set changed")
        for table in current["tables"]:
            if digest(args.fixture / name / table) != digest(
                args.fixture / "v047" / table
            ):
                raise AssertionError(f"Output changed: {name}/{table}")
    result = dict(
        samples=args.samples,
        genes_per_sample=args.genes,
        proteins=args.samples * args.genes,
        synthetic=True,
        platform=sys.platform,
        python=sys.version.split()[0],
        runs=runs,
        shared_tables_identical=True,
        limits="Local storage, warm/uncontrolled OS caches; no native searches or HPC production run.",
    )
    write_json(args.fixture / "benchmark.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
