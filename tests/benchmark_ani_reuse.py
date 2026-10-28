#!/usr/bin/env python3
"""Reproducible local ANI-bookkeeping benchmark; no biological inputs or HPC jobs."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import summarise_checkm2
import test_cohort_update_integration as workflow_fixture
from annotation_common import digest, read_json, write_json
from annotation_source import validate_source
from annotation_workflow_fixture import (
    install_synthetic_container_engine,
    publish_disabled,
)
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord
from prepare_annotation_bundle import build_bundle
from test_prepare_cohort_update import LINEAGE, PrepareCohortUpdateTestCase, table


def rows(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def make_source(root: Path, count: int) -> tuple[Path, list[str], Path]:
    """Create complete publications with one canonical gene per synthetic genome."""
    case = PrepareCohortUpdateTestCase(methodName="runTest")
    case.setUp()
    case.root = root
    case.source = root / "source"
    accessions = [f"SYN{index:06d}" for index in range(count)]
    try:
        case.published_cohort(accessions)
        master_path = case.source / "tables/master_table.tsv"
        status_path = case.source / "tables/sample_status.tsv"
        masters = rows(master_path)
        statuses = rows(status_path)
        sequence = "ATGTGAGCTTAA" + "ACGT" * 397
        genome = ">contig1\n" + sequence + "\n"
        record = SeqRecord(
            Seq(sequence),
            id="contig1",
            name="contig1",
            description="synthetic performance fixture",
        )
        record.annotations["molecule_type"] = "DNA"
        record.features = [
            SeqFeature(
                SimpleLocation(0, 12, strand=1),
                type="CDS",
                qualifiers={
                    "locus_tag": ["gene_1"],
                    "translation": ["MWA"],
                    "transl_table": ["4"],
                    "codon_start": ["1"],
                },
            )
        ]
        bundles = []
        for index, accession in enumerate(accessions):
            sample = case.source / "samples" / accession
            original = case.genomes[accession]
            staged = sample / f"staged/{accession}.fasta"
            original.write_text(genome)
            staged.write_text(genome)
            (sample / "16s/best_16S.fna").write_text(
                ">16S_rRNA::contig1:0-1500(+)\n" + sequence[:1500] + "\n"
            )
            faa, gff, gbk = [
                sample / f"prokka/prokka.{suffix}" for suffix in ("faa", "gff", "gbk")
            ]
            faa.write_text(">gene_1\nMWA\n")
            gff.write_text(
                "##gff-version 3\ncontig1\tProkka\tCDS\t1\t12\t.\t+\t0\tID=gene_1;locus_tag=gene_1\n##FASTA\n"
                + genome
            )
            SeqIO.write(record, gbk, "genbank")
            reports = {}
            for code in (4, 11):
                path = sample / f"checkm2_gcode{code}/quality_report.tsv"
                report = rows(path)
                report[0]["Completeness"] = "20"
                table(path, list(report[0]), report)
                reports[code] = summarise_checkm2.parse_report(path)
            qc = summarise_checkm2.build_output_row(
                accession,
                reports[4],
                reports[11],
                warnings=[],
                gcode_rule="mean_gene_length_ratio",
            )
            summarise_checkm2.write_output(sample / "checkm2/checkm2_summary.tsv", qc)
            masters[index].update(
                {key: value for key, value in qc.items() if key in masters[index]}
            )
            statuses[index]["low_quality"] = "true"
            bundle = sample / "annotation/bundle"
            result = build_bundle(accession, 4, staged, faa, gff, gbk, bundle)
            result["prokka_log_sha256"] = digest(sample / "prokka/prokka.log")
            write_json(bundle / "bundle.json", result)
            bundles.append(bundle)
        table(master_path, list(masters[0]), masters)
        table(status_path, list(statuses[0]), statuses)
        publish_disabled(case.source, accessions, bundles)
        validate_source(case.source)
        metadata = root / "metadata.tsv"
        table(
            metadata,
            [
                "accession",
                "Tax_ID",
                "Organism_Name",
                "Assembly_Level",
                "Atypical_Warnings",
            ],
            [
                {
                    "accession": accession,
                    "Tax_ID": "8",
                    "Organism_Name": "Synthetic bacterium",
                    "Assembly_Level": "Scaffold",
                    "Atypical_Warnings": "NA",
                }
                for accession in accessions
            ],
        )
        return case.source, accessions, metadata
    finally:
        case.doCleanups()


def fixture_project(repo: Path) -> Path:
    """Use real orchestration and bookkeeping; synthetic external callers only."""
    previous = workflow_fixture.ROOT
    workflow_fixture.ROOT = repo
    try:
        workflow_fixture.CohortUpdateIntegrationTestCase.setUpClass()
        project = workflow_fixture.CohortUpdateIntegrationTestCase.project
    finally:
        workflow_fixture.ROOT = previous
    install_synthetic_container_engine(project)
    # All generated inputs are the same 1,600-nt sequence. Avoid measuring the
    # Python fake-seqtk startup cost thousands of times as if it were native C.
    module = project / "modules/local/calculate_assembly_stats.nf"
    original = (repo / "modules/local/calculate_assembly_stats.nf").read_text()
    stub = original.split("\n    stub:", 1)[1]
    sequence = "ATGTGAGCTTAA" + "ACGT" * 397
    gc = 100 * (sequence.count("G") + sequence.count("C")) / len(sequence)
    stub = stub.replace("80", "1600").replace("1600\\\\t50", f"1600\\\\t{gc}")
    module.write_text(module.read_text().rsplit("}", 1)[0] + "\n    stub:" + stub)
    return project


def tree_bytes(root: Path) -> dict:
    apparent = allocated = count = 0
    for path in root.rglob("*"):
        if path.is_symlink() or not path.is_file():
            continue
        info = path.stat()
        apparent += info.st_size
        allocated += info.st_blocks * 512
        count += 1
    return {
        "files": count,
        "apparent_bytes": apparent,
        "reported_allocated_bytes": allocated,
    }


def run_cohort(
    project: Path,
    root: Path,
    source: Path,
    accessions: list[str],
    metadata: Path,
    name: str,
    *,
    resume: bool,
    docker: bool,
) -> dict:
    root.mkdir(parents=True, exist_ok=True)
    samples = root / f"{name}.csv"
    with samples.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["accession", "is_new", "assembly_level", "genome_fasta"])
        for accession in accessions:
            writer.writerow(
                [accession, "false", "NA", source.parent / f"inputs/{accession}.fasta"]
            )
    output = root / name
    command = [
        "nextflow",
        "run",
        str(project),
        "-profile",
        "test,docker" if docker else "test",
        "-stub-run",
        "--domain",
        "bacteria",
        "--pipeline_stage",
        "ani",
        "--busco_lineages",
        LINEAGE,
        "--ani_16s_policy",
        "ignore",
        "--task_attempts",
        "1",
        "--sample_csv",
        str(samples),
        "--metadata",
        str(metadata),
        "--update_from",
        str(source),
        "--outdir",
        str(output),
        "-work-dir",
        str(root / "work"),
    ]
    if resume:
        command.append("-resume")
    start = time.monotonic()
    with (root / f"{name}.log").open("w") as log:
        subprocess.run(
            command,
            cwd=root,
            check=True,
            timeout=3600,
            stdout=log,
            stderr=subprocess.STDOUT,
            env=dict(
                os.environ,
                PATH=str(project / "bin") + os.pathsep + os.environ["PATH"],
                NXF_ANSI_LOG="false",
                NXF_DISABLE_CHECK_LATEST="true",
            ),
        )
    elapsed = time.monotonic() - start
    trace = rows(output / "pipeline_info/trace.tsv")
    names = Counter(
        row["name"].split(" (")[0].split(":")[-1]
        for row in trace
        if row["status"] != "CACHED"
    )
    audits = [
        read_json(path)
        for path in (output / "pipeline_info/ani_import").glob("batch_*.json")
    ]
    metrics = Counter()
    for audit in audits:
        for sample in audit["samples"]:
            for key in (
                "files_copied",
                "files_cloned",
                "logical_bytes_copied",
                "logical_bytes_cloned",
            ):
                metrics[key] += sample[key]
            metrics["bundles_" + sample["bundle"]] += 1
    validate_source(output)
    return {
        "name": name,
        "samples": len(accessions),
        "seconds": elapsed,
        "submitted_tasks": sum(names.values()),
        "tasks": dict(names),
        "cached_tasks": sum(row["status"] == "CACHED" for row in trace),
        "import_metrics": dict(metrics),
        "fingerprinted_bytes": sum(
            record["bytes"]
            for record in read_json(output / "tables/source_validation.json")[
                "fingerprints"
            ]
        )
        if (output / "tables/source_validation.json").exists()
        else None,
        "work": tree_bytes(root / "work"),
        "output": str(output),
    }


def compare(left: Path, right: Path) -> list[str]:
    """Require byte-identical primary reports, matrices, IDs and coordinates."""
    matched = []
    for path in sorted((left / "tables").rglob("*.tsv")):
        if path.name == "tool_and_db_versions.tsv":
            continue
        relative = path.relative_to(left)
        if path.read_bytes() != (right / relative).read_bytes():
            raise AssertionError(f"Scientific table mismatch: {relative}")
        matched.append(str(relative))
    for path in sorted((left / "cohort/16s").iterdir()):
        if path.suffix not in (".fna", ".tsv"):
            continue
        relative = path.relative_to(left)
        if path.read_bytes() != (right / relative).read_bytes():
            raise AssertionError(f"16S output mismatch: {relative}")
        matched.append(str(relative))
    return matched


def code_snapshot(repo: Path) -> dict:
    """Record measured production code without relying on an uncommitted tree label."""
    files = {}
    for directory in ("bin", "lib", "modules", "subworkflows", "conf"):
        for path in sorted((repo / directory).rglob("*")):
            if path.is_file() and path.suffix in (
                ".py",
                ".sh",
                ".nf",
                ".groovy",
                ".config",
            ):
                files[str(path.relative_to(repo))] = digest(path)
    for name in ("main.nf", "nextflow.config", "reannotate.nf"):
        files[name] = digest(repo / name)
    return files


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument(
        "--source-fixture",
        type=Path,
        help="Reuse a previously generated synthetic benchmark source",
    )
    parser.add_argument("--comparison-size", type=int, default=128)
    parser.add_argument("--sizes", type=int, nargs="+", default=[2597, 6015, 6636])
    parser.add_argument("--docker-control", action="store_true")
    args = parser.parse_args()
    if (
        args.output.exists()
        or min([args.comparison_size, *args.sizes]) < 1
        or max(args.sizes) > 10000
    ):
        parser.error("Use a new output directory and sizes between 1 and 10000")
    if shutil.disk_usage(args.output.parent).free < 8 * 1024**3:
        parser.error(
            "The benchmark requires at least 8 GiB free for disposable fixtures"
        )
    args.output.mkdir()
    if args.source_fixture:
        source = args.source_fixture.resolve(strict=True)
        previous = read_json(source.parent / "benchmark.json")
        if previous.get("synthetic") is not True or previous.get(
            "source_manifest_sha256"
        ) != digest(source / "annotation_results.json"):
            parser.error(
                "--source-fixture requires an unchanged synthetic benchmark publication"
            )
        accessions = [
            row["accession"] for row in rows(source / "tables/validated_samples.tsv")
        ]
        metadata = source.parent / "metadata.tsv"
        if max(args.comparison_size, *args.sizes) > len(accessions):
            parser.error("The synthetic source has too few samples")
    else:
        source, accessions, metadata = make_source(
            args.output, max(args.comparison_size, *args.sizes)
        )
    source_manifest = digest(source / "annotation_results.json")
    old = fixture_project(args.baseline)
    new = fixture_project(ROOT)
    report = {
        "synthetic": True,
        "native_assembly_stats": False,
        "native_fastani": False,
        "all_samples_quality_excluded": True,
        "source_manifest_sha256": source_manifest,
        "runs": [],
        "baseline_code": code_snapshot(args.baseline),
        "candidate_code": code_snapshot(ROOT),
        "python": sys.version,
        "fixture_genome_nt": 1600,
        "fixture_proteins_per_genome": 1,
    }
    old_run = run_cohort(
        old,
        args.output / "baseline",
        source,
        accessions[: args.comparison_size],
        metadata,
        "comparison",
        resume=False,
        docker=args.docker_control,
    )
    report["runs"].append(old_run)
    new_run = run_cohort(
        new,
        args.output / "candidate",
        source,
        accessions[: args.comparison_size],
        metadata,
        "comparison",
        resume=False,
        docker=args.docker_control,
    )
    report["runs"].append(new_run)
    report["identical_tables"] = compare(
        Path(old_run["output"]), Path(new_run["output"])
    )
    write_json(args.output / "benchmark.json", report)
    for size in args.sizes:
        run = run_cohort(
            new,
            args.output / "candidate",
            source,
            accessions[:size],
            metadata,
            f"cohort-{size}",
            resume=True,
            docker=args.docker_control,
        )
        report["runs"].append(run)
        assert run["tasks"].get("IMPORT_ANI_BATCH", 0) <= 8
        assert not set(run["tasks"]) & {
            "IMPORT_PUBLISHED_SAMPLE",
            "PREPARE_ANNOTATION_BUNDLE",
            "SUMMARISE_BUSCO",
            "PLAN_ANNOTATION_PART",
        }
        assert run["submitted_tasks"] <= 40
        assert [
            row["accession"]
            for row in rows(Path(run["output"]) / "tables/master_table.tsv")
        ] == accessions[:size]
        write_json(args.output / "benchmark.json", report)
    assert digest(source / "annotation_results.json") == source_manifest
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
