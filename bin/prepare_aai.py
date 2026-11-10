#!/usr/bin/env python3
"""Prepare bounded batches of native marker profiles from existing protein bundles."""

from __future__ import annotations

import argparse
import concurrent.futures
import inspect
import multiprocessing
import shutil
import tempfile
import time
from pathlib import Path

import sample_policy
from annotation_common import AnnotationError, digest, identity, read_json, write_json
from annotation_path_lists import read_path_list
from atypical_warnings import classify_atypical_warnings
from build_fastani_inputs import (
    choose_checkm2_fields,
    detect_metadata_value,
    load_metadata_index,
)
from build_master_table import choose_assembly_level
from independent_copy import copy_independent

RUNTIME = Path("/opt/nf-annotation/provenance/fastaai.json")


def native_init():
    import fastaai.fastaai as native  # ty: ignore[unresolved-import]  # Provided by the pinned native image.

    native.hmm_preproc_initializer(native.find_hmm(), False)


def extract(job):
    import fastaai.fastaai as native  # ty: ignore[unresolved-import]  # Provided by the pinned native image.

    protein, temporary = job
    native.prepare_directories(str(temporary), "protein", "build", make_crystals=True)
    record = native.input_file(str(protein), str(temporary), make_crystal=True)
    record.set_protein(str(protein))
    record.preprocess()
    if record.initial_state != "protein":
        raise AnnotationError("FastAAI must start from existing proteins")
    count = len(record.best_hits_kmers)
    crystal = temporary / "crystals" / (record.basename + "_faai_crystal.txt")
    if count and not crystal.is_file():
        raise AnnotationError("FastAAI omitted a marker profile")
    if not count:
        return 0, None
    selected = {}
    with crystal.open() as handle:
        for line in handle:
            protein_id, accession, sequence = line.rstrip().split("\t")
            selected[protein_id, accession] = sequence
    with (temporary / "markers.faa").open("w") as handle:
        for (protein_id, _), sequence in selected.items():
            handle.write(">" + protein_id + "\n" + sequence + "\n")
    observed = set()
    with (
        (temporary / "hmms" / (record.basename + ".hmm")).open() as source,
        (temporary / "markers.hmm").open("w") as target,
    ):
        for line in source:
            fields = line.split()
            if (
                len(fields) >= 9
                and (fields[0], fields[3].replace(".", "_")) in selected
            ):
                target.write(line)
                observed.add((fields[0], fields[3].replace(".", "_")))
    if observed != set(selected):
        raise AnnotationError("Native HMM evidence omitted a selected marker")
    return count, temporary


def profile_key(accession, protein_sha, method):
    return identity(
        {"accession": accession, "protein_sha256": protein_sha, "method": method}
    )


def read_profile(root, key):
    path = root / key
    if path.is_symlink():
        raise AnnotationError("FastAAI profiles must be independent directories")
    if not path.exists():
        return None
    record = read_json(path / "profile.json")
    if (
        record.get("profile_id") != key
        or profile_key(record["accession"], record["protein_sha256"], record["method"])
        != key
    ):
        raise AnnotationError("FastAAI profile identity mismatch")
    count = record["marker_count"]
    if type(count) is not int or not 0 <= count <= 122:
        raise AnnotationError("Invalid FastAAI marker count")
    expected = {"markers.faa", "markers.hmm"} if count else set()
    if set(record["files"]) != expected or any(
        digest(path / name) != sha for name, sha in record["files"].items()
    ):
        raise AnnotationError("FastAAI marker profile changed")
    return record


def prepare(args):
    start = time.monotonic()
    runtime = read_json(RUNTIME)
    method = {
        key: runtime[key]
        for key in (
            "upstream_commit",
            "hmm_sha256",
            "pyhmmer",
            "numpy",
            "environment_sha256",
        )
    }
    method["profile_format"] = "selected_protein_hmm_v1"
    method["extractor_sha256"] = identity(
        [inspect.getsource(native_init), inspect.getsource(extract)]
    )
    selected = read_json(args.records)
    _, metadata = load_metadata_index(args.metadata)
    bundles = {}
    for path in read_path_list(args.bundles):
        record = read_json(path / "bundle.json")
        bundles[record["accession"]] = (path, record)
    publication = args.destination.resolve()
    for relative in ("cohort", "cohort/fastaai", "cohort/fastaai/profiles"):
        if (publication / relative).is_symlink():
            raise AnnotationError(
                "FastAAI publication directories must not be symlinks"
            )
    destination = publication / "cohort/fastaai/profiles"
    destination.mkdir(parents=True, exist_ok=True)
    source_roots = [
        p.resolve() / "cohort/fastaai/profiles" for p in read_path_list(args.sources)
    ]
    rows, pending = [], []
    reused = 0
    with tempfile.TemporaryDirectory(prefix="aai-preprocessing-") as temp:
        temp = Path(temp)
        for sample in selected:
            acc = sample["accession"]
            q = sample["qc"]
            gcode, completeness, contamination = choose_checkm2_fields(q)
            m = metadata.get(acc, {})
            bundle_path, bundle = bundles.get(acc, (None, {}))
            reasons = []
            if gcode not in ("4", "11"):
                reasons.append("gcode_na")
            if q.get("Low_quality") == "true":
                reasons.append("low_quality")
            elif q.get("Low_quality") != "false":
                reasons.append("missing_low_quality")
            atypical, exception = classify_atypical_warnings(
                detect_metadata_value(m, "Atypical_Warnings")
            )
            if atypical and not exception:
                reasons.append("atypical")
            if bundle.get("status") != "success":
                reasons.append("no_valid_protein_bundle")
            primary = sample_policy.primary_column(sample, args.primary_busco_column)
            quality = {
                "path": "proteins.faa",
                "gcode": gcode,
                "checkm2_completeness": completeness,
                "checkm2_contamination": contamination,
                "assembly_level": choose_assembly_level(sample, m),
                "n50": "NA",
                "scaffolds": "NA",
                "genome_size": "NA",
                "organism_name": detect_metadata_value(m, "Organism_Name"),
                "busco": "NA",
            }
            row = {
                "accession": acc,
                "domain": sample["annotation_domain"]
                if "annotation_domain" in sample
                else sample_policy.domain(args.domain),
                "quality": quality,
                "primary_busco_column": primary,
                "status": "excluded" if reasons else "ready",
                "reason": ";".join(reasons) or "NA",
                "profile": None,
                "marker_count": 0,
            }
            rows.append(row)
            if reasons:
                continue
            protein_sha = bundle["files"]["proteins.faa"]
            key = profile_key(acc, protein_sha, method)
            existing = read_profile(destination, key)
            if existing is None:
                for source in source_roots:
                    existing = read_profile(source, key)
                    if existing is not None:
                        with tempfile.TemporaryDirectory(
                            prefix="." + key + ".", dir=destination
                        ) as temporary:
                            staging = Path(temporary)
                            for name in (
                                ("profile.json", "markers.faa", "markers.hmm")
                                if existing["marker_count"]
                                else ("profile.json",)
                            ):
                                copy_independent(source / key / name, staging / name)
                            staging.rename(destination / key)
                        break
            row["profile"] = str(destination / key)
            if existing is not None:
                row["marker_count"] = existing["marker_count"]
                reused += 1
            else:
                work = temp / key
                work.mkdir()
                assert bundle_path is not None
                pending.append(
                    (row, key, protein_sha, bundle_path / "proteins.faa", work)
                )
        if pending:
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=args.threads,
                mp_context=multiprocessing.get_context("spawn"),
                initializer=native_init,
            ) as pool:
                results = pool.map(extract, ((p[3], p[4]) for p in pending))
                for item, result in zip(pending, results, strict=True):
                    row, key, protein_sha, _, _ = item
                    count, crystal = result
                    record = {
                        "profile_id": key,
                        "accession": row["accession"],
                        "protein_sha256": protein_sha,
                        "method": method,
                        "marker_count": count,
                        "files": {},
                    }
                    staging = Path(tempfile.mkdtemp(prefix="." + key, dir=destination))
                    if count:
                        for name in ("markers.faa", "markers.hmm"):
                            shutil.copyfile(crystal / name, staging / name)
                            record["files"][name] = digest(staging / name)
                    write_json(staging / "profile.json", record)
                    staging.rename(destination / key)
                    row["marker_count"] = count
        for row in rows:
            if row["status"] == "ready" and row["marker_count"] == 0:
                row.update(status="no_markers", reason="no_usable_markers")
    write_json(args.output, {"samples": rows, "method": method})
    write_json(
        args.output.parent / "preprocessing_metrics.json",
        {
            "samples": len(rows),
            "extracted": len(pending),
            "reused": reused,
            "seconds": time.monotonic() - start,
            "workers": args.threads,
        },
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in (
        "records",
        "metadata",
        "bundles",
        "sources",
        "destination",
        "output",
    ):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--primary-busco-column", required=True)
    parser.add_argument(
        "--domain", choices=("bacteria", "archaea", "metadata"), required=True
    )
    parser.add_argument("--threads", type=int, required=True)
    args = parser.parse_args()
    try:
        prepare(args)
    except (AnnotationError, OSError, ValueError, KeyError) as error:
        parser.exit(1, str(error) + "\n")


if __name__ == "__main__":
    main()
