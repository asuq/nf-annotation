#!/usr/bin/env python3
"""Publish a bounded ANI reuse batch once, retaining independent ordinary files."""

from __future__ import annotations

import argparse
import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile
from pathlib import Path

from annotation_common import (
    AnnotationError,
    digest,
    read_json,
    validate_accession,
    write_json,
)
from independent_copy import copy_independent
from publication_references import Publication, qc_boundary, require_qc_provider


def regular_files(root: Path, *, exclude_annotation: bool = False) -> dict[str, Path]:
    """Reject links/special files while retaining cached directory-entry metadata."""
    if root.is_symlink() or not root.is_dir():
        raise AnnotationError(f"Expected an ordinary directory: {root}")
    result = {}
    pending = [root]
    while pending:
        directory = pending.pop()
        with os.scandir(directory) as entries:
            for entry in entries:
                if (
                    directory == root
                    and exclude_annotation
                    and entry.name == "annotation"
                ):
                    continue
                if entry.is_symlink():
                    raise AnnotationError(f"Invalid published entry: {entry.path}")
                if entry.is_dir(follow_symlinks=False):
                    pending.append(Path(entry.path))
                elif stat.S_ISREG(entry.stat(follow_symlinks=False).st_mode):
                    path = Path(entry.path)
                    result[path.relative_to(root).as_posix()] = path
                else:
                    raise AnnotationError(f"Invalid published entry: {entry.path}")
    return dict(sorted(result.items()))


def reusable_bundle(
    sample: Path, meta: dict, gcode: str, bundle: Path | None = None
) -> Path | None:
    """Bind a previously validated successful bundle to the current native inputs."""
    bundle = sample / "annotation/bundle" if bundle is None else bundle
    if not (bundle / "bundle.json").is_file():
        return None
    record = read_json(bundle / "bundle.json")
    code = Path(__file__).parent
    if (
        record.get("status") != "success"
        or record.get("accession") != meta["accession"]
        or str(record.get("genetic_code")) != gcode
        or record.get("producer_sha256")
        != digest(code / "prepare_annotation_bundle.py")
        or record.get("common_contract_sha256") != digest(code / "annotation_common.py")
        or record.get("prokka_log_sha256") != digest(sample / "prokka/prokka.log")
    ):
        return None
    native = {
        "genome.fasta": sample / f"staged/{meta['internal_id']}.fasta",
        "source.faa": sample / "prokka/prokka.faa",
        "source.gff": sample / "prokka/prokka.gff",
        "source.gbk": sample / "prokka/prokka.gbk",
    }
    if any(record["files"].get(name) != digest(path) for name, path in native.items()):
        return None
    return bundle


def build_bundle(sample: Path, meta: dict, gcode: str, output: Path) -> None:
    """Use the unchanged adapter when a matching reusable bundle is unavailable."""
    command = [
        sys.executable,
        str(Path(__file__).with_name("prepare_annotation_bundle.py")),
        "--accession",
        meta["accession"],
        "--gcode",
        gcode,
        "--genome",
        str(sample / f"staged/{meta['internal_id']}.fasta"),
        "--faa",
        str(sample / "prokka/prokka.faa"),
        "--gff",
        str(sample / "prokka/prokka.gff"),
        "--gbk",
        str(sample / "prokka/prokka.gbk"),
        "--prokka-log",
        str(sample / "prokka/prokka.log"),
        "--outdir",
        str(output),
    ]
    subprocess.run(command, check=True)


def import_sample(
    source: Path,
    meta: dict,
    gcode: str,
    destination: Path,
    channels: Path,
    *,
    view: Publication | None = None,
    mode: str = "self_contained",
    source_fingerprint: dict | None = None,
    qc_root: Path | None = None,
) -> dict:
    """Atomically install one sample, or verify an identical previous installation."""
    accession = meta["accession"]
    validate_accession(accession)
    if not re.fullmatch(r"[A-Za-z0-9_]+", meta["internal_id"]) or gcode not in (
        "4",
        "11",
        "NA",
    ):
        raise AnnotationError("Invalid internal ID or genetic code in import batch")
    if mode not in ("self_contained", "reference") or (
        mode == "reference" and (view is None or source_fingerprint is None)
    ):
        raise AnnotationError("Reference import requires a validated source snapshot")
    if mode == "reference":
        if qc_root is None:
            raise AnnotationError(
                "Reference publication requires --qc_publication_root"
            )
        qc_root = qc_boundary(qc_root, destination)
        if view.reference and view.version != 2:
            raise AnnotationError(
                "New reference outputs require original QC or v2 inputs; legacy absolute views are not migrated"
            )
        require_qc_provider(view.native_owner(accession), qc_root)
    source = (
        view.native_root(accession) if view is not None else source.resolve(strict=True)
    )
    if (
        destination == source
        or source in destination.parents
        or destination in source.parents
    ):
        raise AnnotationError("Publication and source must be separate directories")
    # Complete preflight validation precedes this task. Keep publication checks
    # here as well: inputs must remain unchanged throughout a workflow run.
    files = regular_files(source, exclude_annotation=True)
    parents = destination / "samples"
    if parents.is_symlink() or (parents.exists() and not parents.is_dir()):
        raise AnnotationError(
            "Publication samples directory must be an ordinary directory"
        )
    parents.mkdir(parents=True, exist_ok=True)
    target = parents / accession
    metrics = {
        "accession": accession,
        "files_copied": 0,
        "files_cloned": 0,
        "logical_bytes_copied": 0,
        "logical_bytes_cloned": 0,
        "bundle": "not_requested",
        "publication": "created",
        "links_created": 0,
        "logical_bytes_referenced": 0,
        "native_root": str(destination / "samples" / accession)
        if mode == "self_contained"
        else str(source),
        "bundle_root": str(destination / "samples" / accession / "annotation/bundle")
        if gcode in ("4", "11")
        else None,
    }
    with tempfile.TemporaryDirectory(prefix=".ani-batch-", dir=parents) as temporary:
        temporary = Path(temporary)
        bundle = None
        if gcode in ("4", "11"):
            candidate = (
                view.path(f"samples/{accession}/annotation/bundle")
                if view is not None
                else None
            )
            bundle = reusable_bundle(source, meta, gcode, candidate)
            if bundle is None:
                bundle = temporary / "rebuilt"
                build_bundle(source, meta, gcode, bundle)
                metrics["bundle"] = "rebuilt"
            else:
                metrics["bundle"] = "reused"
            files.update(
                {
                    "annotation/bundle/" + name: path
                    for name, path in regular_files(bundle).items()
                }
            )
        if mode == "reference":
            native_owner = view.native_owner(accession)
            entries = sorted(
                path.name for path in source.iterdir() if path.name != "annotation"
            )
            links = {name: source / name for name in entries}
            reference = {
                "native": view.snapshot(native_owner, source_fingerprint),
                "entries": entries,
                "bundle": None,
            }
            if bundle is not None:
                if (
                    metrics["bundle"] == "reused"
                    and not view.bundle_owner(accession).reference
                ):
                    owner = view.bundle_owner(accession)
                    require_qc_provider(owner, qc_root)
                    links["annotation/bundle"] = owner.path(
                        f"samples/{accession}/annotation/bundle"
                    )
                    reference["bundle"] = view.snapshot(owner, source_fingerprint)
                else:
                    reference["bundle"] = {"local": True}
            metrics["reference"] = reference
            metrics["logical_bytes_referenced"] = sum(
                path.stat().st_size
                for name, path in files.items()
                if not name.startswith("annotation/") or "annotation/bundle" in links
            )
            if target.exists() or target.is_symlink():
                if target.is_symlink() or not target.is_dir():
                    raise AnnotationError(
                        f"Invalid reference sample directory: {accession}"
                    )
                expected_entries = set(entries) | (
                    {"annotation"} if bundle is not None else set()
                )
                if {path.name for path in target.iterdir()} != expected_entries:
                    raise AnnotationError(
                        f"Reference publication entries changed for {accession}"
                    )
                for name, provider in links.items():
                    path = target / name
                    if not path.is_symlink() or os.readlink(path) != os.path.relpath(
                        provider, path.parent
                    ):
                        raise AnnotationError(
                            f"Reference publication changed for {accession}: {name}"
                        )
                if bundle is not None and reference["bundle"] == {"local": True}:
                    existing = regular_files(target / "annotation/bundle")
                    expected = regular_files(bundle)
                    if set(existing) != set(expected) or any(
                        digest(existing[name]) != digest(path)
                        for name, path in expected.items()
                    ):
                        raise AnnotationError(
                            f"Local reference bundle changed for {accession}"
                        )
                metrics["publication"] = "verified_existing"
            else:
                prepared = temporary / "sample"
                prepared.mkdir()
                for name, provider in links.items():
                    link = prepared / name
                    link.parent.mkdir(parents=True, exist_ok=True)
                    link.symlink_to(
                        os.path.relpath(provider, (target / name).parent),
                        target_is_directory=provider.is_dir(),
                    )
                    metrics["links_created"] += 1
                if bundle is not None and reference["bundle"] == {"local": True}:
                    (prepared / "annotation").mkdir(exist_ok=True)
                    metrics["local_bundle_bytes"] = sum(
                        path.stat().st_size for path in regular_files(bundle).values()
                    )
                    if metrics["bundle"] == "rebuilt":
                        bundle.rename(prepared / "annotation/bundle")
                    else:
                        shutil.copytree(bundle, prepared / "annotation/bundle")
                        metrics["files_copied"] += len(regular_files(bundle))
                        metrics["logical_bytes_copied"] += metrics["local_bundle_bytes"]
                prepared.rename(target)
        else:
            if target.exists() or target.is_symlink():
                existing = regular_files(target)
                if set(existing) != set(files) or any(
                    digest(existing[name]) != digest(path)
                    for name, path in files.items()
                ):
                    raise AnnotationError(
                        f"Existing publication differs for {accession}; use a fresh --outdir"
                    )
                metrics["publication"] = "verified_existing"
            else:
                prepared = temporary / "sample"
                prepared.mkdir()
                for name, path in files.items():
                    method = copy_independent(path, prepared / name)
                    metrics["files_" + method] += 1
                    metrics["logical_bytes_" + method] += path.stat().st_size
                # An ordinary directory rename on the destination filesystem makes
                # each sample visible only after its complete payload is present.
                prepared.rename(target)

    channels.mkdir(parents=True, exist_ok=True)
    for relative, suffix in (
        ("16s/best_16S.fna", "best_16S.fna"),
        ("16s/16S_status.tsv", "16S_status.tsv"),
        ("checkm2/checkm2_summary.tsv", "checkm2_summary.tsv"),
    ):
        shutil.copyfile(target / relative, channels / f"{meta['internal_id']}_{suffix}")
    return metrics


def run(
    records: list[dict],
    sources: list[Path],
    destination: Path,
    channels: Path,
    audit: Path,
    *,
    mode: str = "self_contained",
    fingerprints: list[dict] | None = None,
    qc_root: Path | None = None,
) -> None:
    """Report individual failures and fail the whole batch before ANI can start."""
    if not records or not sources:
        raise AnnotationError("Batch records and source roots must be nonempty")
    accessions = [record["meta"]["accession"] for record in records]
    if len(set(accessions)) != len(accessions):
        raise AnnotationError("Duplicate accessions in import batch")
    snapshots = {record["root"]: record for record in (fingerprints or [])}
    contexts = {}
    results = []
    for record in records:
        try:
            validate_accession(record["meta"]["accession"])
            index = record["source_index"]
            if not isinstance(index, int) or not 0 <= index < len(sources):
                raise AnnotationError("Invalid source index in import batch")
            root = sources[index].resolve(strict=True)
            if root not in contexts:
                contexts[root] = (
                    Publication(root)
                    if (root / "annotation_results.json").is_file()
                    else None
                )
            source = root / "samples" / record["meta"]["accession"]
            result = import_sample(
                source,
                record["meta"],
                record["gcode"],
                destination,
                channels,
                view=contexts[root],
                mode=mode,
                source_fingerprint=snapshots.get(str(root)),
                qc_root=qc_root,
            )
            result.update(status="done", meta=record["meta"], gcode=record["gcode"])
        except (
            AnnotationError,
            OSError,
            ValueError,
            subprocess.CalledProcessError,
        ) as error:
            result = {
                "accession": record["meta"]["accession"],
                "status": "failed",
                "reason": str(error),
            }
            print(f"{result['accession']}: {error}", file=sys.stderr)
        results.append(result)
    write_json(
        audit,
        {
            "samples": results,
            "complete": all(row["status"] == "done" for row in results),
        },
    )
    if any(row["status"] != "done" for row in results):
        raise AnnotationError(
            "ANI import batch has failed samples; see its per-sample audit"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("records", "sources", "destination", "channels", "audit"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument(
        "--mode", choices=("self_contained", "reference"), default="self_contained"
    )
    parser.add_argument("--source-validation", type=Path)
    parser.add_argument("--qc-publication-root", type=Path)
    args = parser.parse_args()
    try:
        run(
            read_json(args.records),
            [Path(path) for path in read_json(args.sources)],
            args.destination.resolve(strict=True),
            args.channels,
            args.audit,
            mode=args.mode,
            qc_root=args.qc_publication_root,
            fingerprints=read_json(args.source_validation)["fingerprints"]
            if args.source_validation
            else None,
        )
    except (AnnotationError, OSError, ValueError) as error:
        print(error, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
