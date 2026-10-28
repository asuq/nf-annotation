#!/usr/bin/env python3
"""Publish a bounded ANI reuse batch once, retaining independent ordinary files."""

from __future__ import annotations

import argparse
import ctypes
import errno
import functools
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


@functools.cache
def darwin_clonefile():
    """Resolve the optional platform primitive once per allocation."""
    clonefile = ctypes.CDLL(None, use_errno=True).clonefile
    clonefile.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int]
    clonefile.restype = ctypes.c_int
    return clonefile


def copy_independent(source: Path, destination: Path) -> str:
    """Use filesystem copy-on-write when supported; never share writable inodes."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    unsupported = {
        errno.EXDEV,
        errno.EOPNOTSUPP,
        errno.ENOTTY,
        errno.EINVAL,
        errno.ENOSYS,
    }
    try:
        if sys.platform == "linux":
            import fcntl

            with source.open("rb") as src, destination.open("xb") as dst:
                fcntl.ioctl(dst.fileno(), 0x40049409, src.fileno())  # FICLONE
        elif sys.platform == "darwin":
            if darwin_clonefile()(os.fsencode(source), os.fsencode(destination), 0):
                code = ctypes.get_errno()
                raise OSError(code, os.strerror(code))
        else:
            raise OSError(errno.ENOSYS, "Copy-on-write is unavailable")
    except OSError as error:
        if error.errno not in unsupported:
            raise
        destination.unlink(missing_ok=True)
        shutil.copy2(source, destination)
        return "copied"
    shutil.copystat(source, destination)
    return "cloned"


def reusable_bundle(sample: Path, meta: dict, gcode: str) -> Path | None:
    """Bind a previously validated successful bundle to the current native inputs."""
    bundle = sample / "annotation/bundle"
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
    source: Path, meta: dict, gcode: str, destination: Path, channels: Path
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
    source = source.resolve(strict=True)
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
    }
    with tempfile.TemporaryDirectory(prefix=".ani-batch-", dir=parents) as temporary:
        temporary = Path(temporary)
        if gcode in ("4", "11"):
            bundle = reusable_bundle(source, meta, gcode)
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
        if target.exists() or target.is_symlink():
            existing = regular_files(target)
            if set(existing) != set(files) or any(
                digest(existing[name]) != digest(path) for name, path in files.items()
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
) -> None:
    """Report individual failures and fail the whole batch before ANI can start."""
    if not records or not sources:
        raise AnnotationError("Batch records and source roots must be nonempty")
    accessions = [record["meta"]["accession"] for record in records]
    if len(set(accessions)) != len(accessions):
        raise AnnotationError("Duplicate accessions in import batch")
    results = []
    for record in records:
        try:
            validate_accession(record["meta"]["accession"])
            index = record["source_index"]
            if not isinstance(index, int) or not 0 <= index < len(sources):
                raise AnnotationError("Invalid source index in import batch")
            source = sources[index] / "samples" / record["meta"]["accession"]
            result = import_sample(
                source, record["meta"], record["gcode"], destination, channels
            )
            result["status"] = "done"
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
    args = parser.parse_args()
    try:
        run(
            read_json(args.records),
            [Path(path) for path in read_json(args.sources)],
            args.destination.resolve(strict=True),
            args.channels,
            args.audit,
        )
    except (AnnotationError, OSError, ValueError) as error:
        print(error, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
