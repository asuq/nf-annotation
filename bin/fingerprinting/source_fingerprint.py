"""Bounded content readers with deterministic v0.4.5 publication identities."""

from __future__ import annotations

import hashlib
import os
import stat
import sys
import time
from collections import deque
from concurrent.futures import ThreadPoolExecutor, TimeoutError
from pathlib import Path

from annotation_common import AnnotationError, identity, read_json, write_json
from publication_references import REFERENCE_CONTRACT, REFERENCE_FILE, Publication

MAX_WORKERS = 8


def signature(info: os.stat_result) -> tuple:
    return (
        info.st_dev,
        info.st_ino,
        info.st_mode,
        info.st_size,
        info.st_mtime_ns,
        info.st_ctime_ns,
    )


def checked_digest(path: Path, before: os.stat_result) -> str:
    """Never follow a replaced leaf or accept a file modified during its read."""
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(descriptor, "rb") as handle:
        opened = os.fstat(handle.fileno())
        if not stat.S_ISREG(opened.st_mode) or signature(before) != signature(opened):
            raise AnnotationError(f"Publication changed during fingerprinting: {path}")
        checksum = hashlib.file_digest(handle, "sha256").hexdigest()
        if signature(before) != signature(os.fstat(handle.fileno())):
            raise AnnotationError(f"Publication changed during fingerprinting: {path}")
    if signature(before) != signature(path.lstat()):
        raise AnnotationError(f"Publication changed during fingerprinting: {path}")
    return checksum


class FingerprintScanner:
    """One pool per invocation; at most two queued reads per worker, no nesting."""

    def __init__(self, workers=1, *, cache=None, progress=None, requested_roots=1):
        if type(workers) is not int or not 1 <= workers <= MAX_WORKERS:
            raise AnnotationError(
                f"Fingerprint workers must be an integer from 1 to {MAX_WORKERS}"
            )
        self.workers = workers
        self.cache = {} if cache is None else cache
        self.progress = progress
        self.pool = None
        self.start = time.monotonic()
        self.last_report = 0.0
        self.metrics = {
            "operation": "fingerprint",
            "workers": workers,
            "phase": "starting",
            "requested_roots": requested_roots,
            "completed_roots": 0,
            "unique_roots_completed": 0,
            "shared_root_hits": 0,
            "files_hashed": 0,
            "bytes_hashed": 0,
            "scan_seconds": 0.0,
            "hash_seconds": 0.0,
            "peak_pending_reads": 0,
        }

    def __enter__(self):
        self.pool = ThreadPoolExecutor(max_workers=self.workers)
        self.report(force=True)
        return self

    def __exit__(self, kind, value, traceback):
        self.pool.shutdown(wait=True, cancel_futures=True)
        self.metrics["phase"] = "failed" if kind else "complete"
        self.report(force=True)

    def report(self, *, force=False):
        now = time.monotonic()
        self.metrics["seconds"] = now - self.start
        if not force and now - self.last_report < 10:
            return
        self.last_report = now
        if self.progress is not None:
            temporary = self.progress.with_name(self.progress.name + ".tmp")
            write_json(temporary, self.metrics)
            temporary.replace(self.progress)
            print(
                "Fingerprint: "
                + " ".join(
                    f"{key}={self.metrics[key]}"
                    for key in (
                        "phase",
                        "workers",
                        "completed_roots",
                        "requested_roots",
                        "files_hashed",
                        "bytes_hashed",
                        "seconds",
                    )
                ),
                file=sys.stderr,
                flush=True,
            )

    def hashes(self, paths):
        """Submit a bounded window and consume in sorted order, regardless of I/O order."""
        pending = deque()
        items = iter(sorted(paths, key=lambda item: item[0]))
        exhausted = False
        while pending or not exhausted:
            while not exhausted and len(pending) < 2 * self.workers:
                try:
                    path, before = next(items)
                except StopIteration:
                    exhausted = True
                    break
                pending.append(
                    (path, before, self.pool.submit(checked_digest, path, before))
                )
                self.metrics["peak_pending_reads"] = max(
                    self.metrics["peak_pending_reads"], len(pending)
                )
            if not pending:
                break
            path, before, future = pending.popleft()
            while True:
                try:
                    checksum = future.result(timeout=1)
                    break
                except TimeoutError:
                    if future.done():
                        raise
                    self.report()
            self.metrics["files_hashed"] += 1
            self.metrics["bytes_hashed"] += before.st_size
            self.metrics["current_files_hashed"] += 1
            self.metrics["current_bytes_hashed"] += before.st_size
            self.report()
            yield path, before, checksum

    def fingerprint(self, root: Path, *, active=()) -> dict:
        root = root.resolve(strict=True)
        if root in active or len(active) >= 64:
            raise AnnotationError(
                f"Circular or excessively deep publication reference: {root}"
            )
        if root in self.cache:
            self.metrics["shared_root_hits"] += 1
            return self.cache[root]
        manifest = root / "annotation_results.json"
        if manifest.is_symlink() or not manifest.is_file():
            raise AnnotationError(f"Invalid publication manifest: {manifest}")
        paths = [(manifest, manifest.stat())]
        reference = read_json(manifest).get("pipeline_contract") == REFERENCE_CONTRACT
        publication = Publication(root) if reference else None
        dependencies = []
        if publication:
            ref_path = root / REFERENCE_FILE
            paths.append((ref_path, ref_path.stat()))
            for source in publication.document["sources"]:
                current = self.fingerprint(Path(source["root"]), active=(*active, root))
                if current["sha256"] != source["sha256"]:
                    raise AnnotationError(
                        f"Referenced source changed: {source['root']}"
                    )
                dependencies.append(
                    {"root": current["root"], "sha256": current["sha256"]}
                )
        elif (root / REFERENCE_FILE).exists():
            raise AnnotationError(
                "Reference manifest lacks the reference publication contract"
            )
        self.metrics.update(
            current_root=str(root),
            phase="scanning",
            current_files_hashed=0,
            current_bytes_hashed=0,
            current_files_total=0,
            current_bytes_total=0,
        )
        self.report(force=True)
        start = time.monotonic()
        links, directories, pending = [], [], []
        for name in ("tables", "samples", "annotation_batches"):
            directory = root / name
            if name == "annotation_batches" and not directory.exists():
                continue
            if directory.is_symlink() or not directory.is_dir():
                raise AnnotationError(f"Invalid publication directory: {directory}")
            pending.append(directory)
        while pending:
            directory = pending.pop()
            directories.append((directory, directory.lstat()))
            with os.scandir(directory) as entries:
                for entry in entries:
                    if entry.is_symlink():
                        relative = Path(entry.path).relative_to(root).as_posix()
                        if publication is None or relative not in publication.links:
                            raise AnnotationError(
                                f"Publication fingerprint cannot follow a symlink: {entry.path}"
                            )
                        links.append([relative, "reference", os.readlink(entry.path)])
                    elif entry.is_dir(follow_symlinks=False):
                        pending.append(Path(entry.path))
                    else:
                        info = entry.stat(follow_symlinks=False)
                        if not stat.S_ISREG(info.st_mode):
                            raise AnnotationError(
                                f"Special publication file: {entry.path}"
                            )
                        paths.append((Path(entry.path), info))
                        self.metrics["current_files_total"] += 1
                        self.metrics["current_bytes_total"] += info.st_size
                    self.report()
        self.metrics["scan_seconds"] += time.monotonic() - start
        self.metrics.update(
            phase="hashing",
            current_files_total=len(paths),
            current_bytes_total=sum(info.st_size for _, info in paths),
        )
        self.report(force=True)
        start = time.monotonic()
        tree = hashlib.sha256()
        files = total = 0
        for path, before, checksum in self.hashes(paths):
            tree.update(
                identity(
                    [path.relative_to(root).as_posix(), before.st_size, checksum]
                ).encode()
            )
            files += 1
            total += before.st_size
        self.metrics["hash_seconds"] += time.monotonic() - start
        for directory, before in directories:
            if signature(before) != signature(directory.lstat()):
                raise AnnotationError(
                    f"Publication changed during fingerprinting: {directory}"
                )
        if publication:
            publication.check_links()
        for link in sorted(links):
            tree.update(identity(link).encode())
        result = {
            "root": str(root),
            "sha256": tree.hexdigest(),
            "files": files,
            "bytes": total,
        }
        if dependencies:
            result["references"] = dependencies
        self.cache[root] = result
        self.metrics["unique_roots_completed"] += 1
        return result


def fingerprint(
    root: Path,
    *,
    cache: dict | None = None,
    active: tuple[Path, ...] = (),
    workers: int = 1,
) -> dict:
    with FingerprintScanner(workers, cache=cache) as scanner:
        return scanner.fingerprint(root, active=active)


def fingerprint_sources(
    roots: list[Path], *, workers: int = 1, progress: Path | None = None
) -> tuple[list[dict], dict]:
    records = []
    with FingerprintScanner(
        workers, progress=progress, requested_roots=len(roots)
    ) as scanner:
        for root in roots:
            records.append(scanner.fingerprint(root))
            scanner.metrics["completed_roots"] += 1
            scanner.report(force=True)
    return records, scanner.metrics
