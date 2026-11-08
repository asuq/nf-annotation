"""Stage a prepared annotation database once per native eggNOG task."""

import time
from pathlib import Path

from annotation_common import AnnotationError

ANNOTATION_FILES = (
    "eggnog.db",
    "eggnog.db.fieldpresence.bin",
    "eggnog.db.taxids.bin",
    "eggnog.taxa.db",
    "eggnog.taxa.db.traverse.pkl",
    "go-basic.obo",
)


def stage_database(resource: Path, prepared: dict, destination: Path) -> dict:
    """Copy annotation files once; hash while copying, never copy the DIAMOND index."""
    import hashlib

    records = {row["path"]: row for row in prepared["contract"]["files"]}
    if not set(ANNOTATION_FILES) <= records.keys():
        raise AnnotationError("Missing annotation database files in prepared resource")
    destination.mkdir()
    copied = {}
    for name in ANNOTATION_FILES:
        started = time.monotonic()
        checksum, size = hashlib.md5(), 0
        target = destination / name
        with (resource / name).open("rb") as source, target.open("xb") as output:
            while block := source.read(8 * 1024 * 1024):
                output.write(block)
                checksum.update(block)
                size += len(block)
        if (
            size != records[name]["bytes"]
            or checksum.hexdigest() != records[name]["md5"]
        ):
            raise AnnotationError(f"Changed annotation database while staging: {name}")
        copied[name] = {
            "bytes": size,
            "md5": checksum.hexdigest(),
            "wall_seconds": time.monotonic() - started,
        }
    (destination / "eggnog_proteins.dmnd").symlink_to(resource / "eggnog_proteins.dmnd")
    return copied
