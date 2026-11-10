"""Bind published cohort evidence to the existing annotation manifest."""

from pathlib import Path
from annotation_common import AnnotationError, digest, read_json
from publication_references import strict_path


def inventory(entries):
    result = {}
    for entry in entries:
        prefix = entry["prefix"]
        if not (prefix == "cohort" or prefix.startswith("cohort/")) or any(
            p in ("", ".", "..") for p in prefix.split("/")
        ):
            raise AnnotationError("Invalid cohort evidence path")
        root = Path(entry["path"]).resolve(strict=True)
        paths = sorted(root.rglob("*")) if root.is_dir() else [root]
        for path in paths:
            if path.is_symlink():
                raise AnnotationError("Cohort evidence must not contain symlinks")
            if path.is_dir():
                continue
            if not path.is_file():
                raise AnnotationError("Cohort evidence must be an ordinary file")
            name = (
                prefix + "/" + path.relative_to(root).as_posix()
                if root.is_dir()
                else prefix
            )
            if name in result:
                raise AnnotationError("Duplicate cohort evidence path: " + name)
            result[name] = digest(path)
    return result


def verify(root, records, checksum=digest):
    if not isinstance(records, dict):
        raise AnnotationError("Invalid cohort evidence manifest")
    for name, expected in records.items():
        if (
            not name.startswith("cohort/")
            or checksum(strict_path(root, name)) != expected
        ):
            raise AnnotationError("Published cohort evidence changed: " + name)


def read_entries(path):
    return read_json(path) if path is not None else []
