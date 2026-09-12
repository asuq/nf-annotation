"""Read ordered annotation input paths without expanding the process argv."""

from __future__ import annotations

from pathlib import Path

from annotation_common import AnnotationError, read_json


def read_path_list(path: Path) -> list[Path]:
    """Validate an explicit JSON path array, retaining its declared order."""
    values = read_json(path)
    if not isinstance(values, list) or any(
        not isinstance(value, str) or not value or "\0" in value for value in values
    ):
        raise AnnotationError(
            f"Expected a JSON array of nonempty path strings without NUL: {path}"
        )
    paths = [Path(value) for value in values]
    if len(paths) != len(set(paths)):
        raise AnnotationError(f"Duplicate path in annotation input list: {path}")
    return paths
