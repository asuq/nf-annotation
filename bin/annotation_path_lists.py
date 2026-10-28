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


def read_publication_paths(
    path: Path, publication: Path | None = None, mode: str = "self_contained"
) -> list[Path]:
    """Bind exact bundle paths to one mounted publication without thousands of mounts."""
    paths = read_path_list(path)
    if publication is None:
        return paths
    root = publication.resolve(strict=True)
    if mode not in ("self_contained", "reference"):
        raise AnnotationError("Unknown publication mode")
    view = None
    if mode == "reference":
        from publication_references import Publication

        view = Publication(root, pending=True)
    checked = []
    for bundle in paths:
        declared = publication.absolute()
        lexical_root = next(
            (
                candidate
                for candidate in (declared, root)
                if bundle.is_relative_to(candidate)
            ),
            None,
        )
        if not bundle.is_absolute() or ".." in bundle.parts or lexical_root is None:
            raise AnnotationError("Bundle path escapes its declared publication")
        if view is not None:
            resolved = view.path(bundle.relative_to(lexical_root).as_posix())
        else:
            resolved = bundle.resolve(strict=True)
            if not resolved.is_relative_to(root) or not resolved.is_dir():
                raise AnnotationError("Bundle path escapes its declared publication")
            current = lexical_root
            for part in bundle.relative_to(lexical_root).parts:
                current = current / part
                if current.is_symlink():
                    raise AnnotationError(
                        "Publication bundle must not traverse symlinks"
                    )
        if not (resolved / "bundle.json").is_file():
            raise AnnotationError("Missing declared publication bundle")
        checked.append(resolved)
    if len(checked) != len(set(checked)):
        raise AnnotationError("Duplicate resolved publication bundle")
    return checked
