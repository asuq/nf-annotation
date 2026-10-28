"""Explicit, immutable dependencies for reusable ANI reference publications."""

from __future__ import annotations

import os
import re
from pathlib import Path

from annotation_common import (
    AnnotationError,
    digest,
    identity,
    read_json,
    validate_accession,
)

REFERENCE_FILE = "ani_references.json"
REFERENCE_CONTRACT = "nf-annotation-ani-reference-v1"
PORTABLE_CONTRACT = "nf-annotation-v0.4"


def strict_path(root: Path, name: str) -> Path:
    """Resolve an ordinary publication path without following undeclared links."""
    path = Path(name)
    if (
        not name
        or path.is_absolute()
        or any(part in ("", ".", "..") for part in name.split("/"))
    ):
        raise AnnotationError(f"Invalid published artefact path: {name!r}")
    current = root
    for part in path.parts:
        current = current / part
        if current.is_symlink():
            raise AnnotationError(
                f"Published artefact must not be a symlink: {current}"
            )
    return current


class Publication:
    """Resolve only manifest-declared links, flattening physical payload owners."""

    def __init__(
        self,
        root: Path,
        *,
        pending: bool = False,
        document: dict | None = None,
        cache: dict | None = None,
        active: tuple[Path, ...] = (),
    ):
        self.root = root.resolve(strict=True)
        if self.root in active or len(active) >= 64:
            raise AnnotationError(
                f"Circular or excessively deep publication reference: {self.root}"
            )
        cache = {} if cache is None else cache
        manifest_path = strict_path(self.root, "annotation_results.json")
        self.manifest = (
            read_json(manifest_path) if manifest_path.is_file() and not pending else {}
        )
        declared = self.manifest.get("accessions", [])
        if not isinstance(declared, list) or any(
            not isinstance(value, str) for value in declared
        ):
            raise AnnotationError("Invalid publication accession list")
        self.accessions = set(declared)
        contract = self.manifest.get("pipeline_contract")
        self.reference = pending or contract == REFERENCE_CONTRACT
        self.document = None
        self.sources: list[Publication] = []
        self.links: dict[str, Path] = {}
        if not self.reference:
            if contract != PORTABLE_CONTRACT or (self.root / REFERENCE_FILE).exists():
                raise AnnotationError(
                    f"Annotation reuse requires native v0.4 published results: {self.root}"
                )
            return
        if not pending and contract != REFERENCE_CONTRACT:
            raise AnnotationError(
                "Reference publication requires its explicit contract"
            )
        if pending and contract not in (None, REFERENCE_CONTRACT):
            raise AnnotationError(
                "Cannot change an existing publication's storage mode"
            )
        ref_path = strict_path(self.root, REFERENCE_FILE)
        doc = document if document is not None else read_json(ref_path)
        if (
            not isinstance(doc, dict)
            or doc.get("schema_version") != 1
            or doc.get("mode") != "reference"
            or doc.get("reference_id")
            != identity(
                {key: value for key, value in doc.items() if key != "reference_id"}
            )
            or not isinstance(doc.get("sources"), list)
            or not doc["sources"]
            or not isinstance(doc.get("samples"), dict)
            or not doc["samples"]
        ):
            raise AnnotationError("Invalid ANI reference manifest")
        if not pending and self.manifest.get("reference_manifest_sha256") != digest(
            ref_path
        ):
            raise AnnotationError("ANI reference manifest changed")
        if not pending and (
            set(doc["samples"]) != set(self.manifest.get("accessions", []))
            or self.manifest.get("results")
            or self.manifest.get("native_batches")
        ):
            raise AnnotationError(
                "Reference membership or functional-result contract differs"
            )
        self.document = doc
        seen = set()
        for source in doc["sources"]:
            if (
                not isinstance(source, dict)
                or set(source) != {"root", "sha256"}
                or not isinstance(source["root"], str)
                or not Path(source["root"]).is_absolute()
                or not isinstance(source["sha256"], str)
                or not re.fullmatch(r"[a-f0-9]{64}", source["sha256"])
            ):
                raise AnnotationError("Invalid reference source snapshot")
            provider = Path(source["root"]).resolve(strict=True)
            if (
                str(provider) != source["root"]
                or provider in seen
                or provider in (*active, self.root)
            ):
                raise AnnotationError("Duplicate, moved or circular reference source")
            seen.add(provider)
            if provider not in cache:
                cache[provider] = Publication(
                    provider, cache=cache, active=(*active, self.root)
                )
            self.sources.append(cache[provider])
        used = set()
        for accession, sample in doc["samples"].items():
            validate_accession(accession)
            if not isinstance(sample, dict) or set(sample) != {
                "native_source",
                "entries",
                "bundle",
            }:
                raise AnnotationError(f"Invalid reference sample: {accession}")
            provider, index = self.provider(sample["native_source"], accession)
            if provider.reference:
                raise AnnotationError(
                    "Native reference links must be flattened to a self-contained source"
                )
            entries = sample["entries"]
            if (
                not isinstance(entries, list)
                or any(
                    not isinstance(name, str)
                    or name in ("", ".", "..", "annotation")
                    or Path(name).name != name
                    for name in entries
                )
                or len(entries) != len(set(entries))
            ):
                raise AnnotationError(f"Invalid native reference entries: {accession}")
            native = provider.native_root(accession)
            if set(entries) != {
                path.name for path in native.iterdir() if path.name != "annotation"
            }:
                raise AnnotationError(
                    f"Native reference entries changed for {accession}"
                )
            used.add(index)
            for name in entries:
                relative = f"samples/{accession}/{name}"
                self.links[relative] = provider.path(relative)
            bundle = sample["bundle"]
            if bundle is not None:
                if bundle == {"local": True}:
                    strict_path(self.root, f"samples/{accession}/annotation/bundle")
                elif isinstance(bundle, dict) and set(bundle) == {"source"}:
                    owner, index = self.provider(bundle["source"], accession)
                    if owner.reference and owner.document["samples"][accession][
                        "bundle"
                    ] != {"local": True}:
                        raise AnnotationError(
                            "Bundle references must be flattened to their physical owner"
                        )
                    relative = f"samples/{accession}/annotation/bundle"
                    self.links[relative] = owner.path(relative)
                    used.add(index)
                else:
                    raise AnnotationError(f"Invalid bundle reference for {accession}")
        if used != set(range(len(self.sources))):
            raise AnnotationError("Unused reference source snapshot")
        self.check_links()

    def provider(self, index: int, accession: str) -> tuple[Publication, int]:
        if type(index) is not int or not 0 <= index < len(self.sources):
            raise AnnotationError("Invalid reference source index")
        provider = self.sources[index]
        if accession not in provider.accessions:
            raise AnnotationError(f"Reference source does not declare {accession}")
        return provider, index

    def check_links(self) -> None:
        """Reject retargeted, missing and undeclared links in a completed view."""
        samples = strict_path(self.root, "samples")
        if {path.name for path in samples.iterdir()} != set(self.document["samples"]):
            raise AnnotationError(
                "Reference sample directories differ from declared membership"
            )
        for accession, sample in self.document["samples"].items():
            directory = strict_path(self.root, f"samples/{accession}")
            expected = set(sample["entries"]) | (
                {"annotation"} if sample["bundle"] is not None else set()
            )
            if {path.name for path in directory.iterdir()} != expected:
                raise AnnotationError(f"Reference entries differ for {accession}")
            if sample["bundle"] is not None and {
                path.name for path in (directory / "annotation").iterdir()
            } != {"bundle"}:
                raise AnnotationError(
                    f"Unexpected annotation payload in reference view: {accession}"
                )
        for relative, target in self.links.items():
            link = self.root / relative
            if (
                not link.is_symlink()
                or os.readlink(link) != str(target)
                or not target.exists()
            ):
                raise AnnotationError(
                    f"Missing or changed publication reference: {relative}"
                )
        for path in samples.rglob("*"):
            if (
                path.is_symlink()
                and path.relative_to(self.root).as_posix() not in self.links
            ):
                raise AnnotationError(f"Undeclared publication reference: {path}")

    def path(self, name: str) -> Path:
        parts = Path(name).parts
        if (
            not name
            or Path(name).is_absolute()
            or any(part in ("", ".", "..") for part in name.split("/"))
        ):
            raise AnnotationError(f"Invalid published artefact path: {name!r}")
        for end in range(1, len(parts) + 1):
            prefix = "/".join(parts[:end])
            if prefix in self.links:
                link = self.root / prefix
                target = self.links[prefix]
                if not link.is_symlink() or os.readlink(link) != str(target):
                    raise AnnotationError(f"Changed publication reference: {prefix}")
                remainder = "/".join(parts[end:])
                return strict_path(target, remainder) if remainder else target
        return strict_path(self.root, name)

    def native_owner(self, accession: str) -> Publication:
        return (
            self.sources[self.document["samples"][accession]["native_source"]]
            if self.reference
            else self
        )

    def native_root(self, accession: str) -> Path:
        owner = self.native_owner(accession)
        return strict_path(owner.root, f"samples/{accession}")

    def bundle_owner(self, accession: str) -> Publication:
        if self.reference:
            record = self.document["samples"][accession]["bundle"]
            if record and "source" in record:
                return self.sources[record["source"]]
        return self

    def snapshot(self, owner: Publication, current: dict) -> dict:
        if owner.root == self.root:
            if current.get("root") != str(self.root):
                raise AnnotationError("Missing current source fingerprint")
            return {"root": str(self.root), "sha256": current["sha256"]}
        for source in self.document["sources"]:
            if source["root"] == str(owner.root):
                return dict(source)
        raise AnnotationError("Missing reference-provider snapshot")
