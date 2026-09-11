#!/usr/bin/env python3
"""Prepare and validate immutable, tool-consumable annotation resources."""

from __future__ import annotations

import array
import csv
import gzip
import hashlib
import importlib.metadata
import json
import logging
import math
import shutil
import sqlite3
import subprocess
import tarfile
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from download_annotation_resources import (
    ResourceDownloadError,
    checksum_file,
    write_json,
)

LOGGER = logging.getLogger(__name__)
RESOURCE_SCHEMA = 1
RESOURCE_FILE = "annotation_resource.json"
RESOURCE_TOOLS = ("eggnog", "cogclassifier", "pfam", "kofam")


class AnnotationResourceError(RuntimeError):
    """A resource does not satisfy the selected tool's scientific contract."""


def identity(payload: dict[str, Any]) -> str:
    """Hash canonical scientific resource content and preparation settings."""
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def require_files(root: Path, names: list[str] | tuple[str, ...]) -> None:
    """Require non-empty regular inputs within the resource root."""
    for name in names:
        path = root / name
        if path.is_symlink() or not path.is_file() or path.stat().st_size == 0:
            raise AnnotationResourceError(
                f"Missing, empty or linked resource input: {path}"
            )


def file_records(root: Path) -> list[dict[str, Any]]:
    """Inventory each prepared resource file once after successful validation."""
    result = []
    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            raise AnnotationResourceError(
                f"Linked resource input is unsupported: {path}"
            )
        if (
            not path.is_file()
            or path.name == RESOURCE_FILE
            or "validation" in path.relative_to(root).parts
        ):
            continue
        if path.name == "acquisition.json":
            continue
        result.append(
            {
                "path": path.relative_to(root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": checksum_file(path),
            }
        )
    return result


def validate_resource(
    root: Path, component: str, *, verify_checksums: bool = True
) -> dict[str, Any]:
    """Validate the prepared contract; checksum verification belongs in preflight.

    Per-proteome jobs consume the identity already checked once for the run and
    may inspect metadata without reading whole databases again.
    """
    path = root / RESOURCE_FILE
    try:
        payload = json.loads(path.read_text())
    except (OSError, ValueError) as error:
        raise AnnotationResourceError(
            f"Missing or invalid prepared resource record: {path}"
        ) from error
    if (
        payload.get("schema_version") != RESOURCE_SCHEMA
        or payload.get("component") != component
    ):
        raise AnnotationResourceError(
            f"Unsupported or mismatched resource contract: {path}"
        )
    contract = payload.get("contract")
    if not isinstance(contract, dict) or payload.get("resource_id") != identity(
        contract
    ):
        raise AnnotationResourceError(f"Resource identity mismatch: {path}")
    if (
        contract.get("component") != component
        or not contract.get("version")
        or not isinstance(contract.get("settings"), dict)
    ):
        raise AnnotationResourceError(f"Incomplete resource identity: {path}")
    files = contract.get("files")
    if not isinstance(files, list) or not files:
        raise AnnotationResourceError(f"Missing resource inventory: {path}")
    seen = set()
    for entry in files:
        name = entry.get("path")
        if (
            not isinstance(name, str)
            or not name
            or Path(name).is_absolute()
            or any(part in (".", "..") for part in name.split("/"))
            or name in seen
        ):
            raise AnnotationResourceError(
                f"Invalid or duplicate resource path: {name!r}"
            )
        seen.add(name)
        candidate = root / name
        if (
            candidate.is_symlink()
            or not candidate.is_file()
            or candidate.stat().st_size != entry.get("bytes")
        ):
            raise AnnotationResourceError(
                f"Missing or changed resource file: {candidate}"
            )
        if verify_checksums and checksum_file(candidate) != entry.get("sha256"):
            raise AnnotationResourceError(f"Resource checksum mismatch: {candidate}")
    return payload


def run(command: list[str], root: Path, name: str) -> str:
    """Run a required preparation check and preserve command/output evidence."""
    LOGGER.info("Running %s", command[0])
    directory = root / "validation"
    directory.mkdir(exist_ok=True)
    result = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    (directory / f"{name}.log").write_text(result.stdout)
    (directory / f"{name}.command.json").write_text(json.dumps(command) + "\n")
    if result.returncode:
        raise AnnotationResourceError(
            f"{command[0]} exited {result.returncode}; see {directory / (name + '.log')}"
        )
    return result.stdout


def decompress(source: Path, target: Path) -> None:
    """Validate gzip framing while materializing one uncompressed input."""
    with gzip.open(source, "rb") as infile, target.open("wb") as outfile:
        shutil.copyfileobj(infile, outfile, length=8 * 1024 * 1024)


def validate_hmm_profiles(root: Path, profiles: list[Path]) -> None:
    """Pass every selected profile through the native HMM parser in one stream."""
    directory = root / "validation"
    directory.mkdir(exist_ok=True)
    command = ["hmmstat", "-"]
    (directory / "hmmstat.command.json").write_text(
        json.dumps(
            {
                "command": command,
                "stdin_profiles": [
                    path.relative_to(root).as_posix() for path in profiles
                ],
            }
        )
        + "\n"
    )
    with (directory / "hmmstat.log").open("wb") as output:
        with subprocess.Popen(
            command, stdin=subprocess.PIPE, stdout=output, stderr=subprocess.STDOUT
        ) as process:
            try:
                for path in profiles:
                    with path.open("rb") as handle:
                        shutil.copyfileobj(handle, process.stdin, length=1024 * 1024)
                process.stdin.close()
            except (OSError, BrokenPipeError) as error:
                process.terminate()
                raise AnnotationResourceError(
                    f"Could not validate the complete HMM profile subset: {error}"
                ) from error
            if process.wait() != 0:
                raise AnnotationResourceError(
                    f"HMM profile validation failed; see {directory / 'hmmstat.log'}"
                )


def extract(source: Path, target: Path) -> None:
    """Extract data files, rejecting links, devices and escaping archive paths."""
    with tarfile.open(source, "r:gz") as archive:
        members = archive.getmembers()
        seen = set()
        for member in members:
            path = Path(member.name)
            if (
                path.is_absolute()
                or ".." in path.parts
                or not (member.isdir() or member.isfile())
                or path in seen
            ):
                raise AnnotationResourceError(
                    f"Unsafe or duplicate archive entry: {member.name}"
                )
            seen.add(path)
        archive.extractall(target, members=members, filter="data")


def validate_acquisition(source: Path, component: str, version: str) -> dict[str, Any]:
    """Verify the frozen download receipt before resource transformation."""
    try:
        record = json.loads((source / "acquisition.json").read_text())
    except (ValueError, OSError) as error:
        raise AnnotationResourceError(
            f"Missing acquisition record: {source}"
        ) from error
    if record.get("component") != component or record.get("version") != version:
        raise AnnotationResourceError(
            f"Acquisition component/version mismatch: {source}"
        )
    files = record.get("files")
    if not isinstance(files, list) or not files:
        raise AnnotationResourceError(f"Empty acquisition: {source}")
    seen = set()
    for item in files:
        name = item.get("file")
        if not isinstance(name, str) or Path(name).name != name or name in seen:
            raise AnnotationResourceError(f"Invalid acquisition file: {name!r}")
        seen.add(name)
        require_files(source, [name])
        if (source / name).stat().st_size != item.get("bytes") or checksum_file(
            source / name
        ) != item.get("sha256"):
            raise AnnotationResourceError(f"Acquired input changed: {name}")
    return record


def prepare_pfam(source: Path, root: Path) -> dict[str, Any]:
    """Press the complete Pfam library and validate gathering thresholds."""
    for name in ("Pfam-A.hmm", "Pfam-A.hmm.dat", "Pfam-A.clans.tsv", "Pfam.version"):
        decompress(source / f"{name}.gz", root / name)
    models = {}
    current = {}
    with (root / "Pfam-A.hmm").open() as handle:
        for line in handle:
            if line.startswith(("NAME ", "ACC ", "GA ")):
                key, value = line.split(maxsplit=1)
                current[key] = value.strip()
            elif line.strip() == "//":
                try:
                    name, accession = current["NAME"], current["ACC"]
                    scores = [
                        float(value.rstrip(";")) for value in current["GA"].split()
                    ]
                except (KeyError, ValueError) as error:
                    raise AnnotationResourceError(
                        "Pfam model lacks a valid name, accession or gathering thresholds"
                    ) from error
                if (
                    name in models
                    or len(scores) != 2
                    or not all(math.isfinite(score) for score in scores)
                ):
                    raise AnnotationResourceError(
                        f"Invalid or duplicate Pfam model: {name}"
                    )
                models[name] = accession
                current = {}
    if current or not models:
        raise AnnotationResourceError("Empty or truncated Pfam HMM library")
    metadata = {}
    current = {}
    with (root / "Pfam-A.hmm.dat").open() as handle:
        for line in handle:
            if line.startswith(("#=GF ID ", "#=GF AC ")):
                _, key, value = line.split(maxsplit=2)
                current[key] = value.strip()
            elif line.strip() == "//":
                if (
                    "ID" not in current
                    or "AC" not in current
                    or current["ID"] in metadata
                ):
                    raise AnnotationResourceError(
                        "Invalid or duplicate Pfam model metadata"
                    )
                metadata[current["ID"]] = current["AC"]
                current = {}
    if current or models != metadata:
        raise AnnotationResourceError(
            "Pfam library and metadata model identities disagree"
        )
    run(["hmmpress", str(root / "Pfam-A.hmm")], root, "hmmpress")
    run(["hmmstat", str(root / "Pfam-A.hmm")], root, "hmmstat")
    require_files(
        root, ["Pfam-A.hmm" + suffix for suffix in (".h3f", ".h3i", ".h3m", ".h3p")]
    )
    return {
        "threshold": "native gathering thresholds",
        "models": len(models),
        "clan_resolution": "PfamScan 1.6 native",
        "version_record": (root / "Pfam.version").read_text().strip(),
    }


def prepare_kofam(source: Path, root: Path) -> dict[str, Any]:
    """Freeze the native prokaryotic profile subset and adaptive thresholds."""
    extract(source / "profiles.tar.gz", root)
    decompress(source / "ko_list.gz", root / "ko_list")
    hal = root / "profiles/prokaryote.hal"
    require_files(root, ["ko_list", "profiles/prokaryote.hal"])
    selected = []
    for line in hal.read_text().splitlines():
        name = line.strip()
        if not name or name.startswith("#"):
            continue
        if (
            Path(name).name != name
            or not name.startswith("K")
            or not name.endswith(".hmm")
            or name in selected
        ):
            raise AnnotationResourceError(
                f"Invalid or duplicate KOfam profile selection: {name}"
            )
        selected.append(name)
    if not selected:
        raise AnnotationResourceError("Empty prokaryotic profile subset")
    require_files(root / "profiles", selected)
    with (root / "ko_list").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"knum", "threshold", "score_type", "definition"}
        if not required.issubset(reader.fieldnames or []):
            raise AnnotationResourceError("Invalid KOfam ko_list schema")
        kos = {}
        for row in reader:
            if row["knum"] in kos:
                raise AnnotationResourceError(f"Duplicate KO: {row['knum']}")
            if row["threshold"] != "-":
                try:
                    valid = math.isfinite(float(row["threshold"]))
                except ValueError:
                    valid = False
                if not valid or row["score_type"] not in ("full", "domain"):
                    raise AnnotationResourceError(
                        f"Invalid adaptive threshold: {row['knum']}"
                    )
            kos[row["knum"]] = row
    if any(Path(name).stem not in kos for name in selected):
        raise AnnotationResourceError("Selected KOfam profile absent from ko_list")
    # This resource exposes only the approved subset to the search. Full source
    # archive identity remains in acquisition.json for reconstruction.
    for path in (root / "profiles").iterdir():
        if path.name not in selected and path != hal:
            if not path.is_file():
                raise AnnotationResourceError(f"Unexpected KOfam profile entry: {path}")
            path.unlink()
    with (root / "prokaryote_profiles.tsv").open("w") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("profile", "ko", "sha256", "threshold", "score_type"))
        for name in sorted(selected):
            ko = Path(name).stem
            writer.writerow(
                (
                    name,
                    ko,
                    checksum_file(root / "profiles" / name),
                    kos[ko]["threshold"],
                    kos[ko]["score_type"],
                )
            )
    validate_hmm_profiles(root, [root / "profiles" / name for name in selected])
    return {
        "profile_selection": "profiles/prokaryote.hal",
        "profile_count": len(selected),
        "threshold_scale": 1,
        "threshold_policy": "native adaptive thresholds; unavailable thresholds do not pass",
    }


def prepare_cogclassifier(source: Path, root: Path) -> dict[str, Any]:
    """Validate the CDD index, COG mappings and all functional letters."""
    extract(source / "Cog_LE.tar.gz", root)
    decompress(source / "cddid.tbl.gz", root / "cddid.tbl")
    for name in ("cog_definition.tsv", "cog_func_category.tsv"):
        shutil.copyfile(source / name, root / name)
    from cogclassifier.cog import (
        CogCddIdTable,
        CogDefinitionRecord,
        CogFuncCategoryRecord,
    )

    definitions = CogDefinitionRecord(root / "cog_definition.tsv")
    categories = CogFuncCategoryRecord(root / "cog_func_category.tsv")
    mapping = CogCddIdTable(root / "cddid.tbl")
    if (
        not definitions.get_all()
        or not categories.get_all()
        or not mapping._cdd_id2cog_id
    ):
        raise AnnotationResourceError(
            "Empty COG definitions, categories or CDD mapping"
        )
    if len(definitions.get_id_list()) != len(set(definitions.get_id_list())) or len(
        categories.get_letters()
    ) != len(set(categories.get_letters())):
        raise AnnotationResourceError("Duplicate COG definition or functional letter")
    field_errors = []
    for definition in definitions.get_all():
        if not definition.letter or not set(definition.letter) <= set(
            categories.get_letters()
        ):
            field_errors.append(
                {
                    "cog_id": definition.id,
                    "field": "categories",
                    "raw_value": definition.letter,
                    "error": "empty or outside the declared category vocabulary",
                }
            )
    with (root / "definition_field_errors.tsv").open("w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("cog_id", "field", "raw_value", "error"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(field_errors)
    run(["blastdbcmd", "-db", str(root / "Cog"), "-info"], root, "blastdbcmd")
    return {
        "definition_count": len(definitions),
        "category_count": len(categories),
        "definition_field_errors": len(field_errors),
        "database": "Cog",
        "cogclassifier": importlib.metadata.version("cogclassifier"),
    }


def prepare_eggnog(source: Path, root: Path) -> dict[str, Any]:
    """Validate eggNOG v7 schema, indexes and the GO namespace resource."""
    names = (
        "eggnog.db",
        "eggnog_proteins.dmnd",
        "eggnog.db.fieldpresence.bin",
        "eggnog.db.taxids.bin",
        "eggnog.taxa.db",
        "eggnog.taxa.db.traverse.pkl",
        "go-basic.obo",
    )
    require_files(source, names)
    for name in names:
        shutil.copyfile(source / name, root / name)
    from eggnogmapper.annotator.e7.annotate import (
        AnnotationEngine,
        _load_go_namespace_map,
    )
    from eggnogmapper.annotator.e7.db import EggnogDB

    namespace_map = _load_go_namespace_map(str(root / "go-basic.obo"))
    if not namespace_map or set(namespace_map.values()) != {
        "molecular_function",
        "biological_process",
        "cellular_component",
    }:
        raise AnnotationResourceError(
            "GO OBO must provide all three namespaces; flat-GO fallback is forbidden"
        )
    for name in ("eggnog.db", "eggnog.taxa.db"):
        with sqlite3.connect(f"file:{root / name}?mode=ro", uri=True) as database:
            if database.execute("PRAGMA quick_check").fetchall() != [("ok",)]:
                raise AnnotationResourceError(f"SQLite integrity check failed: {name}")
    database = EggnogDB(str(root / "eggnog.db"), load_taxids=False)
    try:
        version = database.get_db_version()
        if version != "7.0.0":
            raise AnnotationResourceError(
                f"Expected eggNOG v7 database; observed {version!r}"
            )
        for query in (
            "SELECT id, name, taxid FROM protein_names LIMIT 1",
            "SELECT * FROM prots LIMIT 1",
            "SELECT protein_id, events FROM event_index LIMIT 1",
            "SELECT i, name, og, og_lca, ev_lca, sp_overlap, side1, side2 FROM sp_events LIMIT 1",
            "SELECT * FROM ogs LIMIT 1",
            "SELECT * FROM ref_terms LIMIT 1",
        ):
            if database.conn.execute(query).fetchone() is None:
                raise AnnotationResourceError(f"Missing annotation data: {query}")
        n_taxids = database.conn.execute(
            "SELECT MAX(id) + 1 FROM protein_names"
        ).fetchone()[0]
        taxids = array.array("i")
        with (root / "eggnog.db.taxids.bin").open("rb") as handle:
            taxids.fromfile(handle, n_taxids)
            if handle.read(1):
                raise AnnotationResourceError("Unexpected taxid cache size")
        # A complete comparison once at preparation prevents a same-size cache
        # from a different database build from silently supplying wrong taxa.
        for row in database.conn.execute("SELECT id, taxid FROM protein_names"):
            if taxids[row["id"]] != row["taxid"]:
                raise AnnotationResourceError(
                    "Taxid cache disagrees with the annotation database"
                )
        engine = AnnotationEngine(
            database, go_obo_path=str(root / "go-basic.obo"), donor_pool="closest"
        )
        n_masks = database.conn.execute("SELECT MAX(id) + 1 FROM prots").fetchone()[0]
        expected_masks = engine.build_field_presence(
            (dict(row) for row in database.conn.execute("SELECT * FROM prots")), n_masks
        )
        masks = array.array("H")
        with (root / "eggnog.db.fieldpresence.bin").open("rb") as handle:
            masks.fromfile(handle, n_masks)
            if handle.read(1):
                raise AnnotationResourceError("Unexpected field-presence cache size")
        if masks != expected_masks:
            raise AnnotationResourceError(
                "Field-presence cache disagrees with the database and GO OBO"
            )
    finally:
        database.conn.close()
    run(
        ["diamond", "dbinfo", "--db", str(root / "eggnog_proteins.dmnd")],
        root,
        "diamond-dbinfo",
    )
    return {
        "database_version": version,
        "go_namespace_count": len(namespace_map),
        "donor_pool": "closest",
        "cache_validation": "complete native field-presence and taxid comparison",
        "eggnog_mapper": importlib.metadata.version("eggnog-mapper"),
    }


def prepare_resource(
    component: str, version: str, source: Path, destination: Path
) -> dict[str, Any]:
    """Build a new resource atomically from verified acquisition inputs."""
    if component not in RESOURCE_TOOLS:
        raise AnnotationResourceError(f"Unsupported annotation resource: {component}")
    if destination.exists():
        payload = validate_resource(destination, component)
        if payload["contract"]["version"] != version:
            raise AnnotationResourceError(
                "Existing resource version differs from the requested version"
            )
        if source != destination:
            source_record = (
                validate_resource(source, component)["acquisition"]
                if (source / RESOURCE_FILE).is_file()
                else validate_acquisition(source, component, version)
            )
            old_files = [
                (entry["file"], entry["sha256"])
                for entry in payload["acquisition"]["files"]
            ]
            new_files = [
                (entry["file"], entry["sha256"]) for entry in source_record["files"]
            ]
            if old_files != new_files:
                raise AnnotationResourceError(
                    "Existing resource has different source content; use a new destination"
                )
        return payload
    if (source / RESOURCE_FILE).is_file():
        payload = validate_resource(source, component)
        if payload["contract"]["version"] != version:
            raise AnnotationResourceError(
                "Prepared source version differs from the requested version"
            )
        temporary = destination.with_name(destination.name + ".preparing")
        shutil.copytree(source, temporary)
        temporary.rename(destination)
        return payload
    acquisition = validate_acquisition(source, component, version)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".preparing")
    temporary.mkdir(exist_ok=False)
    # Failed preparation remains inspectable; it is never called ready or
    # silently deleted. An operator can move it aside before an explicit retry.
    builders = {
        "eggnog": prepare_eggnog,
        "cogclassifier": prepare_cogclassifier,
        "pfam": prepare_pfam,
        "kofam": prepare_kofam,
    }
    try:
        settings = builders[component](source, temporary)
        contract = {
            "component": component,
            "version": version,
            "files": file_records(temporary),
            "settings": settings,
        }
        payload = {
            "schema_version": RESOURCE_SCHEMA,
            "component": component,
            "resource_id": identity(contract),
            "contract": contract,
            "prepared_at": datetime.now(UTC).isoformat(),
            "acquisition": acquisition,
        }
        write_json(temporary / RESOURCE_FILE, payload)
        temporary.rename(destination)
        return payload
    except (OSError, ValueError, sqlite3.Error, ResourceDownloadError) as error:
        raise AnnotationResourceError(
            f"Resource preparation failed for {component}: {error}"
        ) from error
