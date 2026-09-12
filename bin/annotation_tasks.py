"""Validated run identities, task planning, and native v0.4 result reuse."""

from __future__ import annotations

import importlib
import re
import shutil
from pathlib import Path
from typing import Any

from annotation_commands import commands, shell_script
from annotation_common import (
    SCHEMA_VERSION,
    AnnotationError,
    bundle_proteins,
    digest,
    identity,
    read_json,
    selected_tools,
    write_json,
)
from annotation_resources import RESOURCE_FILE, validate_resource
from annotation_result import inventory, native_exit_code, validate_result
from prepare_padloc_input import prepare_gff

PADLOC_RESOURCE = {
    "version": "2.0.0",
    "commit": "7f99b47b75e232b111c18626badb9ac32e8e0b5a",
    "archive_sha256": "de1c9697ee735c7c8a856295ace5614c9c1959b6010cbc074a29125aa9c04fc9",
}
POLICY = {
    "eggnog": {
        "confidence": ["high", "medium"],
        "go": "native_namespace_confidence_no_expansion",
        "batch": "whole_proteome_tool_id_projection",
    },
    "cogclassifier": {
        "best_hit": "native_first",
        "categories": "full_ordered_native_definition",
    },
    "pfam": {
        "acceptance": "native_gathering_and_clan_resolution",
        "counts": "distinct_gene_family",
    },
    "kofam": {"profiles": "prokaryotic", "acceptance": "native_threshold_marker"},
    "padloc": {
        "counts": "distinct_accession_system",
        "coordinates": "canonical_segments",
    },
}


def code_identity(tool: str) -> dict[str, str]:
    """Identify the actual parser and common interpretation code for this tool."""
    root = Path(__file__).resolve().parent
    names = [
        f"normalize_{tool}.py",
        "annotation_normalization.py",
        "annotation_common.py",
        "annotation_tasks.py",
        "annotation_result.py",
    ]
    if tool == "eggnog":
        names.append("annotation_batch_tasks.py")
    return {name: digest(root / name) for name in names}


def runtime_identity(reference: str) -> str:
    """Require a content-addressed Docker reference or hash the actual local SIF."""
    if not isinstance(reference, str) or not reference:
        raise AnnotationError(
            "Every enabled annotation tool requires an immutable container"
        )
    if reference.endswith(".sif"):
        path = Path(reference)
        if not path.is_absolute() or not path.is_file():
            raise AnnotationError(
                f"Container must be an existing absolute SIF path: {reference}"
            )
        return "sif-sha256:" + digest(path)
    if re.fullmatch(r"(?:[^\s]+@)?sha256:[a-f0-9]{64}", reference):
        return reference.rsplit("@", 1)[-1]
    raise AnnotationError(f"Mutable or unsupported annotation container: {reference}")


def planning_code_identity() -> dict[str, str]:
    """Invalidate cached plans when their validation or selection logic changes."""
    root = Path(__file__).resolve().parent
    return {
        name: digest(root / name)
        for name in (
            "prepare_annotation_tasks.py",
            "annotation_commands.py",
            "annotation_path_lists.py",
            "annotation_tasks.py",
            "annotation_source.py",
            "annotation_summary.py",
            "annotation_common.py",
            "annotation_result.py",
            "annotation_batch_tasks.py",
            "eggnog_batches.py",
            "validate_inputs.py",
        )
    }


def preflight(config: dict[str, Any]) -> dict[str, Any]:
    """Check every requested resource once before any proteome search starts."""
    enabled = selected_tools(config["annotation_tools"])
    helper_runtime = runtime_identity(config["helper_container"]) if enabled else None
    result: dict[str, Any] = {
        "schema_version": SCHEMA_VERSION,
        "enabled_tools": list(enabled),
        "planning_code": planning_code_identity(),
        "tools": {},
    }
    for tool in enabled:
        settings = config["tools"][tool]
        runtime = runtime_identity(settings["container"])
        command = commands(tool, settings["cpus"], settings["memory_gib"])
        if tool == "padloc":
            resource_path = None
            resource = PADLOC_RESOURCE
            resource_id = identity(resource)
            manifest_sha = None
        else:
            resource_path = str(Path(settings["resource"]).resolve(strict=True))
            resource = validate_resource(Path(resource_path), tool)
            resource_id = resource["resource_id"]
            manifest_sha = digest(Path(resource_path) / RESOURCE_FILE)
        native_code = {}
        if tool == "eggnog":
            native_code["eggnog_batches.py"] = digest(
                Path(__file__).with_name("eggnog_batches.py")
            )
        elif tool == "cogclassifier":
            native_code["classify_cog_hits.py"] = digest(
                Path(__file__).with_name("classify_cog_hits.py")
            )
        elif tool == "padloc":
            native_code["prepare_padloc_input.py"] = digest(
                Path(__file__).with_name("prepare_padloc_input.py")
            )
        search_method = dict(
            tool=tool,
            runtime_id=runtime,
            resource_id=resource_id,
            command=command,
            native_code=native_code,
        )
        interpretation = dict(
            schema_version=SCHEMA_VERSION,
            runtime_id=helper_runtime,
            code=code_identity(tool),
            policy=POLICY[tool],
        )
        result["tools"][tool] = dict(
            container=settings["container"],
            resource_path=resource_path,
            resource=resource,
            resource_manifest_sha256=manifest_sha,
            search_method=search_method,
            interpretation=interpretation,
            method_id=identity(
                dict(search=search_method, interpretation=interpretation)
            ),
        )
    result["preflight_id"] = identity(result)
    return result


def validate_preflight(record: dict[str, Any]) -> None:
    """Check the receipt and small resource manifests without rehashing databases."""
    if record.get("schema_version") != SCHEMA_VERSION or record.get(
        "preflight_id"
    ) != identity(
        {key: value for key, value in record.items() if key != "preflight_id"}
    ):
        raise AnnotationError("Invalid annotation preflight receipt")
    if record.get("planning_code") != planning_code_identity():
        raise AnnotationError("Annotation planning code changed after preflight")
    for tool in record["enabled_tools"]:
        entry = record["tools"][tool]
        if (
            tool != "padloc"
            and digest(Path(entry["resource_path"]) / RESOURCE_FILE)
            != entry["resource_manifest_sha256"]
        ):
            raise AnnotationError(f"Resource manifest changed after preflight: {tool}")
        if entry["interpretation"]["code"] != code_identity(tool):
            raise AnnotationError(f"Normalizer code changed after preflight: {tool}")


def task_identity(
    bundle: dict[str, Any], tool: str, entry: dict[str, Any]
) -> dict[str, Any]:
    """Separate native-search dependencies from report interpretation dependencies."""
    search = dict(
        input_id=bundle["input_id"],
        genetic_code=bundle["genetic_code"],
        source_genome_sha256=bundle["source_genome_sha256"],
        method=entry["search_method"],
    )
    if tool == "padloc":
        search["coordinate_id"] = bundle["coordinate_id"]
        search["source_gff_sha256"] = bundle["files"]["source.gff"]
    search_id = identity(search)
    return dict(
        schema_version=SCHEMA_VERSION,
        accession=bundle["accession"],
        tool=tool,
        input_id=bundle["input_id"],
        coordinate_id=bundle["coordinate_id"],
        input_proteins=bundle["protein_count"],
        search=search,
        search_fingerprint=search_id,
        normalization_fingerprint=identity(
            dict(search_fingerprint=search_id, interpretation=entry["interpretation"])
        ),
        method_id=entry["method_id"],
        interpretation=entry["interpretation"],
        resource=(
            entry["resource"]
            if tool == "padloc"
            else {
                "resource_id": entry["resource"]["resource_id"],
                "version": entry["resource"]["contract"]["version"],
                "settings": entry["resource"]["contract"]["settings"],
                "manifest_sha256": entry["resource_manifest_sha256"],
            }
        ),
        container=entry["container"],
    )


def plan_task(
    bundle: Path,
    tool: str,
    entry: dict[str, Any],
    outdir: Path,
    previous: Path | None = None,
    *,
    bundle_data: tuple[dict[str, Any], list[dict[str, str]]],
) -> dict[str, Any]:
    """Plan a tool from the bundle validated once in this planner invocation."""
    manifest, proteins = bundle_data
    task = task_identity(manifest, tool, entry)
    task.update(action="run", reason="no_compatible_native_result")
    if previous is not None and (previous / "result.json").is_file():
        old = read_json(previous / "result.json")
        if (
            old.get("status") == "success"
            and old.get("search_fingerprint") == task["search_fingerprint"]
        ):
            reuse_normalized = (
                old.get("normalization_fingerprint")
                == task["normalization_fingerprint"]
            )
            validate_result(previous, normalized=reuse_normalized)
            task.update(
                action="reuse" if reuse_normalized else "renormalize",
                reason="matching_search_and_interpretation"
                if reuse_normalized
                else "interpretation_changed",
            )
    outdir.mkdir(parents=True)
    if task["action"] == "run":
        shutil.copyfile(
            bundle / ("source.faa" if tool == "padloc" else "proteins.faa"),
            outdir / "input.faa",
        )
        if tool == "padloc":
            prepare_gff(bundle, proteins, outdir / "input.gff")
        if tool == "cogclassifier":
            shutil.copyfile(
                Path(__file__).with_name("classify_cog_hits.py"),
                outdir / "classify_cog_hits.py",
            )
        (outdir / "run.sh").write_text(shell_script(entry["search_method"]["command"]))
    else:
        shutil.copytree(previous / "raw", outdir / "previous" / "raw")
        if task["action"] == "reuse":
            shutil.copytree(previous / "normalized", outdir / "previous" / "normalized")
            shutil.copyfile(
                previous / "result.json", outdir / "previous" / "result.json"
            )
    write_json(outdir / "task.json", task)
    return task


def normalize_task(
    taskdir: Path, raw: Path, bundle: Path, resource: Path | None, outdir: Path
) -> dict[str, Any]:
    """Validate a native execution and publish a result even when that tool fails."""
    task = read_json(taskdir / "task.json")
    manifest, proteins = bundle_proteins(bundle)
    if (
        manifest["input_id"] != task["input_id"]
        or manifest["accession"] != task["accession"]
    ):
        raise AnnotationError("Task/bundle identity mismatch")
    if task["interpretation"]["code"] != code_identity(task["tool"]):
        raise AnnotationError("Normalizer changed after task planning")
    outdir.mkdir()
    shutil.copytree(raw, outdir / "raw")
    record = dict(
        task,
        status="failed",
        reason="native_output_not_validated",
        exit_code=None,
        reported_hits=None,
        accepted_proteins=None,
        evidence=None,
        normalized_files={},
    )
    try:
        exit_code = native_exit_code(raw)
        record["exit_code"] = exit_code
        if record["exit_code"] != 0:
            raise AnnotationError(
                f"Native {task['tool']} exited with status {exit_code}"
            )
        tool = task["tool"]
        if tool == "padloc":
            versions = (raw / "versions.txt").read_text()
            if (
                f"commit={PADLOC_RESOURCE['commit']}\n" not in versions
                or f"{PADLOC_RESOURCE['archive_sha256']}  /opt/nf-annotation/provenance/padloc-db.tar.gz\n"
                not in versions
            ):
                raise AnnotationError(
                    "PADLOC bundled resource differs from the pinned contract"
                )
        else:
            prepared = read_json(resource / RESOURCE_FILE)
            if prepared.get("resource_id") != task["search"]["method"]["resource_id"]:
                raise AnnotationError("Task/resource identity mismatch")
        parser = importlib.import_module("normalize_" + tool)
        evidence = parser.normalize(
            raw / "evidence" if tool == "pfam" else raw,
            proteins,
            bundle if tool == "padloc" else resource,
        )
        summary = evidence.publish(outdir / "normalized")
        record.update(
            status="success",
            reason="validated_native_output",
            evidence=summary,
            reported_hits=evidence.reported_hits,
            accepted_proteins=len(evidence.accepted),
            normalized_files=inventory(outdir / "normalized"),
        )
    except (AnnotationError, OSError) as error:
        record["reason"] = str(error)
        (outdir / "validation_error.txt").write_text(str(error) + "\n")
    record["raw_files"] = inventory(outdir / "raw")
    record["result_id"] = identity(record)
    write_json(outdir / "result.json", record)
    return record
