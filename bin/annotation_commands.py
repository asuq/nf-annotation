"""Explicit native commands and bounded resources for the pinned annotation tools."""

from __future__ import annotations

import math
import shlex
from typing import Any

from annotation_common import AnnotationError


def commands(tool: str, cpus: int, memory_gib: float) -> dict[str, Any]:
    """Return commands using staged paths, independent of resource location.

    DIAMOND memory is bounded using the pinned mapper's documented estimate:
    block_size * 6 + database_GB / index_chunks + threads * 0.5 GB. Reserve
    eight GiB for overhead and annotation; the supported eggNOG 7 index is <24 GB.
    The scheduler remains the authority for enforcing the requested memory.
    """
    if cpus < 1 or not math.isfinite(memory_gib) or memory_gib <= 0:
        raise AnnotationError("Annotation CPU and memory allocations must be positive")
    threads = str(cpus)
    if tool == "eggnog":
        budget = memory_gib - 8 - 24 / 4 - cpus * 0.5
        block = min(8, math.floor(budget / 6 * 10) / 10)
        if block < 0.5:
            raise AnnotationError("eggNOG requires more memory for this CPU allocation")
        steps = [
            [
                "emapper.py",
                "-i",
                "task/input.faa",
                "--itype",
                "proteins",
                "-m",
                "diamond",
                "--dmnd_sensmode",
                "ultra-sensitive",
                "--dmnd_iterate",
                "yes",
                "--dmnd_top",
                "1",
                "--dmnd_block_size",
                str(block),
                "--dmnd_index_chunks",
                "4",
                "--donor_pool",
                "closest",
                "--lazy_cascade",
                "--tax_scope",
                "auto",
                "--target_orthologs",
                "all",
                "--report_orthologs",
                "--annot_evalue",
                "0.001",
                "--pfam_realign",
                "none",
                "--cpu",
                threads,
                "--data_dir",
                "eggnog_data",
                "--temp_dir",
                "scratch",
                "--output_dir",
                "raw",
                "-o",
                "eggnog",
            ]
        ]
        version = [["emapper.py", "--version"], ["diamond", "version"]]
        # The pinned mapper chooses its taxid-cache path based on the directory's
        # writability, even when a validated shipped cache is readable. A small
        # task-local directory of links makes it consume that existing cache;
        # source files stay read-only and no multi-gigabyte data are copied.
        linked_files = (
            "eggnog.db",
            "eggnog_proteins.dmnd",
            "eggnog.db.fieldpresence.bin",
            "eggnog.db.taxids.bin",
            "eggnog.taxa.db",
            "eggnog.taxa.db.traverse.pkl",
            "go-basic.obo",
        )
        steps = [
            ["mkdir", "eggnog_data"],
            *[
                ["ln", "-s", "../resource/" + name, "eggnog_data/" + name]
                for name in linked_files
            ],
            *steps,
        ]
        environment = {"EGGNOG_GO_OBO": "eggnog_data/go-basic.obo"}
    elif tool == "cogclassifier":
        steps = [
            [
                "rpsblast",
                "-query",
                "task/input.faa",
                "-db",
                "resource/Cog",
                "-outfmt",
                "6",
                "-out",
                "raw/rpsblast.tsv",
                "-evalue",
                "0.01",
                "-num_threads",
                threads,
                "-mt_mode",
                "1",
            ],
            [
                "python3",
                "task/classify_cog_hits.py",
                "--input",
                "task/input.faa",
                "--hits",
                "raw/rpsblast.tsv",
                "--database",
                "resource",
                "--output",
                "raw/cogclassifier.native.tsv",
            ],
        ]
        version = [
            ["rpsblast", "-version"],
            [
                "python3",
                "-c",
                "import importlib.metadata; print(importlib.metadata.version('cogclassifier'))",
            ],
        ]
        environment = {}
    elif tool == "pfam":
        steps = [
            [
                "pfam_scan_evidence.pl",
                "--fasta",
                "task/input.faa",
                "--database",
                "resource",
                "--outdir",
                "raw/evidence",
                "--cpus",
                threads,
            ]
        ]
        version = [["hmmscan", "-h"]]
        environment = {}
    elif tool == "kofam":
        steps = [
            [
                "exec_annotation",
                "--profile",
                "resource/profiles/prokaryote.hal",
                "--ko-list",
                "resource/ko_list",
                "--cpu",
                threads,
                "--tmp-dir",
                "scratch",
                "--format",
                "detail-tsv",
                "--no-report-unannotated",
                "-T",
                "1",
                "-o",
                "raw/kofam.tsv",
                "task/input.faa",
            ]
        ]
        version = [["exec_annotation", "--version"], ["hmmsearch", "-h"]]
        environment = {}
    elif tool == "padloc":
        steps = [
            [
                "padloc",
                "--faa",
                "task/input.faa",
                "--gff",
                "task/input.gff",
                "--outdir",
                "raw",
                "--cpu",
                threads,
            ]
        ]
        version = [
            ["padloc", "--version"],
            ["hmmsearch", "-h"],
            ["cat", "/opt/nf-annotation/provenance/padloc-db.txt"],
            ["sha256sum", "/opt/nf-annotation/provenance/padloc-db.tar.gz"],
            ["sha256sum", "-c", "/opt/nf-annotation/provenance/padloc-compiled.sha256"],
        ]
        environment = {}
    else:
        raise AnnotationError(f"Unsupported annotation tool: {tool}")
    return dict(
        steps=steps,
        version_commands=version,
        environment=environment,
        cpus=cpus,
        memory_gib=memory_gib,
    )


def shell_script(command: dict[str, Any]) -> str:
    """Record native failure without preventing independent samples from finishing."""
    exports = [
        f"export {key}={shlex.quote(value)}"
        for key, value in command["environment"].items()
    ]
    versions = [shlex.join(value) for value in command["version_commands"]]
    searches = [shlex.join(value) for value in command["steps"]]
    # The subshell's errexit applies to every native step; its actual exit status
    # is recorded outside it. Nextflow subsequently validates the complete result.
    return "\n".join(
        [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            "mkdir raw scratch",
            "set +e",
            "(",
            "set -e",
            *exports,
            "{",
            *versions,
            "} > raw/versions.txt 2>&1",
            *searches,
            ") > raw/tool.log 2>&1",
            "annotation_exit=$?",
            "printf '%s\\n' \"$annotation_exit\" > raw/exit_code.txt",
            "exit 0",
            "",
        ]
    )
