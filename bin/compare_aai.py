#!/usr/bin/env python3
"""Run native FastAAI comparisons once and publish one record per unordered pair."""

from __future__ import annotations
import argparse
import shutil
import subprocess
import tempfile
import time
from pathlib import Path
from annotation_common import AnnotationError, read_json, write_json
from annotation_path_lists import read_path_list


def compare(records, profiles, output, threads):
    start = time.monotonic()
    records = sorted(records, key=lambda r: r["accession"])
    ready = [r for r in records if r["status"] == "ready"]
    by_id = {read_json(p / "profile.json")["profile_id"]: p for p in profiles}
    output.mkdir()
    aliases = {f"S{i:08d}.faa": r for i, r in enumerate(ready)}
    expected = len(ready) * (len(ready) - 1) // 2
    count = 0
    phase_seconds = {}
    from aai_evidence import create, offset
    import numpy as np

    runtime = read_json(Path("/opt/nf-annotation/provenance/fastaai.json"))
    if runtime.get("export") != "native_condensed_v1":
        raise AnnotationError("FastAAI runtime must provide the v0.5 binary export")
    arrays = create(output, [r["accession"] for r in ready])
    if len(ready) > 1:
        with tempfile.TemporaryDirectory(prefix="fastaai-query-") as directory:
            work = Path(directory)
            proteins = work / "proteins"
            hmms = work / "hmms"
            proteins.mkdir()
            hmms.mkdir()
            for alias, row in aliases.items():
                profile = by_id[Path(row["profile"]).name]
                (proteins / alias).symlink_to((profile / "markers.faa").resolve())
                (hmms / (Path(alias).stem + ".hmm")).symlink_to(
                    (profile / "markers.hmm").resolve()
                )
            native_out = work / "native"
            for phase, command in zip(
                ("database", "query"),
                (
                    [
                        "fastaai",
                        "build_db",
                        "--proteins",
                        str(proteins),
                        "--hmms",
                        str(hmms),
                        "--database",
                        "cohort.db",
                        "--output",
                        str(native_out),
                        "--threads",
                        str(threads),
                    ],
                    [
                        "fastaai",
                        "db_query",
                        "--query",
                        str(native_out / "database/cohort.db"),
                        "--target",
                        str(native_out / "database/cohort.db"),
                        "--output",
                        str(native_out),
                        "--threads",
                        str(threads),
                        "--output_style",
                        "binary",
                    ],
                ),
                strict=True,
            ):
                if phase == "query":
                    import sqlite3

                    with sqlite3.connect(native_out / "database/cohort.db") as database:
                        observed = database.execute(
                            "SELECT genome FROM genome_index ORDER BY gen_id"
                        ).fetchall()
                    if [row[0] for row in observed] != list(aliases):
                        raise AnnotationError(
                            "Native database changed genome identities or ordering"
                        )
                phase_start = time.monotonic()
                with (output / "native.log").open("a") as log:
                    subprocess.run(
                        command, check=True, stdout=log, stderr=subprocess.STDOUT
                    )
                phase_seconds[phase] = time.monotonic() - phase_start
            shutil.copyfile(
                native_out / "logs/FastAAI_preprocessing_log.txt",
                output / "database_inputs.tsv",
            )

            for i in range(len(ready) - 1):
                with np.load(
                    native_out / "results" / f"{i:08d}.npz", allow_pickle=False
                ) as row:
                    length = len(ready) - i - 1
                    if any(row[name].shape != (length,) for name in arrays):
                        raise AnnotationError("Incomplete native FastAAI row")
                    row_start = offset(len(ready), i, i + 1)
                    for name in arrays:
                        arrays[name][row_start : row_start + length] = row[name]
                    count += length
    for array in arrays.values():
        if hasattr(array, "flush"):
            array.flush()
    if count != expected:
        raise AnnotationError("Incomplete FastAAI pair coverage")
    write_json(
        output / "comparison.json",
        dict(
            accessions=[r["accession"] for r in ready],
            pairs=count,
            runtime=runtime,
        ),
    )
    write_json(
        output / "metrics.json",
        dict(
            seconds=time.monotonic() - start,
            pairs=count,
            workers=threads,
            phases=phase_seconds,
        ),
    )
    (output / "native.log").touch(exist_ok=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("records", "profiles", "output"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--threads", type=int, required=True)
    args = p.parse_args()
    try:
        compare(
            read_json(args.records),
            read_path_list(args.profiles),
            args.output,
            args.threads,
        )
    except (
        AnnotationError,
        OSError,
        ValueError,
        KeyError,
        subprocess.CalledProcessError,
    ) as e:
        p.exit(1, str(e) + "\n")


if __name__ == "__main__":
    main()
