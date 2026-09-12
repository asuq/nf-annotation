"""Synthetic native phase records for focused eggNOG archive and parser tests."""

from __future__ import annotations

import csv
import shutil
from pathlib import Path

from annotation_common import write_json
from eggnog_batches import validate_batch
from eggnog_native import (
    EXECUTION_SCHEMA,
    SEED_NATIVE_COLUMNS,
    member_name,
    partition_rows,
    query_lookup,
    read_native_seeds,
    seed_payload,
    stage_plan,
    write_partitions,
)
from normalize_eggnog import FIELDS, GO_HEADER, HEADER


def native_annotations(path: Path, rows: list[dict[str, str]], count: int) -> None:
    """Write explicitly synthetic native-format annotations with complete markers."""
    with path.open("w", newline="") as handle:
        handle.write("## confidence codes: h=high m=medium l=low -=not annotated\n")
        handle.write("## confidence field order: " + " ".join(FIELDS) + "\n")
        handle.write("#" + "\t".join(HEADER) + "\n")
        for row in rows:
            handle.write("\t".join(row[name] for name in HEADER) + "\n")
        handle.write(f"## {count} queries scanned\n")


def phased_raw(inputs: Path, raw: Path, source: Path | None = None) -> None:
    """Build a synthetic shared search and one genuine-format member table each."""
    batch, proteins = validate_batch(inputs)
    raw.mkdir(exist_ok=True)
    search = raw / "search"
    search.mkdir()
    annotations, namespaces = [], []
    if source is None:
        (search / "eggnog.emapper.seed_orthologs").write_text(
            "#" + "\t".join(SEED_NATIVE_COLUMNS) + "\n## 0 queries scanned\n"
        )
    else:
        shutil.copyfile(
            source / "eggnog.emapper.seed_orthologs",
            search / "eggnog.emapper.seed_orthologs",
        )
        with (source / "eggnog.emapper.annotations").open() as handle:
            lines = [line for line in handle if not line.startswith("##")]
        lines[0] = lines[0].removeprefix("#")
        annotations = list(csv.DictReader(lines, delimiter="\t"))
        with (source / "eggnog.emapper.annotations.go_namespaces.tsv").open() as handle:
            namespaces = list(csv.DictReader(handle, delimiter="\t"))
    rows = read_native_seeds(
        search / "eggnog.emapper.seed_orthologs", query_lookup(proteins)
    )
    partitions = partition_rows(batch, proteins, rows)
    write_partitions(inputs, raw, batch, proteins)
    plans = stage_plan(batch)
    (search / "exit_code.txt").write_text("0\n")
    (raw / "annotations").mkdir()
    for index, member in enumerate(batch["members"]):
        name, accession = member_name(index), member["accession"]
        directory = raw / "annotations" / name
        directory.mkdir()
        identifiers = {
            row["tool_id"] for row in proteins if row["accession"] == accession
        }
        selected = partitions[accession]
        native_annotations(
            directory / "eggnog.emapper.annotations",
            [row for row in annotations if row["query"] in identifiers],
            len(selected),
        )
        with (directory / "eggnog.emapper.annotations.go_namespaces.tsv").open(
            "w", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle, fieldnames=GO_HEADER, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(row for row in namespaces if row["query"] in identifiers)
        (directory / "exit_code.txt").write_text("0\n")
        (raw / "derived" / name / "seeds.tsv.sorted").write_bytes(
            seed_payload(
                sorted(selected, key=lambda row: (int(row["sseqid"]), row["qseqid"]))
            )
        )
    (raw / "exit_code.txt").write_text("0\n")
    write_json(
        raw / "execution.json",
        {
            "schema": EXECUTION_SCHEMA,
            "batch_id": batch["batch_id"],
            "stages": [dict(plan, exit_code=0, wall_seconds=0.01) for plan in plans],
        },
    )
