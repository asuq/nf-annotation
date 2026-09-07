"""Input-driven external-tool fixtures for cohort-update integration tests.

Only expensive tool execution is synthetic. The integration harness runs the
production validation, summarisation, clustering and reporting scripts.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Data import CodonTable

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import summarise_codetta  # noqa: E402


def write(path: str | Path, content: str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


def table(path: str | Path, header: list[str], rows: list[list[str]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def fake_seqtk(argv: list[str]) -> int:
    """Exercise the real staging/statistics wrappers on small FASTA records."""
    if not argv:
        print("Version: synthetic-test")
        return 0
    records = list(SeqIO.parse(argv[-1], "fasta"))
    if argv[:2] == ["seq", "-A"]:
        SeqIO.write(records, sys.stdout, "fasta")
    elif argv[0] == "comp":
        for record in records:
            sequence = str(record.seq).upper()
            print(
                "\t".join(
                    [
                        record.id,
                        str(len(sequence)),
                        *(str(sequence.count(base)) for base in "ACGT"),
                    ]
                )
            )
    else:
        raise ValueError(f"Unsupported synthetic seqtk invocation: {argv}")
    return 0


def fake_fastani(argv: list[str]) -> int:
    """Exercise the real FastANI wrapper, including its empty-cohort guard."""
    if argv == ["--version"]:
        print("FastANI synthetic-test")
        return 0
    parser = argparse.ArgumentParser()
    parser.add_argument("--rl", type=Path, required=True)
    parser.add_argument("--ql", type=Path, required=True)
    parser.add_argument("--matrix", action="store_true")
    parser.add_argument("-t")
    parser.add_argument("-o", type=Path, required=True)
    args = parser.parse_args(argv)
    paths = args.ql.read_text().splitlines()
    if not paths:
        raise ValueError("The workflow must skip FastANI for an empty eligible cohort.")
    if args.rl.read_text().splitlines() != paths:
        raise ValueError("This fixture expects all-vs-all FastANI inputs.")
    sequences = [
        "".join(str(record.seq) for record in SeqIO.parse(path, "fasta"))
        for path in paths
    ]
    matrix = [str(len(paths))]
    for index, path in enumerate(paths):
        sequence = sequences[index]
        values = [
            f"{100 * sum(a == b for a, b in zip(sequence, other)) / max(len(sequence), len(other)):.6f}"
            for other in sequences[:index]
        ]
        matrix.append("\t".join([path, *values]))
    write(str(args.o) + ".matrix", "\n".join(matrix) + "\n")
    write(args.o, "")
    print("Synthetic test identities; no real FastANI execution.")
    return 0


def make_stub(args: argparse.Namespace) -> None:
    """Write one process's declared outputs from the actual process inputs."""
    tool = args.tool
    internal_id = args.internal_id
    accession = args.accession
    if tool == "barrnap":
        record = next(SeqIO.parse(args.genome, "fasta"))
        sequence = str(record.seq)[:1500]
        partial = "partial" in accession.lower()
        attributes = "Name=16S_rRNA" + (";partial=01" if partial else "")
        write(
            "rrna.gff",
            f"{record.id}\tbarrnap\trRNA\t1\t{len(sequence)}\t1e-20\t+\t.\t{attributes}\n",
        )
        write("rrna.fa", f">16S_rRNA::{record.id}:0-{len(sequence)}(+)\n{sequence}\n")
        write("barrnap.log", "exit_code=0\n")
    elif tool == "checkm2":
        Path(f"checkm2_gcode{args.translation_table}").mkdir()
        failed = "nogcode" in accession.lower()
        completeness = (
            "30"
            if "low" in accession.lower()
            else ("95" if args.translation_table == "4" else "82")
        )
        if failed:
            write("quality_report.tsv", "")
        else:
            table(
                "quality_report.tsv",
                [
                    "Name",
                    "Completeness",
                    "Contamination",
                    "Coding_Density",
                    "Average_Gene_Length",
                    "Total_Coding_Sequences",
                ],
                [[internal_id, completeness, "1", "0.9", "900", "800"]],
            )
        write("checkm2.log", f"exit_code={1 if failed else 0}\n")
    elif tool == "busco":
        Path(f"busco_{args.lineage}").mkdir()
        write(
            "short_summary.json",
            json.dumps(
                {
                    "results": {
                        "C": 98.0,
                        "S": 98.0,
                        "D": 0.0,
                        "F": 1.0,
                        "M": 1.0,
                        "n": 200,
                    }
                }
            ),
        )
        write("busco.log", "exit_code=0\n")
    elif tool == "codetta":
        genetic_code = summarise_codetta.build_ncbi_table_string(
            CodonTable.generic_by_id[4]
        )
        lines = [
            "# Codon inferences",
            "# codon inference std code diff? N aligned N used",
        ]
        lines.extend(
            f"{codon}\t{amino_acid}\tX\t.\t1\t1"
            for codon, amino_acid in zip(
                summarise_codetta.CODON_ORDER, genetic_code, strict=True
            )
        )
        write("codetta/codetta_inference.txt", "\n".join([*lines, "#", ""]))
        write("codetta/codetta.log", "exit_code=0\n")
    elif tool == "prokka":
        Path("prokka").mkdir()
        failed = "failed" in accession.lower()
        for suffix, content in {
            "gff": f"##gff-version 3\n{internal_id}\tProkka\tCDS\t1\t30\t.\t+\t0\tID={internal_id}_1\n",
            "faa": f">{internal_id}_1\nMAAAAAAAAA\n",
            "gbk": f"LOCUS       {internal_id}\n//\n",
        }.items():
            write(f"prokka.{suffix}", "" if failed else content)
        write("prokka.log", f"exit_code={1 if failed else 0}\n")
    elif tool == "ccfinder":
        Path("ccfinder").mkdir()
        write(
            "result.json",
            json.dumps(
                {"Sequences": [{"Id": internal_id, "Length": 2000, "Crisprs": []}]}
            ),
        )
        write("ccfinder.log", "exit_code=0\n")
    elif tool == "padloc":
        Path("padloc").mkdir()
        write("padloc/results.tsv", "system\n")
        write("padloc.log", "exit_code=0\n")
    elif tool == "eggnog":
        Path("eggnog").mkdir()
        table(
            "eggnog_annotations.tsv",
            ["query", "seed_ortholog", "evalue", "score"],
            [[f"{internal_id}_1", "test_ortholog", "1e-20", "200"]],
        )
        write("eggnog.log", "exit_code=0\n")
    else:
        raise ValueError(f"Unknown synthetic process: {tool}")
    write("versions.yml", f'"{args.process_name}":\n  {tool}: "synthetic-test"\n')


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tool")
    parser.add_argument("--process-name", required=True)
    parser.add_argument("--accession", default="")
    parser.add_argument("--internal-id", default="")
    parser.add_argument("--genome", type=Path)
    parser.add_argument("--translation-table")
    parser.add_argument("--lineage")
    make_stub(parser.parse_args())
    return 0


if __name__ == "__main__":
    sys.exit(main())
