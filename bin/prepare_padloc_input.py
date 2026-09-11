"""Align PADLOC's GFF protein_id attribute with the validated native FAA IDs."""

from pathlib import Path
from urllib.parse import unquote

from annotation_common import AnnotationError, read_tsv


def prepare_gff(bundle: Path, proteins: list[dict[str, str]], output: Path) -> None:
    """Preserve feature order and coordinates while exporting the checked ID join."""
    by_gene = {row["gene_id"]: row["protein_id"] for row in proteins}
    expected = {}
    for row in read_tsv(bundle / "gene_coordinates.tsv"):
        key = tuple(
            row[name]
            for name in ("contig_id", "start", "end", "strand", "phase", "attributes")
        )
        if key in expected:
            raise AnnotationError("Duplicate canonical PADLOC coordinate segment")
        expected[key] = by_gene[row["gene_id"]]
    with (bundle / "source.gff").open() as source, output.open("w") as dest:
        for line in source:
            if line.rstrip("\r\n") == "##FASTA":
                break
            if not line.strip() or line.startswith("#"):
                dest.write(line)
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise AnnotationError("Invalid PADLOC source GFF feature")
            if fields[2] == "CDS":
                key = tuple(fields[index] for index in (0, 3, 4, 6, 7, 8))
                if key not in expected:
                    raise AnnotationError(
                        "PADLOC source CDS differs from canonical coordinates"
                    )
                protein = expected.pop(key)
                if any(char in protein for char in ";=%\t\r\n"):
                    raise AnnotationError(
                        "PADLOC cannot preserve this native protein ID in GFF"
                    )
                attributes = [
                    value
                    for value in fields[8].split(";")
                    if value and unquote(value.split("=", 1)[0]) != "protein_id"
                ]
                fields[8] = ";".join([*attributes, f"protein_id={protein}"])
            dest.write("\t".join(fields) + "\n")
    if expected:
        raise AnnotationError("PADLOC input omits canonical CDS segments")
