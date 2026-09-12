"""Identity and coordinate controls for reusable annotation bundles."""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from annotation_common import AnnotationError, bundle_proteins, read_tsv
from prepare_annotation_bundle import build_bundle


class AnnotationBundleTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.genome = self.root / "genome.fa"
        self.faa = self.root / "proteins.fa"
        self.gff = self.root / "genes.gff"
        self.gbk = self.root / "genes.gbk"
        self.genome.write_text(">original_contig\nATGTGAGCTTAA\n")
        self.faa.write_text(">gene_1 a description with pipes|and_underscores\nMWA\n")
        self.gff.write_text(
            "##gff-version 3\nrenamed\tProkka\tCDS\t1\t12\t.\t+\t0\tID=gene_1;locus_tag=gene_1\n##FASTA\n>renamed\nATGTGAGCTTAA\n"
        )
        self.write_gbk()

    def write_gbk(self, *, location=None, translation="MWA", gcode=4, topology=None):
        record = SeqRecord(
            Seq("ATGTGAGCTTAA"), id="renamed", name="renamed", description="fixture"
        )
        record.annotations["molecule_type"] = "DNA"
        if topology:
            record.annotations["topology"] = topology
        record.features = [
            SeqFeature(
                location if location is not None else SimpleLocation(0, 12, strand=1),
                type="CDS",
                qualifiers={
                    "locus_tag": ["gene_1"],
                    "translation": [translation],
                    "transl_table": [str(gcode)],
                    "codon_start": ["1"],
                },
            )
        ]
        SeqIO.write([record], self.gbk, "genbank")

    def build(self, accession="SAMPLE_A", name="bundle", gcode=4):
        output = self.root / name
        build_bundle(
            accession, gcode, self.genome, self.faa, self.gff, self.gbk, output
        )
        return output

    def test_descriptions_do_not_become_ids_and_sequences_are_unchanged(self):
        output = self.build()
        manifest, proteins = bundle_proteins(output)
        self.assertEqual(manifest["genetic_code"], 4)
        self.assertEqual(proteins[0]["protein_id"], "gene_1")
        self.assertEqual(proteins[0]["source_contig_id"], "original_contig")
        self.assertEqual(proteins[0]["gene_id"], "SAMPLE_A::gene_1")
        self.assertEqual(
            (output / "proteins.faa").read_text(), f">{proteins[0]['tool_id']}\nMWA\n"
        )
        self.assertEqual(proteins[0]["topology"], "NA")
        for path, name in (
            (self.genome, "genome.fasta"),
            (self.faa, "source.faa"),
            (self.gff, "source.gff"),
            (self.gbk, "source.gbk"),
        ):
            self.assertEqual(path.read_bytes(), (output / name).read_bytes())

    def test_repeated_locus_tag_in_different_samples_has_distinct_identity(self):
        a = bundle_proteins(self.build())[1][0]
        b = bundle_proteins(self.build("SAMPLE_B", "bundle_b"))[1][0]
        self.assertNotEqual(a["tool_id"], b["tool_id"])
        self.assertNotEqual(a["gene_id"], b["gene_id"])

    def test_orphan_and_duplicate_proteins_fail(self):
        self.faa.write_text(">other\nMWA\n")
        with self.assertRaises(AnnotationError):
            self.build()
        self.faa.write_text(">gene_1\nMWA\n>gene_1\nMWA\n")
        with self.assertRaises(AnnotationError):
            self.build()

    def test_changed_bundle_bytes_are_rejected(self):
        output = self.build()
        (output / "proteins.faa").write_text(">unknown\nMWA\n")
        with self.assertRaisesRegex(AnnotationError, "content mismatch"):
            bundle_proteins(output)

    def test_negative_strand_coordinates_are_preserved(self):
        # Reverse complement encodes the same code-4 protein on the opposite strand.
        self.genome.write_text(">original_contig\nTTAAGCTCACAT\n")
        self.gff.write_text(
            self.gff.read_text().replace("ATGTGAGCTTAA", "TTAAGCTCACAT")
        )
        self.gff.write_text(self.gff.read_text().replace("\t+\t", "\t-\t"))
        self.write_gbk(location=SimpleLocation(0, 12, strand=-1))
        record = SeqIO.read(self.gbk, "genbank")
        record.seq = Seq("TTAAGCTCACAT")
        SeqIO.write([record], self.gbk, "genbank")
        output = self.build()
        row = read_tsv(output / "gene_coordinates.tsv")[0]
        self.assertEqual((row["start"], row["end"], row["strand"]), ("1", "12", "-"))

    def test_ambiguous_contig_mapping_and_out_of_bounds_features_fail(self):
        self.genome.write_text(
            self.genome.read_text() + ">duplicate_sequence\nATGTGAGCTTAA\n"
        )
        with self.assertRaisesRegex(AnnotationError, "uniquely map"):
            self.build()
        self.genome.write_text(">original_contig\nATGTGAGCTTAA\n")
        self.gff.write_text(self.gff.read_text().replace("\t12\t", "\t13\t"))
        with self.assertRaisesRegex(AnnotationError, "outside source"):
            self.build()

    def test_native_iupac_masking_preserves_sources_and_records_provenance(self):
        source = "ATGTGAYCTTAARYSWKMBDHV"
        native = "ATGTGANCTTAA" + "N" * 10
        self.genome.write_text(f">original_contig\n{source.lower()}\n")
        self.faa.write_text(">gene_1\nMWX\n")
        self.gff.write_text(self.gff.read_text().replace("ATGTGAGCTTAA", native))
        self.write_gbk(translation="MWX")
        record = SeqIO.read(self.gbk, "genbank")
        record.seq = Seq(native)
        SeqIO.write([record], self.gbk, "genbank")
        output = self.build()
        manifest, proteins = bundle_proteins(output)
        self.assertEqual(
            manifest["source_sequence_policy"], "prokka_uppercase_iupac_to_n"
        )
        provenance = manifest["contig_sequence_provenance"][0]
        self.assertEqual(provenance["masked_iupac_bases"], 11)
        self.assertEqual(provenance["length"], 22)
        self.assertEqual(provenance["source_contig_id"], "original_contig")
        self.assertNotEqual(
            provenance["source_sequence_sha256"], provenance["native_sequence_sha256"]
        )
        self.assertEqual(proteins[0]["source_contig_id"], "original_contig")
        self.assertEqual((output / "source.faa").read_bytes(), self.faa.read_bytes())
        self.assertEqual(
            (output / "genome.fasta").read_bytes(), self.genome.read_bytes()
        )
        self.assertEqual((output / "source.gff").read_bytes(), self.gff.read_bytes())
        self.assertEqual((output / "source.gbk").read_bytes(), self.gbk.read_bytes())

    def test_iupac_masking_does_not_resolve_nonunique_contigs(self):
        self.genome.write_text(">first\nATGTGAYCTTAA\n>second\nATGTGANCTTAA\n")
        self.gff.write_text(
            self.gff.read_text().replace("ATGTGAGCTTAA", "ATGTGANCTTAA")
        )
        with self.assertRaisesRegex(AnnotationError, "uniquely map"):
            self.build()

    def test_source_substitutions_indels_and_ambiguity_resolution_fail(self):
        for native in ("ATGAGANCTTAA", "ATGTGANCTTA", "ATGTGACCTTAA"):
            with self.subTest(native=native):
                # Even resolving source Y to a possible base is not Prokka's rule.
                self.genome.write_text(">original_contig\nATGTGAYCTTAA\n")
                self.gff.write_text(
                    "##gff-version 3\nrenamed\tProkka\tCDS\t1\t12\t.\t+\t0\tID=gene_1\n"
                    f"##FASTA\n>renamed\n{native}\n"
                )
                with self.assertRaisesRegex(AnnotationError, "uniquely map"):
                    self.build()

    def test_empty_protein_input_fails(self):
        self.faa.write_text("")
        with self.assertRaisesRegex(AnnotationError, "empty_proteome"):
            self.build()

    def test_invalid_sequence_symbols_and_embedded_whitespace_fail(self):
        for sequence in ("MW*", "M W", "mwA"):
            with self.subTest(sequence=sequence):
                self.faa.write_text(f">gene_1\n{sequence}\n")
                with self.assertRaises(AnnotationError):
                    self.build()
        self.faa.write_text(">gene_1\nMWA\n")
        self.genome.write_text(">original_contig\nATG-GAGCTTAA\n")
        with self.assertRaisesRegex(AnnotationError, "DNA symbols"):
            self.build()

    def test_code_11_native_translation_and_code_mismatch(self):
        self.faa.write_text(">gene_1\nM\n")
        self.gff.write_text(self.gff.read_text().replace("\t12\t", "\t6\t"))
        self.write_gbk(
            location=SimpleLocation(0, 6, strand=1), translation="M", gcode=11
        )
        self.assertEqual(bundle_proteins(self.build(gcode=11))[0]["genetic_code"], 11)
        with self.assertRaisesRegex(AnnotationError, "translation table"):
            self.build(name="wrong_code")

    def test_genbank_translation_and_coordinates_must_match(self):
        self.write_gbk(translation="MMA")
        with self.assertRaisesRegex(AnnotationError, "translation disagrees"):
            self.build()
        self.write_gbk(location=SimpleLocation(0, 9, strand=1))
        with self.assertRaisesRegex(AnnotationError, "coordinates disagree"):
            self.build()

    def test_compound_origin_segments_follow_genbank_biological_order(self):
        self.gff.write_text(
            "##gff-version 3\n"
            "renamed\tProkka\tCDS\t1\t3\t.\t+\t0\tID=gene_1\n"
            "renamed\tProkka\tCDS\t7\t9\t.\t+\t0\tID=gene_1\n"
            "##FASTA\n>renamed\nATGTGAGCTTAA\n"
        )
        location = CompoundLocation(
            [SimpleLocation(6, 9, strand=1), SimpleLocation(0, 3, strand=1)]
        )
        self.faa.write_text(">gene_1\nAM\n")
        self.write_gbk(location=location, translation="AM", topology="circular")
        output = self.build()
        self.assertEqual(
            [r["start"] for r in read_tsv(output / "gene_coordinates.tsv")], ["7", "1"]
        )
        protein = bundle_proteins(output)[1][0]
        self.assertEqual(
            (protein["segment_count"], protein["topology"]), ("2", "circular")
        )

    def test_partial_location_and_non_cds_features_are_preserved(self):
        self.write_gbk(location=SimpleLocation(BeforePosition(0), 12, strand=1))
        self.gff.write_text(
            self.gff.read_text().replace(
                "##FASTA", "renamed\tProkka\trRNA\t1\t10\t.\t+\t.\tID=rRNA_1\n##FASTA"
            )
        )
        output = self.build()
        self.assertIn("<0", bundle_proteins(output)[1][0]["genbank_location"])
        self.assertIn("rRNA_1", (output / "source.gff").read_text())

    def test_manifest_cannot_omit_source_checksums_or_change_identities(self):
        for field, value in (
            ("files", {}),
            ("input_id", "bad"),
            ("coordinate_id", "bad"),
        ):
            with self.subTest(field=field):
                output = self.build(name=field)
                manifest = json.loads((output / "bundle.json").read_text())
                manifest[field] = value
                (output / "bundle.json").write_text(json.dumps(manifest))
                with self.assertRaises(AnnotationError):
                    bundle_proteins(output)

    def test_cli_retains_explicit_upstream_failure(self):
        log = self.root / "prokka.log"
        log.write_text("exit_code=1\n")
        output = self.root / "failed"
        result = subprocess.run(
            [
                sys.executable,
                str(
                    Path(__file__).resolve().parents[1]
                    / "bin/prepare_annotation_bundle.py"
                ),
                "--accession",
                "SAMPLE_A",
                "--gcode",
                "4",
                "--genome",
                str(self.genome),
                "--faa",
                str(self.faa),
                "--gff",
                str(self.gff),
                "--gbk",
                str(self.gbk),
                "--prokka-log",
                str(log),
                "--outdir",
                str(output),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(
            json.loads((output / "bundle.json").read_text())["status"],
            "upstream_failed",
        )
        with self.assertRaises(AnnotationError):
            bundle_proteins(output)


if __name__ == "__main__":
    unittest.main()
