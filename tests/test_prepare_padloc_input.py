"""Protect PADLOC's native FAA/GFF join without changing gene geometry."""

import unittest

import test_annotation_bundle as bundle_fixture
from annotation_common import AnnotationError, bundle_proteins
from prepare_padloc_input import prepare_gff


class PadlocInputTests(unittest.TestCase):
    def setUp(self):
        self.fixture = bundle_fixture.AnnotationBundleTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)

    def test_prokka_protein_attribute_matches_faa_without_changing_other_fields(self):
        original = self.fixture.gff.read_text().replace(
            "locus_tag=gene_1", "locus_tag=gene_1;protein_id=gnl|Prokka|gene_1"
        )
        self.fixture.gff.write_text(original)
        bundle = self.fixture.build()
        _, proteins = bundle_proteins(bundle)
        output = self.fixture.root / "padloc.gff"
        prepare_gff(bundle, proteins, output)
        self.assertEqual(
            output.read_text(),
            original.split("##FASTA")[0].replace(
                "protein_id=gnl|Prokka|gene_1", "protein_id=gene_1"
            ),
        )
        self.assertEqual((bundle / "source.gff").read_text(), original)

    def test_compound_segments_keep_original_gff_order_and_positions(self):
        self.fixture.test_compound_origin_segments_follow_genbank_biological_order()
        bundle = self.fixture.root / "bundle"
        _, proteins = bundle_proteins(bundle)
        output = self.fixture.root / "padloc.gff"
        prepare_gff(bundle, proteins, output)
        rows = [
            line.split("\t")
            for line in output.read_text().splitlines()
            if not line.startswith("#")
        ]
        self.assertEqual([row[3:5] for row in rows], [["1", "3"], ["7", "9"]])
        self.assertEqual([row[8] for row in rows], ["ID=gene_1;protein_id=gene_1"] * 2)

    def test_unknown_or_missing_source_segments_are_rejected(self):
        bundle = self.fixture.build()
        _, proteins = bundle_proteins(bundle)
        source = bundle / "source.gff"
        original = source.read_text()
        for altered in (
            original.replace("\t12\t", "\t11\t"),
            "\n".join(line for line in original.splitlines() if "\tCDS\t" not in line),
        ):
            with self.subTest(source=altered):
                source.write_text(altered)
                with self.assertRaises(AnnotationError):
                    prepare_gff(bundle, proteins, self.fixture.root / "padloc.gff")
