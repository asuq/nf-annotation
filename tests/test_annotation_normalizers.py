"""Scientific regression controls for native v0.4 annotation normalization."""

from __future__ import annotations

import csv
import json
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import normalize_cogclassifier as cog
import normalize_eggnog as egg
import normalize_kofam as kofam
import normalize_padloc as padloc
import normalize_pfam as pfam
from annotation_common import AnnotationError, write_tsv


class NormalizerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.raw, self.resource, self.bundle = [
            self.root / name for name in ("raw", "resource", "bundle")
        ]
        for path in (self.raw, self.resource, self.bundle):
            path.mkdir()
        self.proteins = [
            dict(
                accession="genome.1",
                gene_id=f"genome.1::gene{i}",
                protein_id=f"gene{i}",
                tool_id=f"p{i}",
                length="200",
            )
            for i in range(1, 4)
        ]

    def tearDown(self):
        self.temporary.cleanup()

    def write_eggnog(self, rows=None, go_rows=None):
        (self.resource / "go-basic.obo").write_text(
            "format-version: 1.2\n\n"
            + "\n".join(
                f"[Term]\nid: GO:{index:07d}\nname: Control {index}\nnamespace: {namespace}\n"
                for index, namespace in enumerate(
                    ("molecular_function", "biological_process", "cellular_component"),
                    1,
                )
            )
        )
        if rows is None:
            row = dict.fromkeys(egg.HEADER, "-")
            row.update(
                query="p1",
                seed_ortholog="1234",
                evalue="0",
                score="100",
                GOs="GO:0000001,GO:0000002,GO:0000003",
                EC="ec:1.2.3.4",
                KEGG_ko="K00001",
                COG_category="RS",
                annotation_confidence="-hml" + "-" * 9,
            )
            rows = [row]
        if go_rows is None:
            go_rows = [
                dict(
                    query="p1",
                    gos_mf="GO:0000001",
                    gos_mf_confidence="high",
                    gos_bp="GO:0000002",
                    gos_bp_confidence="medium",
                    gos_cc="GO:0000003",
                    gos_cc_confidence="low",
                )
            ]
        with (self.raw / "eggnog.emapper.annotations").open("w", newline="") as handle:
            handle.write("## confidence codes: h=high m=medium l=low -=not annotated\n")
            handle.write("## confidence field order: " + " ".join(egg.FIELDS) + "\n")
            handle.write("#" + "\t".join(egg.HEADER) + "\n")
            writer = csv.DictWriter(
                handle, fieldnames=egg.HEADER, delimiter="\t", lineterminator="\n"
            )
            writer.writerows(rows)
            handle.write(f"## {len(rows)} queries scanned\n")
        with (self.raw / "eggnog.emapper.seed_orthologs").open("w") as handle:
            handle.write(
                "#qseqid\tsseqid\tevalue\tbitscore\tqstart\tqend\tsstart\tsend\tpident\tqcov\tscov\n"
            )
            for row in rows:
                handle.write(
                    f"{row['query']}\t1\t{row['evalue']}\t{row['score']}\t1\t100\t1\t100\t90\t50\t50\n"
                )
            handle.write(f"## {len(rows)} queries scanned\n")
        with sqlite3.connect(self.resource / "eggnog.db") as database:
            database.execute(
                "CREATE TABLE IF NOT EXISTS protein_names (id INTEGER PRIMARY KEY, name TEXT)"
            )
            database.execute("INSERT OR REPLACE INTO protein_names VALUES (1, '1234')")
        write_tsv(
            self.raw / "eggnog.emapper.annotations.go_namespaces.tsv",
            egg.GO_HEADER,
            go_rows,
        )

    def test_eggnog_filters_each_namespace_and_field(self):
        self.write_eggnog()
        result = egg.normalize(self.raw, self.proteins, self.resource)
        gene = self.proteins[0]["gene_id"]
        self.assertEqual(result.features["go"][gene], {"GO:0000001", "GO:0000002"})
        self.assertEqual(result.features["ec"][gene], {"1.2.3.4"})
        self.assertNotIn("ko", result.features)
        self.assertEqual(result.errors, [])
        self.assertIn(
            "GO:0000003", result.tables["eggnog_annotations.tsv"][1][0]["go_namespaces"]
        )

    def test_eggnog_optional_field_error_keeps_independent_assignments(self):
        self.write_eggnog()
        path = self.raw / "eggnog.emapper.annotations"
        path.write_text(
            path.read_text().replace("1.2.3.4", "bad_ec").replace("\tRS\t", "\tO \t")
        )
        result = egg.normalize(self.raw, self.proteins, self.resource)
        self.assertEqual({row["field"] for row in result.errors}, {"ec", "categories"})
        self.assertIn("go", result.features)
        self.assertNotIn("ec", result.features)
        self.assertIn("O ", [row["raw_value"] for row in result.errors])

    def test_eggnog_native_ec_namespace_and_cog_identifier_quarantine(self):
        self.write_eggnog()
        path = self.raw / "eggnog.emapper.annotations"
        path.write_text(
            path.read_text()
            .replace("ec:1.2.3.4", "ec:5.6.2.2,ec:2.7.4.9")
            .replace("\tRS\t", "\tCOG0484\t")
        )
        result = egg.normalize(self.raw, self.proteins, self.resource)
        gene = self.proteins[0]["gene_id"]
        self.assertEqual(result.features["ec"][gene], {"5.6.2.2", "2.7.4.9"})
        self.assertEqual(
            [(row["field"], row["raw_value"]) for row in result.errors],
            [("categories", "COG0484")],
        )
        evidence = result.tables["eggnog_annotations.tsv"][1][0]
        self.assertEqual(evidence["EC"], "ec:5.6.2.2,ec:2.7.4.9")
        self.assertEqual(
            json.loads(evidence["accepted_fields"])["EC"], ["5.6.2.2", "2.7.4.9"]
        )
        self.assertEqual(evidence["categories"], "null")

    def test_eggnog_ec_namespace_must_be_explicit_and_well_formed(self):
        for value in (
            "5.6.2.2",
            "EC:5.6.2.2",
            "ec:",
            "ec:ec:5.6.2.2",
            "ec:5.6.2.2,2.7.4.9",
        ):
            with self.subTest(value=value), self.assertRaises(AnnotationError):
                egg.values(value, "EC")

    def test_eggnog_namespace_must_match_prepared_ontology(self):
        self.write_eggnog()
        path = self.resource / "go-basic.obo"
        path.write_text(
            path.read_text().replace(
                "namespace: molecular_function", "namespace: cellular_component"
            )
        )
        result = egg.normalize(self.raw, self.proteins, self.resource)
        self.assertEqual({error["field"] for error in result.errors}, {"go"})
        self.assertNotIn("go", result.features)
        self.assertIn("ec", result.features)

    def test_eggnog_native_ko_identifiers_and_mapped_without_assignment(self):
        self.write_eggnog()
        path = self.raw / "eggnog.emapper.annotations"
        path.write_text(path.read_text().replace("-hml", "-hmm"))
        result = egg.normalize(self.raw, self.proteins, self.resource)
        self.assertEqual(result.features["ko"][self.proteins[0]["gene_id"]], {"K00001"})
        seeds = (self.raw / "eggnog.emapper.seed_orthologs").read_text()
        self.write_eggnog(rows=[], go_rows=[])
        (self.raw / "eggnog.emapper.seed_orthologs").write_text(seeds)
        result = egg.normalize(self.raw, self.proteins, self.resource)
        self.assertEqual(result.mapped, {self.proteins[0]["gene_id"]})
        self.assertEqual(result.accepted, set())
        self.assertEqual(result.reported_hits, 1)

    def test_eggnog_rejects_truncated_seed_evidence(self):
        self.write_eggnog()
        path = self.raw / "eggnog.emapper.seed_orthologs"
        path.write_text(path.read_text().replace("## 1 queries scanned\n", ""))
        with self.assertRaisesRegex(AnnotationError, "completion count"):
            egg.normalize(self.raw, self.proteins, self.resource)

    def test_eggnog_rejects_missing_legend_schema_ids_and_sidecar(self):
        for old, new in (
            ("## confidence codes:", "## removed:"),
            ("#query", "#unknown"),
            ("\np1\t", "\nunknown\t"),
        ):
            with self.subTest(old=old):
                self.write_eggnog()
                path = self.raw / "eggnog.emapper.annotations"
                path.write_text(path.read_text().replace(old, new))
                with self.assertRaises(AnnotationError):
                    egg.normalize(self.raw, self.proteins, self.resource)
        self.write_eggnog(go_rows=[])
        with self.assertRaises(AnnotationError):
            egg.normalize(self.raw, self.proteins, self.resource)

    def test_eggnog_valid_zero_and_malformed_confidence(self):
        self.write_eggnog(rows=[], go_rows=[])
        result = egg.normalize(self.raw, self.proteins, self.resource)
        self.assertEqual((result.reported_hits, result.features), (0, {}))
        self.write_eggnog()
        path = self.raw / "eggnog.emapper.annotations"
        path.write_text(path.read_text().replace("-hml" + "-" * 9, "?"))
        result = egg.normalize(self.raw, self.proteins, self.resource)
        self.assertEqual(
            len({error["field"] for error in result.errors}), len(egg.FIELDS)
        )
        self.assertFalse(result.features)

    def write_cog(self, category="RS"):
        (self.resource / "cog_func_category.tsv").write_text(
            "R\tG\t#fff\tGeneral\nS\tG\t#fff\tUnknown\nO\tG\t#fff\tProcessing\n"
        )
        (self.resource / "cog_definition.tsv").write_text(
            f"COG0001\t{category}\tDefinition one\tx\nCOG0002\tS\tDefinition two\ty\n"
        )
        (self.resource / "cddid.tbl").write_text("1\tCOG0001\n2\tCOG0002\n")
        (self.raw / "rpsblast.tsv").write_text(
            "p1\tCDD:1\t90\t100\t1\t0\t1\t100\t1\t100\t0\t200\n"
            "p1\tCDD:2\t90\t100\t1\t0\t1\t100\t1\t100\t0\t200\n"
        )
        write_tsv(
            self.raw / "cogclassifier.native.tsv",
            cog.NATIVE_COLUMNS,
            [
                dict(
                    QUERY_ID="p1",
                    COG_ID="COG0001",
                    CDD_ID="1",
                    COG_LETTER=category[0],
                    EVALUE="0",
                    IDENTITY="90",
                    GENE_NAME="x",
                    COG_NAME="Definition one",
                    COG_DESCRIPTION="native description",
                )
            ],
        )

    def test_cog_native_first_hit_zero_evalue_and_multicategory(self):
        self.write_cog()
        result = cog.normalize(self.raw, self.proteins, self.resource)
        gene = self.proteins[0]["gene_id"]
        self.assertEqual(result.features["cog"][gene], {"COG0001"})
        self.assertEqual(result.features["categories"][gene], {"R", "S"})
        self.assertEqual(result.reported_hits, 2)
        path = self.raw / "cogclassifier.native.tsv"
        path.write_text(path.read_text().replace("COG0001\t1", "COG0002\t2"))
        with self.assertRaisesRegex(AnnotationError, "first reported hit"):
            cog.normalize(self.raw, self.proteins, self.resource)

    def test_real_cog6144_optional_category_error_does_not_drop_cog(self):
        self.write_cog(category="O ")
        result = cog.normalize(self.raw, self.proteins, self.resource)
        self.assertIn("cog", result.features)
        self.assertNotIn("categories", result.features)
        self.assertEqual(result.errors[0]["raw_value"], "O ")
        self.assertEqual(
            result.tables["cog_assignments.tsv"][1][0]["native_primary_category"], "O"
        )

    def write_kofam(self):
        (self.resource / "ko_list").write_text(
            "knum\tthreshold\tscore_type\tdefinition\n"
            "K00001\t100\tfull\tone\nK00002\t100\tdomain\ttwo\nK00003\t-\t-\tthree\n"
        )
        (self.resource / "prokaryote_profiles.tsv").write_text(
            "ko\tprofile\tthreshold\tscore_type\n"
            "K00001\tK00001.hmm\t100\tfull\nK00002\tK00002.hmm\t100\tdomain\nK00003\tK00003.hmm\t-\t-\n"
        )
        (self.raw / "kofam.tsv").write_text(
            '#\tgene name\tKO\tthrshld\tscore\tE-value\t"KO definition"\n'
            "#\t---------\t------\t-------\t------\t---------\t-------------\n"
            '*\tp1\tK00001\t100.00\t100.0\t0\t"one"\n'
            '\tp1\tK00002\t100.00\t100.0\t1e-10\t"two"\n'
            '\tp2\tK00003\t\t100.0\t1e-10\t"three"\n'
        )

    def test_kofam_marker_rounding_and_unavailable_threshold(self):
        self.write_kofam()
        result = kofam.normalize(self.raw, self.proteins, self.resource)
        self.assertEqual(result.features["ko"][self.proteins[0]["gene_id"]], {"K00001"})
        self.assertEqual(
            [row["state"] for row in result.tables["kofam_hits.tsv"][1]],
            ["threshold_pass", "threshold_fail", "threshold_unavailable"],
        )
        self.assertEqual(result.reported_hits, 3)

    def test_kofam_rejects_unselected_profile_and_unavailable_pass(self):
        for old, new in (
            ("K00001\t100.00", "K99999\t100.00"),
            ("\tp2\t", "\tunknown\t"),
            ("\tp2\tK00003", "*\tp2\tK00003"),
        ):
            with self.subTest(old=old):
                self.write_kofam()
                path = self.raw / "kofam.tsv"
                path.write_text(path.read_text().replace(old, new))
                with self.assertRaises(AnnotationError):
                    kofam.normalize(self.raw, self.proteins, self.resource)

    def write_pfam(self):
        (self.resource / "Pfam-A.hmm.dat").write_text(
            "#=GF ID Model\n#=GF AC PF00001.1\n#=GF DE Definition\n#=GF CL CL0001\n//\n"
        )
        scans = [
            f"p1 {start} {end} {start} {end} PF00001.1 Model Domain 1 50 50 100.0 1e-20 1 CL0001"
            for start, end in ((1, 50), (101, 150))
        ]
        for mode in ("raw", "resolved"):
            (self.raw / f"pfam.{mode}.tsv").write_text(
                "# <seq id> <alignment start>\n#    resolve clan overlaps: "
                + ("off" if mode == "raw" else "on")
                + "\n"
                + "\n".join(scans)
                + "\n"
            )
        rows = [
            f"Model PF00001.1 50 p1 - 200 1e-30 200 0 {i} 2 1e-22 1e-20 100 0 1 50 {start} {end} {start} {end} 0.99 Definition"
            for i, (start, end) in enumerate(((1, 50), (101, 150)), 1)
        ]
        (self.raw / "hmmscan.domtblout").write_text(
            "# Program: hmmscan\n" + "\n".join(rows) + "\n# [ok]\n"
        )

    def test_pfam_repeated_domains_count_once_per_gene(self):
        self.write_pfam()
        result = pfam.normalize(self.raw, self.proteins, self.resource)
        self.assertEqual(result.reported_hits, 2)
        self.assertEqual(
            result.features["pfam"], {self.proteins[0]["gene_id"]: {"PF00001"}}
        )
        path = self.raw / "pfam.resolved.tsv"
        path.write_text("\n".join(path.read_text().splitlines()[:-1]) + "\n")
        result = pfam.normalize(self.raw, self.proteins, self.resource)
        self.assertEqual(
            [row["state"] for row in result.tables["pfam_domains.tsv"][1]],
            ["accepted", "clan_overlap"],
        )

    def test_pfam_truncation_and_impossible_coordinates_fail(self):
        for old, new in (("# [ok]", ""), ("101 150 101 150", "101 250 101 250")):
            with self.subTest(old=old):
                self.write_pfam()
                path = self.raw / "hmmscan.domtblout"
                path.write_text(path.read_text().replace(old, new))
                with self.assertRaises(AnnotationError):
                    pfam.normalize(self.raw, self.proteins, self.resource)

    def write_padloc(self, *, positive=True):
        coords = [
            dict(
                gene_id=protein["gene_id"],
                contig_id="contig",
                start=str(i * 600 + 1),
                end=str(i * 600 + 600),
                strand="+",
                segment="1",
            )
            for i, protein in enumerate(self.proteins)
        ]
        write_tsv(self.bundle / "gene_coordinates.tsv", tuple(coords[0]), coords)
        domains = (
            [
                f"gene{i + 1} - 200 model PDLC0001 200 0 200 0 1 1 0 1e-10 100 0 1 150 1 150 1 150 0.99 Description"
                for i in range(2)
            ]
            if positive
            else []
        )
        (self.raw / "input.domtblout").write_text(
            "# Program: hmmsearch\n" + "\n".join(domains) + "\n# [ok]\n"
        )
        if not positive:
            (self.raw / "tool.log").write_text("[time] Nothing found for input\n")
            return
        rows = []
        for i in range(2):
            row = dict.fromkeys(padloc.HEADER, "evidence")
            row.update(
                {
                    "system.number": "1",
                    "seqid": "contig",
                    "system": "defence_type",
                    "target.name": f"gene{i + 1}",
                    "hmm.accession": "PDLC0001",
                    "hmm.name": "model",
                    "protein.name": "defence_gene",
                    "full.seq.E.value": "0",
                    "domain.iE.value": "1e-10",
                    "target.coverage": "0.9",
                    "hmm.coverage": "0.8",
                    "start": coords[i]["start"],
                    "end": coords[i]["end"],
                    "strand": "+",
                    "relative.position": str(i + 1),
                    "contig.end": "3",
                }
            )
            rows.append(row)
        with (self.raw / "input_padloc.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=padloc.HEADER)
            writer.writeheader()
            writer.writerows(rows)

    def test_padloc_system_not_member_count_and_coordinate_identity(self):
        self.write_padloc()
        result = padloc.normalize(self.raw, self.proteins, self.bundle)
        self.assertEqual(len(result.tables["defence_systems.tsv"][1]), 1)
        self.assertEqual(len(result.tables["defence_genes.tsv"][1]), 2)
        self.assertEqual(result.tables["defence_systems.tsv"][1][0]["member_genes"], 2)
        path = self.raw / "input_padloc.csv"
        path.write_text(path.read_text().replace(",601,1200,", ",602,1200,"))
        with self.assertRaisesRegex(AnnotationError, "coordinate segment"):
            padloc.normalize(self.raw, self.proteins, self.bundle)

    def test_padloc_native_zero_requires_complete_evidence_and_declaration(self):
        self.write_padloc(positive=False)
        result = padloc.normalize(self.raw, self.proteins, self.bundle)
        self.assertEqual(result.reported_hits, 0)
        self.assertEqual(result.tables["defence_systems.tsv"][1], [])
        (self.raw / "tool.log").write_text("incomplete run\n")
        with self.assertRaises(AnnotationError):
            padloc.normalize(self.raw, self.proteins, self.bundle)


if __name__ == "__main__":
    unittest.main()
