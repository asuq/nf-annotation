"""Tests for final master-table and sample-status assembly."""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIN_DIR = ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import build_master_table  # noqa: E402
import master_table_contract  # noqa: E402

APPEND_COLUMNS_ASSET = master_table_contract.DEFAULT_APPEND_COLUMNS_ASSET


def read_tsv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV file into a header and row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames is not None
        return reader.fieldnames, list(reader)


class BuildMasterTableTestCase(unittest.TestCase):
    """Cover manifest-driven joins, stable order, and NA propagation."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write a UTF-8 text file and return its path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def test_main_builds_master_table_with_stable_order(self) -> None:
        """Preserve metadata order and append the locked derived columns."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "\n".join(
                    [
                        "accession\tis_new\tassembly_level\tgenome_fasta\torganism_name\tinternal_id",
                        "ACC1\tfalse\tNA\t/path/one.fna\tKnown one\tid_1",
                        "ACC2\ttrue\tScaffold\t/path/two.fna\tCandidate two\tid_2",
                    ]
                )
                + "\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "\n".join(
                    [
                        "Accession\tTax_ID\tOrganism_Name\tAtypical_Warnings",
                        "ACC1\t123\tKnown one\tNA",
                    ]
                )
                + "\n",
            )
            taxonomy = self.write_text_file(
                tmpdir / "taxonomy.tsv",
                "\n".join(
                    [
                        "Tax_ID\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies",
                        "123\tBacteria\tFirmicutes\tMollicutes\tMycoplasmatales\tMycoplasmataceae\tMycoplasma\tsp1",
                    ]
                )
                + "\n",
            )
            checkm2 = self.write_text_file(
                tmpdir / "checkm2.tsv",
                "\n".join(
                    [
                        "accession\tCompleteness_gcode4\tCompleteness_gcode11\tContamination_gcode4\tContamination_gcode11\tCoding_Density_gcode4\tCoding_Density_gcode11\tAverage_Gene_Length_gcode4\tAverage_Gene_Length_gcode11\tTotal_Coding_Sequences_gcode4\tTotal_Coding_Sequences_gcode11\tGcode\tGcode_Rule\tGcode_Length_Ratio\tGcode_Length_Ratio_Threshold\tGcode_Selection_Reason\tLow_quality\tcheckm2_status\twarnings",
                        "ACC1\t95\t80\t2\t1\t0.9\t0.8\t900\t450\t800\t780\t4\tmean_gene_length_ratio\t2\t1.5\tlength_ratio_above_threshold\tfalse\tdone\t",
                        "ACC2\tNA\t72\tNA\t4\tNA\t0.71\tNA\t710\tNA\t610\tNA\tmean_gene_length_ratio\tNA\t1.5\tinvalid_report_pair\tNA\tfailed\tcheckm2_gcode4_failed",
                    ]
                )
                + "\n",
            )
            status_16s = self.write_text_file(
                tmpdir / "16s.tsv",
                "\n".join(
                    [
                        "accession\t16S\tbest_16S_header\tbest_16S_length\twarnings",
                        "ACC1\tYes\th1\t1500\t",
                        "ACC2\tNo\tNA\tNA\t",
                    ]
                )
                + "\n",
            )
            busco_one = self.write_text_file(
                tmpdir / "busco_b.tsv",
                "\n".join(
                    [
                        "accession\tlineage\tBUSCO_bacillota_odb12\tbusco_status\twarnings",
                        "ACC1\tbacillota_odb12\tC:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200\tdone\t",
                        "ACC2\tbacillota_odb12\tNA\tfailed\tbusco_summary_failed",
                    ]
                )
                + "\n",
            )
            busco_two = self.write_text_file(
                tmpdir / "busco_m.tsv",
                "\n".join(
                    [
                        "accession\tlineage\tBUSCO_mycoplasmatota_odb12\tbusco_status\twarnings",
                        "ACC1\tmycoplasmatota_odb12\tC:97.0%[S:97.0%,D:0.0%],F:2.0%,M:1.0%,n:180\tdone\t",
                        "ACC2\tmycoplasmatota_odb12\tNA\tfailed\tbusco_summary_failed",
                    ]
                )
                + "\n",
            )
            ccfinder = self.write_text_file(
                tmpdir / "ccfinder.tsv",
                "\n".join(
                    [
                        "accession\tCRISPRS\tSPACERS_SUM\tCRISPR_FRAC\tccfinder_status\twarnings",
                        "ACC1\t2\t7\t0.1\tdone\t",
                    ]
                )
                + "\n",
            )
            codetta = self.write_text_file(
                tmpdir / "codetta.tsv",
                "\n".join(
                    [
                        "accession\tCodetta_Genetic_Code\tCodetta_NCBI_Table_Candidates\tcodetta_status\twarnings",
                        "ACC1\tFFLLSSSSYY??CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\t11\tdone\t",
                    ]
                )
                + "\n",
            )
            ani = self.write_text_file(
                tmpdir / "ani.tsv",
                "\n".join(
                    [
                        "Accession\tCluster_ID\tIs_Representative\tANI_to_Representative\tScore",
                        "ACC1\tcluster_1\tyes\t100\t0.95",
                    ]
                )
                + "\n",
            )
            master_output = tmpdir / "master_table.tsv"

            exit_code = build_master_table.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--append-columns",
                    str(APPEND_COLUMNS_ASSET),
                    "--taxonomy",
                    str(taxonomy),
                    "--checkm2",
                    str(checkm2),
                    "--16s-status",
                    str(status_16s),
                    "--busco",
                    str(busco_one),
                    "--busco",
                    str(busco_two),
                    "--codetta-summary",
                    str(codetta),
                    "--ccfinder-strains",
                    str(ccfinder),
                    "--ani",
                    str(ani),
                    "--output",
                    str(master_output),
                ]
            )

            self.assertEqual(exit_code, 0)

            master_header, master_rows = read_tsv_rows(master_output)
            expected_master_header = master_table_contract.build_master_table_columns(
                ["Accession", "Tax_ID", "Organism_Name", "Atypical_Warnings"]
            )

            self.assertEqual(master_header, expected_master_header)
            self.assertNotIn("PADLOC", master_header)
            self.assertNotIn("eggNOG", master_header)
            self.assertEqual(len(master_rows), 2)

            master_by_accession = {row["Accession"]: row for row in master_rows}

            self.assertEqual(master_by_accession["ACC1"]["Tax_ID"], "123")
            self.assertEqual(master_by_accession["ACC1"]["is_new"], "false")
            self.assertEqual(master_by_accession["ACC1"]["superkingdom"], "Bacteria")
            self.assertEqual(master_by_accession["ACC1"]["16S"], "Yes")
            self.assertEqual(master_by_accession["ACC1"]["GC_Content"], "NA")
            self.assertEqual(
                master_by_accession["ACC1"]["Codetta_Genetic_Code"],
                "FFLLSSSSYY??CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
            )
            self.assertEqual(master_by_accession["ACC1"]["Codetta_NCBI_Table_Candidates"], "11")
            self.assertEqual(master_by_accession["ACC1"]["CRISPRS"], "2")
            self.assertEqual(master_by_accession["ACC1"]["Cluster_ID"], "cluster_1")

            self.assertEqual(master_by_accession["ACC2"]["Accession"], "ACC2")
            self.assertEqual(master_by_accession["ACC2"]["is_new"], "true")
            self.assertEqual(master_by_accession["ACC2"]["Tax_ID"], "NA")
            self.assertEqual(master_by_accession["ACC2"]["Organism_Name"], "Candidate two")
            self.assertEqual(master_by_accession["ACC2"]["Atypical_Warnings"], "NA")
            self.assertEqual(master_by_accession["ACC2"]["superkingdom"], "NA")
            self.assertEqual(master_by_accession["ACC2"]["GC_Content"], "NA")
            self.assertEqual(master_by_accession["ACC2"]["Codetta_Genetic_Code"], "NA")
            self.assertEqual(master_by_accession["ACC2"]["Codetta_NCBI_Table_Candidates"], "NA")
            self.assertEqual(master_by_accession["ACC2"]["CRISPRS"], "NA")
            self.assertEqual(master_by_accession["ACC2"]["Cluster_ID"], "NA")

    def test_main_populates_custom_busco_lineage_columns(self) -> None:
        """Preserve non-default BUSCO lineage columns from configured lineages."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id\n"
                "ACC1\tfalse\tNA\t/path/one.fna\tid_1\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\tOrganism_Name\n"
                "ACC1\t123\tKnown one\n",
            )
            busco_custom = self.write_text_file(
                tmpdir / "busco_custom.tsv",
                "\n".join(
                    [
                        "accession\tlineage\tBUSCO_custom_odb12\tbusco_status\twarnings",
                        "ACC1\tcustom_odb12\tC:99.0%[S:99.0%,D:0.0%],F:0.0%,M:1.0%,n:180\tdone\t",
                    ]
                )
                + "\n",
            )
            output = tmpdir / "master_table.tsv"

            exit_code = build_master_table.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--busco-lineage",
                    "custom_odb12",
                    "--busco",
                    str(busco_custom),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            header, rows = read_tsv_rows(output)
            self.assertIn("BUSCO_custom_odb12", header)
            self.assertNotIn("BUSCO_bacillota_odb12", header)
            self.assertEqual(rows[0]["BUSCO_custom_odb12"], "C:99.0%[S:99.0%,D:0.0%],F:0.0%,M:1.0%,n:180")

    def test_main_uses_keyed_joins_and_marks_missing_joins(self) -> None:
        """Join by accession only and propagate missing derived rows as NA."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "\n".join(
                    [
                        "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id",
                        "ACC1\tfalse\tNA\t/path/one.fna\tid_1",
                        "ACC2\tfalse\tNA\t/path/two.fna\tid_2",
                    ]
                )
                + "\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "\n".join(
                    [
                        "Accession\tTax_ID\tOrganism_Name\tAtypical_Warnings",
                        "ACC1\t123\tOne\tNA",
                        "ACC2\t456\tTwo\tNA",
                    ]
                )
                + "\n",
            )
            checkm2 = self.write_text_file(
                tmpdir / "checkm2.tsv",
                "\n".join(
                    [
                        "accession\tCompleteness_gcode4\tCompleteness_gcode11\tContamination_gcode4\tContamination_gcode11\tCoding_Density_gcode4\tCoding_Density_gcode11\tAverage_Gene_Length_gcode4\tAverage_Gene_Length_gcode11\tTotal_Coding_Sequences_gcode4\tTotal_Coding_Sequences_gcode11\tGcode\tGcode_Rule\tGcode_Length_Ratio\tGcode_Length_Ratio_Threshold\tGcode_Selection_Reason\tLow_quality\tcheckm2_status\twarnings",
                        "ACC2\t60\t92\t6\t2\t0.6\t0.92\t600\t800\t500\t900\t11\tmean_gene_length_ratio\t0.75\t1.5\tlength_ratio_at_or_below_threshold\tfalse\tdone\t",
                        "ACC1\t94\t80\t1\t2\t0.95\t0.8\t940\t470\t910\t790\t4\tmean_gene_length_ratio\t2\t1.5\tlength_ratio_above_threshold\tfalse\tdone\t",
                    ]
                )
                + "\n",
            )
            status_16s = self.write_text_file(
                tmpdir / "16s.tsv",
                "\n".join(
                    [
                        "accession\t16S\tbest_16S_header\tbest_16S_length\twarnings",
                        "ACC2\tpartial\th2\t1200\t",
                        "ACC1\tYes\th1\t1500\t",
                    ]
                )
                + "\n",
            )
            busco_one = self.write_text_file(
                tmpdir / "busco_b.tsv",
                "\n".join(
                    [
                        "accession\tlineage\tBUSCO_bacillota_odb12\tbusco_status\twarnings",
                        "ACC2\tbacillota_odb12\tC:96.0%[S:96.0%,D:0.0%],F:2.0%,M:2.0%,n:200\tdone\t",
                    ]
                )
                + "\n",
            )
            busco_two = self.write_text_file(
                tmpdir / "busco_m.tsv",
                "\n".join(
                    [
                        "accession\tlineage\tBUSCO_mycoplasmatota_odb12\tbusco_status\twarnings",
                        "ACC1\tmycoplasmatota_odb12\tC:97.0%[S:97.0%,D:0.0%],F:2.0%,M:1.0%,n:180\tdone\t",
                    ]
                )
                + "\n",
            )
            master_output = tmpdir / "master_table.tsv"

            exit_code = build_master_table.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--append-columns",
                    str(APPEND_COLUMNS_ASSET),
                    "--checkm2",
                    str(checkm2),
                    "--16s-status",
                    str(status_16s),
                    "--busco",
                    str(busco_one),
                    "--busco",
                    str(busco_two),
                    "--output",
                    str(master_output),
                ]
            )

            self.assertEqual(exit_code, 0)
            _, master_rows = read_tsv_rows(master_output)
            master_by_accession = {row["Accession"]: row for row in master_rows}

            self.assertEqual(master_by_accession["ACC1"]["Gcode"], "4")
            self.assertEqual(master_by_accession["ACC1"]["Gcode_Rule"], "mean_gene_length_ratio")
            self.assertEqual(master_by_accession["ACC1"]["Gcode_Length_Ratio"], "2")
            self.assertEqual(master_by_accession["ACC1"]["Gcode_Length_Ratio_Threshold"], "1.5")
            self.assertEqual(master_by_accession["ACC1"]["16S"], "Yes")
            self.assertEqual(master_by_accession["ACC1"]["BUSCO_bacillota_odb12"], "NA")
            self.assertEqual(
                master_by_accession["ACC1"]["BUSCO_mycoplasmatota_odb12"],
                "C:97.0%[S:97.0%,D:0.0%],F:2.0%,M:1.0%,n:180",
            )

            self.assertEqual(master_by_accession["ACC2"]["Gcode"], "11")
            self.assertEqual(
                master_by_accession["ACC2"]["Gcode_Selection_Reason"],
                "length_ratio_at_or_below_threshold",
            )
            self.assertEqual(master_by_accession["ACC2"]["16S"], "partial")
            self.assertEqual(
                master_by_accession["ACC2"]["BUSCO_bacillota_odb12"],
                "C:96.0%[S:96.0%,D:0.0%],F:2.0%,M:2.0%,n:200",
            )
            self.assertEqual(master_by_accession["ACC2"]["BUSCO_mycoplasmatota_odb12"], "NA")

    def test_main_consumes_precomputed_ani_summary_during_final_build(self) -> None:
        """Join ANI fields from a precomputed ANI summary TSV."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "\n".join(
                    [
                        "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id",
                        "ACC1\tfalse\tNA\t/path/one.fna\tid_1",
                        "ACC2\tfalse\tNA\t/path/two.fna\tid_2",
                    ]
                )
                + "\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "\n".join(
                    [
                        "Accession\tTax_ID\tOrganism_Name\tAssembly_Level\tN50\tScaffolds\tGenome_Size\tAtypical_Warnings",
                        "ACC1\t123\tOne\tComplete Genome\t100000\t1\t900000\tNA",
                        "ACC2\t456\tTwo\tScaffold\t50000\t5\t850000\tNA",
                    ]
                )
                + "\n",
            )
            checkm2 = self.write_text_file(
                tmpdir / "checkm2.tsv",
                "\n".join(
                    [
                        "accession\tCompleteness_gcode4\tCompleteness_gcode11\tContamination_gcode4\tContamination_gcode11\tCoding_Density_gcode4\tCoding_Density_gcode11\tAverage_Gene_Length_gcode4\tAverage_Gene_Length_gcode11\tTotal_Coding_Sequences_gcode4\tTotal_Coding_Sequences_gcode11\tGcode\tGcode_Rule\tGcode_Length_Ratio\tGcode_Length_Ratio_Threshold\tGcode_Selection_Reason\tLow_quality\tcheckm2_status\twarnings",
                        "ACC1\t90\t97\t2\t1\t0.8\t0.95\t850\t1000\t780\t910\t11\tmean_gene_length_ratio\t0.85\t1.5\tlength_ratio_at_or_below_threshold\tfalse\tdone\t",
                        "ACC2\t88\t93\t3\t2\t0.75\t0.90\t800\t1000\t760\t880\t11\tmean_gene_length_ratio\t0.8\t1.5\tlength_ratio_at_or_below_threshold\tfalse\tdone\t",
                    ]
                )
                + "\n",
            )
            status_16s = self.write_text_file(
                tmpdir / "16s.tsv",
                "\n".join(
                    [
                        "accession\t16S\tbest_16S_header\tbest_16S_length\twarnings",
                        "ACC1\tYes\th1\t1500\t",
                        "ACC2\tYes\th2\t1490\t",
                    ]
                )
                + "\n",
            )
            busco = self.write_text_file(
                tmpdir / "busco.tsv",
                "\n".join(
                    [
                        "accession\tlineage\tBUSCO_bacillota_odb12\tbusco_status\twarnings",
                        "ACC1\tbacillota_odb12\tC:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200\tdone\t",
                        "ACC2\tbacillota_odb12\tC:96.0%[S:96.0%,D:0.0%],F:2.0%,M:2.0%,n:200\tdone\t",
                    ]
                )
                + "\n",
            )
            ani_summary = self.write_text_file(
                tmpdir / "ani_summary.tsv",
                "\n".join(
                    [
                        "Accession\tCluster_ID\tIs_Representative\tANI_to_Representative\tScore",
                        "ACC1\tC000001\tyes\t100.0000\t7.500000",
                        "ACC2\tC000001\tno\t97.2500\t3.750000",
                    ]
                )
                + "\n",
            )
            master_output = tmpdir / "master_table.tsv"

            exit_code = build_master_table.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--append-columns",
                    str(APPEND_COLUMNS_ASSET),
                    "--checkm2",
                    str(checkm2),
                    "--16s-status",
                    str(status_16s),
                    "--busco",
                    str(busco),
                    "--ani",
                    str(ani_summary),
                    "--output",
                    str(master_output),
                ]
            )

            self.assertEqual(exit_code, 0)
            _, master_rows = read_tsv_rows(master_output)
            master_by_accession = {row["Accession"]: row for row in master_rows}

            self.assertEqual(master_by_accession["ACC1"]["Cluster_ID"], "C000001")
            self.assertEqual(master_by_accession["ACC1"]["Is_Representative"], "yes")
            self.assertEqual(master_by_accession["ACC1"]["ANI_to_Representative"], "100.0000")
            self.assertNotEqual(master_by_accession["ACC1"]["Score"], "NA")

            self.assertEqual(master_by_accession["ACC2"]["Cluster_ID"], "C000001")
            self.assertEqual(master_by_accession["ACC2"]["Is_Representative"], "no")
            self.assertEqual(master_by_accession["ACC2"]["ANI_to_Representative"], "97.2500")
            self.assertNotEqual(master_by_accession["ACC2"]["Score"], "NA")

    def test_main_handles_new_genome_with_sparse_metadata(self) -> None:
        """Fill missing metadata with NA while preserving supplemental new-sample values."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "\n".join(
                    [
                        "accession\tis_new\tassembly_level\tgenome_fasta\torganism_name\tinternal_id",
                        "ACC_NEW\ttrue\tContig\t/path/new.fna\tNovel isolate\tnew_id",
                    ]
                )
                + "\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\tOrganism_Name\tAtypical_Warnings\n",
            )
            master_output = tmpdir / "master_table.tsv"

            exit_code = build_master_table.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--append-columns",
                    str(APPEND_COLUMNS_ASSET),
                    "--output",
                    str(master_output),
                ]
            )

            self.assertEqual(exit_code, 0)
            _, master_rows = read_tsv_rows(master_output)

            self.assertEqual(len(master_rows), 1)
            self.assertEqual(master_rows[0]["Accession"], "ACC_NEW")
            self.assertEqual(master_rows[0]["is_new"], "true")
            self.assertEqual(master_rows[0]["Tax_ID"], "NA")
            self.assertEqual(master_rows[0]["Organism_Name"], "Novel isolate")
            self.assertEqual(master_rows[0]["Atypical_Warnings"], "NA")
            self.assertEqual(master_rows[0]["Gcode"], "NA")
            self.assertEqual(master_rows[0]["Codetta_Genetic_Code"], "NA")
            self.assertEqual(master_rows[0]["Codetta_NCBI_Table_Candidates"], "NA")
            self.assertEqual(master_rows[0]["BUSCO_bacillota_odb12"], "NA")
            self.assertEqual(master_rows[0]["CRISPRS"], "NA")

    def test_main_preserves_new_sample_assembly_level_in_metadata_block(self) -> None:
        """Map manifest assembly_level into Assembly_Level when metadata row is absent."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "\n".join(
                    [
                        "accession\tis_new\tassembly_level\tgenome_fasta\torganism_name\tinternal_id",
                        "ACC_NEW\ttrue\tComplete Genome\t/path/new.fna\tNovel isolate\tnew_id",
                    ]
                )
                + "\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\tOrganism_Name\tAssembly_Level\tAtypical_Warnings\n",
            )
            master_output = tmpdir / "master_table.tsv"

            exit_code = build_master_table.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--append-columns",
                    str(APPEND_COLUMNS_ASSET),
                    "--output",
                    str(master_output),
                ]
            )

            self.assertEqual(exit_code, 0)
            _, master_rows = read_tsv_rows(master_output)

            self.assertEqual(master_rows[0]["Accession"], "ACC_NEW")
            self.assertEqual(master_rows[0]["is_new"], "true")
            self.assertEqual(master_rows[0]["Organism_Name"], "Novel isolate")
            self.assertEqual(master_rows[0]["Assembly_Level"], "Complete Genome")
            self.assertEqual(master_rows[0]["Atypical_Warnings"], "NA")

    def test_main_rejects_metadata_header_overlap_with_appended_is_new(self) -> None:
        """Fail when metadata already declares the appended is_new column."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "\n".join(
                    [
                        "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id",
                        "ACC1\tfalse\tNA\t/path/one.fna\tid_1",
                    ]
                )
                + "\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "\n".join(
                    [
                        "Accession\tis_new\tTax_ID",
                        "ACC1\tfalse\t123",
                    ]
                )
                + "\n",
            )
            master_output = tmpdir / "master_table.tsv"

            exit_code = build_master_table.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--append-columns",
                    str(APPEND_COLUMNS_ASSET),
                    "--output",
                    str(master_output),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse(master_output.exists())

    def test_main_rejects_duplicate_accessions_in_validated_manifest(self) -> None:
        """Fail when validated_samples.tsv contains duplicate accession rows."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "\n".join(
                    [
                        "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id",
                        "ACC1\tfalse\tNA\t/path/one.fna\tid_1",
                        "ACC1\ttrue\tScaffold\t/path/two.fna\tid_2",
                    ]
                )
                + "\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\nACC1\t123\n",
            )
            master_output = tmpdir / "master_table.tsv"

            exit_code = build_master_table.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--append-columns",
                    str(APPEND_COLUMNS_ASSET),
                    "--output",
                    str(master_output),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse(master_output.exists())

    def test_main_rejects_duplicate_rows_in_derived_tables(self) -> None:
        """Fail when a derived accession-keyed table contains ambiguous duplicates."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "\n".join(
                    [
                        "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id",
                        "ACC1\tfalse\tNA\t/path/one.fna\tid_1",
                    ]
                )
                + "\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\nACC1\t123\n",
            )
            duplicate_checkm2 = self.write_text_file(
                tmpdir / "checkm2.tsv",
                "\n".join(
                    [
                        "accession\tCompleteness_gcode4\tCompleteness_gcode11\tContamination_gcode4\tContamination_gcode11\tCoding_Density_gcode4\tCoding_Density_gcode11\tAverage_Gene_Length_gcode4\tAverage_Gene_Length_gcode11\tTotal_Coding_Sequences_gcode4\tTotal_Coding_Sequences_gcode11\tGcode\tGcode_Rule\tGcode_Length_Ratio\tGcode_Length_Ratio_Threshold\tGcode_Selection_Reason\tLow_quality\tcheckm2_status\twarnings",
                        "ACC1\t95\t80\t2\t1\t0.9\t0.8\t900\t450\t800\t780\t4\tmean_gene_length_ratio\t2\t1.5\tlength_ratio_above_threshold\tfalse\tdone\t",
                        "ACC1\t95\t80\t2\t1\t0.9\t0.8\t900\t450\t800\t780\t4\tmean_gene_length_ratio\t2\t1.5\tlength_ratio_above_threshold\tfalse\tdone\t",
                    ]
                )
                + "\n",
            )
            master_output = tmpdir / "master_table.tsv"

            exit_code = build_master_table.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--append-columns",
                    str(APPEND_COLUMNS_ASSET),
                    "--checkm2",
                    str(duplicate_checkm2),
                    "--output",
                    str(master_output),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse(master_output.exists())

    def test_main_rejects_append_contract_drift(self) -> None:
        """Fail when the append-column asset drifts from the locked contract."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "\n".join(
                    [
                        "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id",
                        "ACC1\tfalse\tNA\t/path/one.fna\tid_1",
                    ]
                )
                + "\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\nACC1\t123\n",
            )
            append_columns = self.write_text_file(
                tmpdir / "append_columns.txt",
                "\n".join(
                    master_table_contract.build_append_columns()[:-1]
                    + ["Unexpected_Column"]
                )
                + "\n",
            )
            master_output = tmpdir / "master_table.tsv"

            exit_code = build_master_table.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--append-columns",
                    str(append_columns),
                    "--output",
                    str(master_output),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse(master_output.exists())

    def test_main_backfills_missing_metadata_stats_from_in_house_values(self) -> None:
        """Backfill missing metadata metrics from computed assembly stats only."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "\n".join(
                    [
                        "accession\tis_new\tassembly_level\tgenome_fasta\torganism_name\tinternal_id",
                        "ACC1\tfalse\tNA\t/path/one.fna\tKnown one\tid_1",
                        "ACC2\ttrue\tComplete Genome\t/path/two.fna\tCandidate two\tid_2",
                    ]
                )
                + "\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "\n".join(
                    [
                        "Accession\tTax_ID\tOrganism_Name\tN50\tScaffolds\tGenome_Size\tAtypical_Warnings",
                        "ACC1\t123\tKnown one\t50000\tNA\tNA\tNA",
                    ]
                )
                + "\n",
            )
            assembly_stats = self.write_text_file(
                tmpdir / "assembly_stats.tsv",
                "\n".join(
                    [
                        "accession\tn50\tscaffolds\tgenome_size\tgc_content",
                        "ACC1\t90000\t2\t800000\t37.5",
                        "ACC2\t120000\t1\t120000\t50",
                    ]
                )
                + "\n",
            )
            master_output = tmpdir / "master_table.tsv"

            exit_code = build_master_table.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--assembly-stats",
                    str(assembly_stats),
                    "--append-columns",
                    str(APPEND_COLUMNS_ASSET),
                    "--output",
                    str(master_output),
                ]
            )

            self.assertEqual(exit_code, 0)
            _, master_rows = read_tsv_rows(master_output)
            master_by_accession = {row["Accession"]: row for row in master_rows}

            self.assertEqual(master_by_accession["ACC1"]["N50"], "50000")
            self.assertEqual(master_by_accession["ACC1"]["Scaffolds"], "2")
            self.assertEqual(master_by_accession["ACC1"]["Genome_Size"], "800000")
            self.assertEqual(master_by_accession["ACC1"]["GC_Content"], "37.5")
            self.assertEqual(master_by_accession["ACC2"]["N50"], "120000")
            self.assertEqual(master_by_accession["ACC2"]["Scaffolds"], "1")
            self.assertEqual(master_by_accession["ACC2"]["Genome_Size"], "120000")
            self.assertEqual(master_by_accession["ACC2"]["GC_Content"], "50")

    def test_main_does_not_insert_missing_metadata_metric_columns(self) -> None:
        """Avoid adding new metadata columns when the input header omits them."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "accession\tis_new\tassembly_level\tgenome_fasta\torganism_name\tinternal_id\nACC1\ttrue\tComplete Genome\t/path/one.fna\tCandidate one\tid_1\n",
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\tOrganism_Name\tAtypical_Warnings\n",
            )
            assembly_stats = self.write_text_file(
                tmpdir / "assembly_stats.tsv",
                "accession\tn50\tscaffolds\tgenome_size\tgc_content\nACC1\t120000\t1\t120000\t50\n",
            )
            master_output = tmpdir / "master_table.tsv"

            exit_code = build_master_table.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--metadata",
                    str(metadata),
                    "--assembly-stats",
                    str(assembly_stats),
                    "--append-columns",
                    str(APPEND_COLUMNS_ASSET),
                    "--output",
                    str(master_output),
                ]
            )

            self.assertEqual(exit_code, 0)
            master_header, master_rows = read_tsv_rows(master_output)
            self.assertNotIn("N50", master_header)
            self.assertNotIn("Scaffolds", master_header)
            self.assertNotIn("Genome_Size", master_header)
            self.assertEqual(master_rows[0]["Accession"], "ACC1")


if __name__ == "__main__":
    unittest.main()
