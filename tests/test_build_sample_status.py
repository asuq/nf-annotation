"""Tests for the final authoritative sample-status writer."""

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

import build_sample_status  # noqa: E402

SAMPLE_STATUS_COLUMNS_ASSET = (
    ROOT / "assets" / "tables" / "contracts" / "sample_status_columns.txt"
)


def read_tsv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    """Read a TSV file into a header and row dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames is not None
        return reader.fieldnames, list(reader)


class BuildSampleStatusTestCase(unittest.TestCase):
    """Cover the authoritative status-writer contract and failure modes."""

    def write_text_file(self, path: Path, content: str) -> Path:
        """Write a UTF-8 text file and return its path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        return path

    def write_tsv_rows(
        self,
        path: Path,
        header: list[str],
        rows: list[dict[str, str]],
    ) -> Path:
        """Write TSV rows with a fixed header and return the path."""
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
        return path

    def make_initial_status_row(
        self,
        *,
        accession: str,
        internal_id: str,
        is_new: str,
        columns: list[str] | None = None,
        warnings: str = "",
        notes: str = "",
        **overrides: str,
    ) -> dict[str, str]:
        """Build one seed status row matching the locked sample-status contract."""
        if columns is None:
            columns = [
                line.strip()
                for line in SAMPLE_STATUS_COLUMNS_ASSET.read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]
        row: dict[str, str] = {}
        for column in columns:
            if column.endswith("_status") or column == "ani_included":
                row[column] = "na"
            elif column in {"gcode", "low_quality"}:
                row[column] = "NA"
            else:
                row[column] = ""
        row["accession"] = accession
        row["internal_id"] = internal_id
        row["is_new"] = is_new
        row["validation_status"] = "done"
        row["warnings"] = warnings
        row["notes"] = notes
        row.update(overrides)
        return row

    def make_ani_metrics(self) -> build_sample_status.EffectiveAssemblyMetrics:
        """Build complete assembly metrics for direct ANI-decision tests."""
        return build_sample_status.EffectiveAssemblyMetrics(
            n50="50000",
            scaffolds="2",
            genome_size="800000",
        )

    def test_ani_decision_defaults_exclude_no_and_partial_16s(self) -> None:
        """Exclude incomplete 16S samples from ANI by default."""
        for sixteen_s_value, expected_reason in (
            ("No", "no_16s"),
            ("partial", "partial_16s"),
        ):
            with self.subTest(sixteen_s_value=sixteen_s_value):
                included, reason = build_sample_status.derive_ani_decision(
                    "ACC1",
                    gcode_value="11",
                    low_quality_value="false",
                    sixteen_s_value=sixteen_s_value,
                    assembly_level_value="Scaffold",
                    metadata_row={"Atypical_Warnings": "NA"},
                    assembly_metrics=self.make_ani_metrics(),
                    ani_index={},
                    ani_requested=True,
                    primary_busco_value="C:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200",
                )

                self.assertEqual(included, "false")
                self.assertEqual(reason, expected_reason)

    def test_ani_decision_three_16s_policies_and_independent_quality_gates(
        self,
    ) -> None:
        """Report the same inclusion policy as FastANI input selection."""
        for policy, allowed in (
            ("complete", {"Yes"}),
            ("allow_incomplete", {"Yes", "partial"}),
            ("ignore", {"Yes", "partial", "No", "NA"}),
        ):
            for status in ("Yes", "partial", "No", "NA"):
                for low_quality, atypical in (
                    (False, False),
                    (True, False),
                    (False, True),
                ):
                    with self.subTest(
                        policy=policy,
                        status=status,
                        low_quality=low_quality,
                        atypical=atypical,
                    ):
                        eligible = (
                            status in allowed and not low_quality and not atypical
                        )
                        included, reason = build_sample_status.derive_ani_decision(
                            "ACC1",
                            gcode_value="11",
                            low_quality_value=str(low_quality).lower(),
                            sixteen_s_value=status,
                            assembly_level_value="Scaffold",
                            metadata_row={
                                "Atypical_Warnings": "contaminated"
                                if atypical
                                else "NA"
                            },
                            assembly_metrics=self.make_ani_metrics(),
                            ani_index={"ACC1": {}} if eligible else {},
                            ani_requested=True,
                            primary_busco_value="C:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200",
                            ani_16s_policy=policy,
                        )
                        self.assertEqual(included, str(eligible).lower())
                        if low_quality:
                            self.assertIn("low_quality", reason)
                        if atypical:
                            self.assertIn("atypical", reason)
                        if policy == "ignore":
                            self.assertNotIn("16s", reason)

    def test_main_builds_authoritative_status_with_stable_order(self) -> None:
        """Preserve seed warnings while overlaying downstream status columns."""
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
            initial_status_rows = [
                self.make_initial_status_row(
                    accession="ACC1",
                    internal_id="id_1",
                    is_new="false",
                    warnings="internal_id_collision_resolved",
                ),
                self.make_initial_status_row(
                    accession="ACC2",
                    internal_id="id_2",
                    is_new="true",
                    warnings="missing_metadata_for_new_sample",
                    notes="New sample metadata missing; unavailable fields filled with NA.",
                ),
            ]
            initial_status = self.write_tsv_rows(
                tmpdir / "initial_status.tsv",
                [row for row in initial_status_rows[0]],
                initial_status_rows,
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "\n".join(
                    [
                        "Accession\tTax_ID\tOrganism_Name\tAssembly_Level\tN50\tScaffolds\tGenome_Size\tAtypical_Warnings",
                        "ACC1\t123\tKnown one\tScaffold\tNA\tNA\tNA\tNA",
                    ]
                )
                + "\n",
            )
            assembly_stats = self.write_text_file(
                tmpdir / "assembly_stats.tsv",
                "\n".join(
                    [
                        "accession\tn50\tscaffolds\tgenome_size\tgc_content",
                        "ACC1\t100000\t1\t100000\t50",
                        "ACC2\t80000\t2\t900000\t48",
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
            codetta_summary = self.write_text_file(
                tmpdir / "codetta.tsv",
                "\n".join(
                    [
                        "accession\tCodetta_Genetic_Code\tCodetta_NCBI_Table_Candidates\tcodetta_status\twarnings",
                        "ACC1\tFFLLSSSSYY??CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG\t11\tdone\t",
                        "ACC2\tNA\tNA\tfailed\tcodetta_failed",
                    ]
                )
                + "\n",
            )
            prokka_manifest = self.write_text_file(
                tmpdir / "prokka_manifest.tsv",
                "accession\texit_code\tgff_size\tfaa_size\nACC1\t0\t100\t50\n",
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
            output = tmpdir / "sample_status.tsv"

            exit_code = build_sample_status.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--initial-status",
                    str(initial_status),
                    "--metadata",
                    str(metadata),
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
                    str(codetta_summary),
                    "--ccfinder-strains",
                    str(ccfinder),
                    "--prokka-manifest",
                    str(prokka_manifest),
                    "--ani",
                    str(ani),
                    "--assembly-stats",
                    str(assembly_stats),
                    "--primary-busco-column",
                    "BUSCO_bacillota_odb12",
                    "--columns",
                    str(SAMPLE_STATUS_COLUMNS_ASSET),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            header, rows = read_tsv_rows(output)
            self.assertEqual(
                header,
                [
                    line.strip()
                    for line in SAMPLE_STATUS_COLUMNS_ASSET.read_text(encoding="utf-8").splitlines()
                    if line.strip()
                ],
            )
            by_accession = {row["accession"]: row for row in rows}

            self.assertEqual(by_accession["ACC1"]["taxonomy_status"], "done")
            self.assertEqual(by_accession["ACC1"]["barrnap_status"], "done")
            self.assertEqual(by_accession["ACC1"]["checkm2_gcode4_status"], "done")
            self.assertEqual(by_accession["ACC1"]["gcode"], "4")
            self.assertEqual(by_accession["ACC1"]["codetta_status"], "done")
            self.assertEqual(by_accession["ACC1"]["ccfinder_status"], "done")
            self.assertEqual(by_accession["ACC1"]["prokka_status"], "done")
            self.assertEqual(by_accession["ACC1"]["ani_included"], "true")
            self.assertEqual(
                by_accession["ACC1"]["warnings"],
                "internal_id_collision_resolved",
            )

            self.assertEqual(by_accession["ACC2"]["gcode_status"], "failed")
            self.assertEqual(by_accession["ACC2"]["gcode"], "NA")
            self.assertEqual(by_accession["ACC2"]["codetta_status"], "failed")
            self.assertEqual(by_accession["ACC2"]["prokka_status"], "skipped")
            self.assertEqual(by_accession["ACC2"]["ccfinder_status"], "skipped")
            self.assertEqual(by_accession["ACC2"]["ani_included"], "false")
            self.assertEqual(
                by_accession["ACC2"]["ani_exclusion_reason"],
                "gcode_na;no_16s;missing_primary_busco",
            )
            self.assertEqual(
                by_accession["ACC2"]["warnings"],
                "missing_metadata_for_new_sample;checkm2_gcode4_failed;busco_summary_failed;codetta_failed",
            )
            self.assertEqual(
                by_accession["ACC2"]["notes"],
                "New sample metadata missing; unavailable fields filled with NA.",
            )

    def test_main_accepts_case_insensitive_is_new_in_initial_status_seed(self) -> None:
        """Treat legacy TRUE/FALSE seed values as valid boolean equivalents."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id\nACC1\tfalse\tNA\t/path/one.fna\tid_1\n",
            )
            initial_status = self.write_tsv_rows(
                tmpdir / "initial_status.tsv",
                [
                    line.strip()
                    for line in SAMPLE_STATUS_COLUMNS_ASSET.read_text(encoding="utf-8").splitlines()
                    if line.strip()
                ],
                [
                    self.make_initial_status_row(
                        accession="ACC1",
                        internal_id="id_1",
                        is_new="FALSE",
                    )
                ],
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\tOrganism_Name\tAtypical_Warnings\nACC1\t123\tOne\tNA\n",
            )
            output = tmpdir / "sample_status.tsv"

            exit_code = build_sample_status.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--initial-status",
                    str(initial_status),
                    "--metadata",
                    str(metadata),
                    "--columns",
                    str(SAMPLE_STATUS_COLUMNS_ASSET),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            _header, rows = read_tsv_rows(output)
            self.assertEqual(rows[0]["is_new"], "false")

    def test_main_rejects_invalid_is_new_token_in_initial_status_seed(self) -> None:
        """Fail with an accession-specific error for invalid seed booleans."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_rows = [
                {
                    "accession": "ACC1",
                    "is_new": "false",
                    "assembly_level": "NA",
                    "genome_fasta": "/path/one.fna",
                    "internal_id": "id_1",
                }
            ]
            initial_status = self.write_tsv_rows(
                tmpdir / "initial_status.tsv",
                [
                    line.strip()
                    for line in SAMPLE_STATUS_COLUMNS_ASSET.read_text(encoding="utf-8").splitlines()
                    if line.strip()
                ],
                [
                    self.make_initial_status_row(
                        accession="ACC1",
                        internal_id="id_1",
                        is_new="maybe",
                    )
                ],
            )
            output_columns = [
                line.strip()
                for line in SAMPLE_STATUS_COLUMNS_ASSET.read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]

            with self.assertRaises(build_sample_status.SampleStatusError) as excinfo:
                build_sample_status.load_initial_status(
                    initial_status,
                    output_columns=output_columns,
                    validated_samples=validated_rows,
                )

            self.assertIn(
                "initial sample-status table contains invalid is_new value 'maybe' for accession 'ACC1'.",
                str(excinfo.exception),
            )

    def test_main_marks_missing_annotation_manifests_as_failed(self) -> None:
        """Fail gcode-qualified annotation statuses when manifest rows are missing."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id\nACC1\tfalse\tNA\t/path/one.fna\tid_1\n",
            )
            initial_status_rows = [
                self.make_initial_status_row(
                    accession="ACC1",
                    internal_id="id_1",
                    is_new="false",
                )
            ]
            initial_status = self.write_tsv_rows(
                tmpdir / "initial_status.tsv",
                [row for row in initial_status_rows[0]],
                initial_status_rows,
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\tOrganism_Name\tAtypical_Warnings\nACC1\t123\tOne\tNA\n",
            )
            checkm2 = self.write_text_file(
                tmpdir / "checkm2.tsv",
                "accession\tCompleteness_gcode4\tCompleteness_gcode11\tContamination_gcode4\tContamination_gcode11\tCoding_Density_gcode4\tCoding_Density_gcode11\tAverage_Gene_Length_gcode4\tAverage_Gene_Length_gcode11\tTotal_Coding_Sequences_gcode4\tTotal_Coding_Sequences_gcode11\tGcode\tGcode_Rule\tGcode_Length_Ratio\tGcode_Length_Ratio_Threshold\tGcode_Selection_Reason\tLow_quality\tcheckm2_status\twarnings\nACC1\t80\t95\t3\t1\t0.8\t0.95\t800\t1000\t700\t920\t11\tmean_gene_length_ratio\t0.8\t1.5\tlength_ratio_at_or_below_threshold\tfalse\tdone\t\n",
            )
            status_16s = self.write_text_file(
                tmpdir / "16s.tsv",
                "accession\t16S\tbest_16S_header\tbest_16S_length\twarnings\nACC1\tYes\th1\t1500\t\n",
            )
            busco = self.write_text_file(
                tmpdir / "busco.tsv",
                "accession\tlineage\tBUSCO_bacillota_odb12\tbusco_status\twarnings\nACC1\tbacillota_odb12\tC:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200\tdone\t\n",
            )
            ccfinder = self.write_text_file(
                tmpdir / "ccfinder.tsv",
                "accession\tCRISPRS\tSPACERS_SUM\tCRISPR_FRAC\tccfinder_status\twarnings\nACC1\t2\t7\t0.1\tdone\t\n",
            )
            prokka_manifest = self.write_text_file(
                tmpdir / "prokka_manifest.tsv",
                "accession\texit_code\tgff_size\tfaa_size\n",
            )
            output = tmpdir / "sample_status.tsv"

            exit_code = build_sample_status.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--initial-status",
                    str(initial_status),
                    "--metadata",
                    str(metadata),
                    "--checkm2",
                    str(checkm2),
                    "--16s-status",
                    str(status_16s),
                    "--busco",
                    str(busco),
                    "--ccfinder-strains",
                    str(ccfinder),
                    "--prokka-manifest",
                    str(prokka_manifest),
                    "--columns",
                    str(SAMPLE_STATUS_COLUMNS_ASSET),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            _, rows = read_tsv_rows(output)
            row = rows[0]
            self.assertEqual(row["prokka_status"], "failed")
            self.assertEqual(
                row["warnings"],
                "missing_prokka_result",
            )

    def test_main_marks_annotation_failures_from_logs_and_outputs(self) -> None:
        """Use exit codes and output presence to mark annotation failures."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id\nACC1\tfalse\tNA\t/path/one.fna\tid_1\n",
            )
            initial_status_rows = [
                self.make_initial_status_row(
                    accession="ACC1",
                    internal_id="id_1",
                    is_new="false",
                )
            ]
            initial_status = self.write_tsv_rows(
                tmpdir / "initial_status.tsv",
                [row for row in initial_status_rows[0]],
                initial_status_rows,
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\tOrganism_Name\tAtypical_Warnings\nACC1\t123\tOne\tNA\n",
            )
            checkm2 = self.write_text_file(
                tmpdir / "checkm2.tsv",
                "accession\tCompleteness_gcode4\tCompleteness_gcode11\tContamination_gcode4\tContamination_gcode11\tCoding_Density_gcode4\tCoding_Density_gcode11\tAverage_Gene_Length_gcode4\tAverage_Gene_Length_gcode11\tTotal_Coding_Sequences_gcode4\tTotal_Coding_Sequences_gcode11\tGcode\tGcode_Rule\tGcode_Length_Ratio\tGcode_Length_Ratio_Threshold\tGcode_Selection_Reason\tLow_quality\tcheckm2_status\twarnings\nACC1\t80\t95\t3\t1\t0.8\t0.95\t800\t1000\t700\t920\t11\tmean_gene_length_ratio\t0.8\t1.5\tlength_ratio_at_or_below_threshold\tfalse\tdone\t\n",
            )
            status_16s = self.write_text_file(
                tmpdir / "16s.tsv",
                "accession\t16S\tbest_16S_header\tbest_16S_length\twarnings\nACC1\tYes\th1\t1500\t\n",
            )
            busco = self.write_text_file(
                tmpdir / "busco.tsv",
                "accession\tlineage\tBUSCO_bacillota_odb12\tbusco_status\twarnings\nACC1\tbacillota_odb12\tC:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200\tdone\t\n",
            )
            ccfinder = self.write_text_file(
                tmpdir / "ccfinder.tsv",
                "accession\tCRISPRS\tSPACERS_SUM\tCRISPR_FRAC\tccfinder_status\twarnings\nACC1\t2\t7\t0.1\tdone\t\n",
            )
            prokka_manifest = self.write_text_file(
                tmpdir / "prokka_manifest.tsv",
                "accession\texit_code\tgff_size\tfaa_size\nACC1\t0\t0\t50\n",
            )
            output = tmpdir / "sample_status.tsv"

            exit_code = build_sample_status.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--initial-status",
                    str(initial_status),
                    "--metadata",
                    str(metadata),
                    "--checkm2",
                    str(checkm2),
                    "--16s-status",
                    str(status_16s),
                    "--busco",
                    str(busco),
                    "--ccfinder-strains",
                    str(ccfinder),
                    "--prokka-manifest",
                    str(prokka_manifest),
                    "--columns",
                    str(SAMPLE_STATUS_COLUMNS_ASSET),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            _, rows = read_tsv_rows(output)
            row = rows[0]
            self.assertEqual(row["prokka_status"], "failed")
            self.assertEqual(
                row["warnings"],
                "missing_prokka_outputs",
            )

    def test_main_uses_explicit_primary_busco_column_and_atypical_exception(self) -> None:
        """Use the configured primary BUSCO lineage and allow the locked exception."""
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
            initial_status_rows = [
                self.make_initial_status_row(accession="ACC1", internal_id="id_1", is_new="false"),
                self.make_initial_status_row(accession="ACC2", internal_id="id_2", is_new="false"),
            ]
            initial_status = self.write_tsv_rows(
                tmpdir / "initial_status.tsv",
                [row for row in initial_status_rows[0]],
                initial_status_rows,
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "\n".join(
                    [
                        "Accession\tTax_ID\tOrganism_Name\tAssembly_Level\tN50\tScaffolds\tGenome_Size\tAtypical_Warnings",
                        "ACC1\t123\tOne\tScaffold\tNA\tNA\tNA\tunverified source organism",
                        "ACC2\t456\tTwo\tContig\tNA\tNA\tNA\tmultiple source names",
                    ]
                )
                + "\n",
            )
            assembly_stats = self.write_text_file(
                tmpdir / "assembly_stats.tsv",
                "\n".join(
                    [
                        "accession\tn50\tscaffolds\tgenome_size\tgc_content",
                        "ACC1\t100000\t1\t100000\t50",
                        "ACC2\t90000\t2\t950000\t47.2",
                    ]
                )
                + "\n",
            )
            checkm2 = self.write_text_file(
                tmpdir / "checkm2.tsv",
                "\n".join(
                    [
                        "accession\tCompleteness_gcode4\tCompleteness_gcode11\tContamination_gcode4\tContamination_gcode11\tCoding_Density_gcode4\tCoding_Density_gcode11\tAverage_Gene_Length_gcode4\tAverage_Gene_Length_gcode11\tTotal_Coding_Sequences_gcode4\tTotal_Coding_Sequences_gcode11\tGcode\tGcode_Rule\tGcode_Length_Ratio\tGcode_Length_Ratio_Threshold\tGcode_Selection_Reason\tLow_quality\tcheckm2_status\twarnings",
                        "ACC1\t80\t95\t3\t1\t0.8\t0.95\t800\t1000\t700\t920\t11\tmean_gene_length_ratio\t0.8\t1.5\tlength_ratio_at_or_below_threshold\tfalse\tdone\t",
                        "ACC2\t80\t94\t3\t1\t0.8\t0.94\t800\t1000\t700\t910\t11\tmean_gene_length_ratio\t0.8\t1.5\tlength_ratio_at_or_below_threshold\tfalse\tdone\t",
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
                        "ACC2\tYes\th2\t1500\t",
                    ]
                )
                + "\n",
            )
            busco_bacillota = self.write_text_file(
                tmpdir / "busco_b.tsv",
                "accession\tlineage\tBUSCO_bacillota_odb12\tbusco_status\twarnings\nACC2\tbacillota_odb12\tC:97.0%[S:97.0%,D:0.0%],F:1.0%,M:2.0%,n:200\tdone\t\n",
            )
            busco_mycoplasmatota = self.write_text_file(
                tmpdir / "busco_m.tsv",
                "\n".join(
                    [
                        "accession\tlineage\tBUSCO_mycoplasmatota_odb12\tbusco_status\twarnings",
                        "ACC1\tmycoplasmatota_odb12\tC:99.0%[S:99.0%,D:0.0%],F:0.0%,M:1.0%,n:180\tdone\t",
                        "ACC2\tmycoplasmatota_odb12\tC:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:180\tdone\t",
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
                        "ACC2\t1\t2\t0.05\tdone\t",
                    ]
                )
                + "\n",
            )
            prokka_manifest = self.write_text_file(
                tmpdir / "prokka_manifest.tsv",
                "accession\texit_code\tgff_size\tfaa_size\nACC1\t0\t100\t50\nACC2\t0\t100\t50\n",
            )
            ani = self.write_text_file(
                tmpdir / "ani.tsv",
                "Accession\tCluster_ID\tIs_Representative\tANI_to_Representative\tScore\nACC1\tcluster_1\tyes\t100\t0.95\n",
            )
            output = tmpdir / "sample_status.tsv"

            exit_code = build_sample_status.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--initial-status",
                    str(initial_status),
                    "--metadata",
                    str(metadata),
                    "--checkm2",
                    str(checkm2),
                    "--16s-status",
                    str(status_16s),
                    "--busco",
                    str(busco_bacillota),
                    "--busco",
                    str(busco_mycoplasmatota),
                    "--ccfinder-strains",
                    str(ccfinder),
                    "--prokka-manifest",
                    str(prokka_manifest),
                    "--ani",
                    str(ani),
                    "--assembly-stats",
                    str(assembly_stats),
                    "--primary-busco-column",
                    "BUSCO_mycoplasmatota_odb12",
                    "--columns",
                    str(SAMPLE_STATUS_COLUMNS_ASSET),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            _, rows = read_tsv_rows(output)
            by_accession = {row["accession"]: row for row in rows}
            self.assertEqual(by_accession["ACC1"]["ani_included"], "true")
            self.assertEqual(by_accession["ACC1"]["ani_exclusion_reason"], "")
            self.assertEqual(by_accession["ACC2"]["ani_included"], "false")
            self.assertEqual(by_accession["ACC2"]["ani_exclusion_reason"], "atypical")

    def test_main_excludes_mixed_atypical_warnings_from_ani(self) -> None:
        """Exclude ANI samples when unverified source organism is not the sole reason."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id\n"
                "ACC1\tfalse\tNA\t/path/one.fna\tid_1\n",
            )
            initial_status_rows = [
                self.make_initial_status_row(accession="ACC1", internal_id="id_1", is_new="false"),
            ]
            initial_status = self.write_tsv_rows(
                tmpdir / "initial_status.tsv",
                [row for row in initial_status_rows[0]],
                initial_status_rows,
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "\n".join(
                    [
                        "Accession\tTax_ID\tOrganism_Name\tAssembly_Level\tN50\tScaffolds\tGenome_Size\tAtypical_Warnings",
                        "ACC1\t123\tOne\tScaffold\tNA\tNA\tNA\tgenome length too small, unverified source organism",
                    ]
                )
                + "\n",
            )
            assembly_stats = self.write_text_file(
                tmpdir / "assembly_stats.tsv",
                "accession\tn50\tscaffolds\tgenome_size\tgc_content\nACC1\t100000\t1\t100000\t50\n",
            )
            checkm2 = self.write_text_file(
                tmpdir / "checkm2.tsv",
                "\n".join(
                    [
                        "accession\tCompleteness_gcode4\tCompleteness_gcode11\tContamination_gcode4\tContamination_gcode11\tCoding_Density_gcode4\tCoding_Density_gcode11\tAverage_Gene_Length_gcode4\tAverage_Gene_Length_gcode11\tTotal_Coding_Sequences_gcode4\tTotal_Coding_Sequences_gcode11\tGcode\tGcode_Rule\tGcode_Length_Ratio\tGcode_Length_Ratio_Threshold\tGcode_Selection_Reason\tLow_quality\tcheckm2_status\twarnings",
                        "ACC1\t80\t95\t3\t1\t0.8\t0.95\t800\t1000\t700\t920\t11\tmean_gene_length_ratio\t0.8\t1.5\tlength_ratio_at_or_below_threshold\tfalse\tdone\t",
                    ]
                )
                + "\n",
            )
            status_16s = self.write_text_file(
                tmpdir / "16s.tsv",
                "accession\t16S\tbest_16S_header\tbest_16S_length\twarnings\nACC1\tYes\th1\t1500\t\n",
            )
            busco = self.write_text_file(
                tmpdir / "busco.tsv",
                "accession\tlineage\tBUSCO_bacillota_odb12\tbusco_status\twarnings\nACC1\tbacillota_odb12\tC:99.0%[S:99.0%,D:0.0%],F:0.0%,M:1.0%,n:180\tdone\t\n",
            )
            ani = self.write_text_file(
                tmpdir / "ani.tsv",
                "Accession\tCluster_ID\tIs_Representative\tANI_to_Representative\tScore\n",
            )
            output = tmpdir / "sample_status.tsv"

            exit_code = build_sample_status.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--initial-status",
                    str(initial_status),
                    "--metadata",
                    str(metadata),
                    "--checkm2",
                    str(checkm2),
                    "--16s-status",
                    str(status_16s),
                    "--busco",
                    str(busco),
                    "--ani",
                    str(ani),
                    "--assembly-stats",
                    str(assembly_stats),
                    "--primary-busco-column",
                    "BUSCO_bacillota_odb12",
                    "--columns",
                    str(SAMPLE_STATUS_COLUMNS_ASSET),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            _, rows = read_tsv_rows(output)
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["ani_included"], "false")
            self.assertEqual(rows[0]["ani_exclusion_reason"], "atypical")

    def test_main_rejects_initial_status_identity_mismatch(self) -> None:
        """Hard-fail when the seed status table disagrees with validated samples."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id\nACC1\tfalse\tNA\t/path/one.fna\tid_1\n",
            )
            initial_status_rows = [
                self.make_initial_status_row(
                    accession="ACC1",
                    internal_id="other_id",
                    is_new="false",
                )
            ]
            initial_status = self.write_tsv_rows(
                tmpdir / "initial_status.tsv",
                [row for row in initial_status_rows[0]],
                initial_status_rows,
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\tOrganism_Name\tAtypical_Warnings\nACC1\t123\tOne\tNA\n",
            )
            output = tmpdir / "sample_status.tsv"

            exit_code = build_sample_status.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--initial-status",
                    str(initial_status),
                    "--metadata",
                    str(metadata),
                    "--columns",
                    str(SAMPLE_STATUS_COLUMNS_ASSET),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse(output.exists())

    def test_main_populates_custom_busco_status_columns(self) -> None:
        """Use configured BUSCO lineages for non-default sample-status columns."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            custom_columns = [
                "accession",
                "internal_id",
                "is_new",
                "validation_status",
                "taxonomy_status",
                "barrnap_status",
                "checkm2_gcode4_status",
                "checkm2_gcode11_status",
                "gcode_status",
                "gcode",
                "low_quality",
                "busco_custom_odb12_status",
                "codetta_status",
                "prokka_status",
                "ccfinder_status",
                "padloc_status",
                "eggnog_status",
                "cogclassifier_status",
                "pfam_status",
                "kofam_status",
                "annotation_complete",
                "ani_included",
                "ani_exclusion_reason",
                "warnings",
                "notes",
            ]
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id\n"
                "ACC1\tfalse\tNA\t/path/one.fna\tid_1\n",
            )
            initial_status_rows = [
                self.make_initial_status_row(
                    accession="ACC1",
                    internal_id="id_1",
                    is_new="false",
                    columns=custom_columns,
                )
            ]
            initial_status = self.write_tsv_rows(
                tmpdir / "initial_status.tsv",
                custom_columns,
                initial_status_rows,
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\tOrganism_Name\tAssembly_Level\tN50\tScaffolds\tGenome_Size\tAtypical_Warnings\n"
                "ACC1\t123\tKnown one\tScaffold\t50000\t2\t800000\tNA\n",
            )
            checkm2 = self.write_text_file(
                tmpdir / "checkm2.tsv",
                "accession\tCompleteness_gcode4\tCompleteness_gcode11\tContamination_gcode4\tContamination_gcode11\tCoding_Density_gcode4\tCoding_Density_gcode11\tAverage_Gene_Length_gcode4\tAverage_Gene_Length_gcode11\tTotal_Coding_Sequences_gcode4\tTotal_Coding_Sequences_gcode11\tGcode\tGcode_Rule\tGcode_Length_Ratio\tGcode_Length_Ratio_Threshold\tGcode_Selection_Reason\tLow_quality\tcheckm2_status\twarnings\n"
                "ACC1\t80\t95\t3\t1\t0.8\t0.95\t800\t1000\t700\t920\t11\tmean_gene_length_ratio\t0.8\t1.5\tlength_ratio_at_or_below_threshold\tfalse\tdone\t\n",
            )
            status_16s = self.write_text_file(
                tmpdir / "16s.tsv",
                "accession\t16S\tbest_16S_header\tbest_16S_length\twarnings\n"
                "ACC1\tYes\th1\t1500\t\n",
            )
            busco_custom = self.write_text_file(
                tmpdir / "busco_custom.tsv",
                "accession\tlineage\tBUSCO_custom_odb12\tbusco_status\twarnings\n"
                "ACC1\tcustom_odb12\tC:99.0%[S:99.0%,D:0.0%],F:0.0%,M:1.0%,n:180\tdone\t\n",
            )
            ani = self.write_text_file(
                tmpdir / "ani.tsv",
                "Accession\tCluster_ID\tIs_Representative\tANI_to_Representative\tScore\n"
                "ACC1\tcluster_1\tyes\t100\t0.95\n",
            )
            output = tmpdir / "sample_status.tsv"

            exit_code = build_sample_status.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--initial-status",
                    str(initial_status),
                    "--metadata",
                    str(metadata),
                    "--checkm2",
                    str(checkm2),
                    "--16s-status",
                    str(status_16s),
                    "--busco",
                    str(busco_custom),
                    "--ani",
                    str(ani),
                    "--primary-busco-column",
                    "BUSCO_custom_odb12",
                    "--busco-lineage",
                    "custom_odb12",
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            header, rows = read_tsv_rows(output)
            self.assertIn("busco_custom_odb12_status", header)
            self.assertNotIn("busco_bacillota_odb12_status", header)
            self.assertEqual(rows[0]["busco_custom_odb12_status"], "done")
            self.assertEqual(rows[0]["ani_included"], "true")

    def test_main_uses_in_house_assembly_stats_for_ani_decision(self) -> None:
        """Treat computed stats as authoritative when metadata metrics are missing."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id\nACC1\ttrue\tComplete Genome\t/path/one.fna\tid_1\n",
            )
            initial_status_rows = [
                self.make_initial_status_row(
                    accession="ACC1",
                    internal_id="id_1",
                    is_new="true",
                    warnings="missing_metadata_for_new_sample",
                )
            ]
            initial_status = self.write_tsv_rows(
                tmpdir / "initial_status.tsv",
                [row for row in initial_status_rows[0]],
                initial_status_rows,
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\tOrganism_Name\tN50\tScaffolds\tGenome_Size\tAtypical_Warnings\n"
                "ACC1\t123\tOne\tNA\tNA\tNA\tNA\n",
            )
            assembly_stats = self.write_text_file(
                tmpdir / "assembly_stats.tsv",
                "accession\tn50\tscaffolds\tgenome_size\tgc_content\nACC1\t100000\t1\t100000\t50\n",
            )
            checkm2 = self.write_text_file(
                tmpdir / "checkm2.tsv",
                "accession\tCompleteness_gcode4\tCompleteness_gcode11\tContamination_gcode4\tContamination_gcode11\tCoding_Density_gcode4\tCoding_Density_gcode11\tAverage_Gene_Length_gcode4\tAverage_Gene_Length_gcode11\tTotal_Coding_Sequences_gcode4\tTotal_Coding_Sequences_gcode11\tGcode\tGcode_Rule\tGcode_Length_Ratio\tGcode_Length_Ratio_Threshold\tGcode_Selection_Reason\tLow_quality\tcheckm2_status\twarnings\nACC1\t80\t95\t3\t1\t0.8\t0.95\t800\t1000\t700\t920\t11\tmean_gene_length_ratio\t0.8\t1.5\tlength_ratio_at_or_below_threshold\tfalse\tdone\t\n",
            )
            status_16s = self.write_text_file(
                tmpdir / "16s.tsv",
                "accession\t16S\tbest_16S_header\tbest_16S_length\twarnings\nACC1\tYes\th1\t1500\t\n",
            )
            busco = self.write_text_file(
                tmpdir / "busco.tsv",
                "accession\tlineage\tBUSCO_bacillota_odb12\tbusco_status\twarnings\nACC1\tbacillota_odb12\tC:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200\tdone\t\n",
            )
            ani = self.write_text_file(
                tmpdir / "ani.tsv",
                "Accession\tCluster_ID\tIs_Representative\tANI_to_Representative\tScore\nACC1\tcluster_1\tyes\t100\t0.95\n",
            )
            output = tmpdir / "sample_status.tsv"

            exit_code = build_sample_status.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--initial-status",
                    str(initial_status),
                    "--metadata",
                    str(metadata),
                    "--checkm2",
                    str(checkm2),
                    "--16s-status",
                    str(status_16s),
                    "--busco",
                    str(busco),
                    "--ani",
                    str(ani),
                    "--assembly-stats",
                    str(assembly_stats),
                    "--primary-busco-column",
                    "BUSCO_bacillota_odb12",
                    "--columns",
                    str(SAMPLE_STATUS_COLUMNS_ASSET),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            _, rows = read_tsv_rows(output)
            self.assertEqual(rows[0]["ani_included"], "true")
            self.assertEqual(rows[0]["ani_exclusion_reason"], "")

    def test_main_still_fails_when_eligible_sample_is_missing_from_ani_summary(self) -> None:
        """Keep the hard-fail for truly eligible samples missing from ANI output."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id\nACC1\ttrue\tComplete Genome\t/path/one.fna\tid_1\n",
            )
            initial_status_rows = [
                self.make_initial_status_row(
                    accession="ACC1",
                    internal_id="id_1",
                    is_new="true",
                )
            ]
            initial_status = self.write_tsv_rows(
                tmpdir / "initial_status.tsv",
                [row for row in initial_status_rows[0]],
                initial_status_rows,
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\tOrganism_Name\tN50\tScaffolds\tGenome_Size\tAtypical_Warnings\n"
                "ACC1\t123\tOne\tNA\tNA\tNA\tNA\n",
            )
            assembly_stats = self.write_text_file(
                tmpdir / "assembly_stats.tsv",
                "accession\tn50\tscaffolds\tgenome_size\tgc_content\nACC1\t100000\t1\t100000\t50\n",
            )
            checkm2 = self.write_text_file(
                tmpdir / "checkm2.tsv",
                "accession\tCompleteness_gcode4\tCompleteness_gcode11\tContamination_gcode4\tContamination_gcode11\tCoding_Density_gcode4\tCoding_Density_gcode11\tAverage_Gene_Length_gcode4\tAverage_Gene_Length_gcode11\tTotal_Coding_Sequences_gcode4\tTotal_Coding_Sequences_gcode11\tGcode\tGcode_Rule\tGcode_Length_Ratio\tGcode_Length_Ratio_Threshold\tGcode_Selection_Reason\tLow_quality\tcheckm2_status\twarnings\nACC1\t80\t95\t3\t1\t0.8\t0.95\t800\t1000\t700\t920\t11\tmean_gene_length_ratio\t0.8\t1.5\tlength_ratio_at_or_below_threshold\tfalse\tdone\t\n",
            )
            status_16s = self.write_text_file(
                tmpdir / "16s.tsv",
                "accession\t16S\tbest_16S_header\tbest_16S_length\twarnings\nACC1\tYes\th1\t1500\t\n",
            )
            busco = self.write_text_file(
                tmpdir / "busco.tsv",
                "accession\tlineage\tBUSCO_bacillota_odb12\tbusco_status\twarnings\nACC1\tbacillota_odb12\tC:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200\tdone\t\n",
            )
            ani = self.write_text_file(
                tmpdir / "ani.tsv",
                "Accession\tCluster_ID\tIs_Representative\tANI_to_Representative\tScore\n",
            )
            output = tmpdir / "sample_status.tsv"

            exit_code = build_sample_status.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--initial-status",
                    str(initial_status),
                    "--metadata",
                    str(metadata),
                    "--checkm2",
                    str(checkm2),
                    "--16s-status",
                    str(status_16s),
                    "--busco",
                    str(busco),
                    "--ani",
                    str(ani),
                    "--assembly-stats",
                    str(assembly_stats),
                    "--primary-busco-column",
                    "BUSCO_bacillota_odb12",
                    "--columns",
                    str(SAMPLE_STATUS_COLUMNS_ASSET),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 1)
            self.assertFalse(output.exists())

    def test_main_excludes_missing_assembly_level_without_hard_fail(self) -> None:
        """Treat missing ANI assembly level as an exclusion, not a missing-ANI crash."""
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmpdir = Path(tmpdir_name)
            validated_samples = self.write_text_file(
                tmpdir / "validated_samples.tsv",
                "accession\tis_new\tassembly_level\tgenome_fasta\tinternal_id\nACC1\tfalse\tNA\t/path/one.fna\tid_1\n",
            )
            initial_status_rows = [
                self.make_initial_status_row(
                    accession="ACC1",
                    internal_id="id_1",
                    is_new="false",
                )
            ]
            initial_status = self.write_tsv_rows(
                tmpdir / "initial_status.tsv",
                [row for row in initial_status_rows[0]],
                initial_status_rows,
            )
            metadata = self.write_text_file(
                tmpdir / "metadata.tsv",
                "Accession\tTax_ID\tOrganism_Name\tAssembly_Level\tN50\tScaffolds\tGenome_Size\tAtypical_Warnings\n"
                "ACC1\t123\tOne\tNA\t50000\t2\t800000\tNA\n",
            )
            checkm2 = self.write_text_file(
                tmpdir / "checkm2.tsv",
                "accession\tCompleteness_gcode4\tCompleteness_gcode11\tContamination_gcode4\tContamination_gcode11\tCoding_Density_gcode4\tCoding_Density_gcode11\tAverage_Gene_Length_gcode4\tAverage_Gene_Length_gcode11\tTotal_Coding_Sequences_gcode4\tTotal_Coding_Sequences_gcode11\tGcode\tGcode_Rule\tGcode_Length_Ratio\tGcode_Length_Ratio_Threshold\tGcode_Selection_Reason\tLow_quality\tcheckm2_status\twarnings\nACC1\t80\t95\t3\t1\t0.8\t0.95\t800\t1000\t700\t920\t11\tmean_gene_length_ratio\t0.8\t1.5\tlength_ratio_at_or_below_threshold\tfalse\tdone\t\n",
            )
            status_16s = self.write_text_file(
                tmpdir / "16s.tsv",
                "accession\t16S\tbest_16S_header\tbest_16S_length\twarnings\nACC1\tYes\th1\t1500\t\n",
            )
            busco = self.write_text_file(
                tmpdir / "busco.tsv",
                "accession\tlineage\tBUSCO_bacillota_odb12\tbusco_status\twarnings\nACC1\tbacillota_odb12\tC:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200\tdone\t\n",
            )
            ani = self.write_text_file(
                tmpdir / "ani.tsv",
                "Accession\tCluster_ID\tIs_Representative\tANI_to_Representative\tScore\n",
            )
            output = tmpdir / "sample_status.tsv"

            exit_code = build_sample_status.main(
                [
                    "--validated-samples",
                    str(validated_samples),
                    "--initial-status",
                    str(initial_status),
                    "--metadata",
                    str(metadata),
                    "--checkm2",
                    str(checkm2),
                    "--16s-status",
                    str(status_16s),
                    "--busco",
                    str(busco),
                    "--ani",
                    str(ani),
                    "--primary-busco-column",
                    "BUSCO_bacillota_odb12",
                    "--columns",
                    str(SAMPLE_STATUS_COLUMNS_ASSET),
                    "--output",
                    str(output),
                ]
            )

            self.assertEqual(exit_code, 0)
            _, rows = read_tsv_rows(output)
            self.assertEqual(rows[0]["ani_included"], "false")
            self.assertEqual(rows[0]["ani_exclusion_reason"], "missing_assembly_level")


if __name__ == "__main__":
    unittest.main()
