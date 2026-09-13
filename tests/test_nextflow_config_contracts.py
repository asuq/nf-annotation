"""Regression checks for Nextflow runtime configuration contracts."""

from __future__ import annotations

import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW_CONFIG = ROOT / "nextflow.config"
MAIN_WORKFLOW = ROOT / "main.nf"
PREP_WORKFLOW = ROOT / "prepare_databases.nf"
DOCKER_CONFIG = ROOT / "conf" / "docker.config"
SINGULARITY_CONFIG = ROOT / "conf" / "singularity.config"
SLURM_CONFIG = ROOT / "conf" / "slurm.config"
OIST_CONFIG = ROOT / "conf" / "oist.config"
OIST_20K_STORAGE_CONFIG = ROOT / "conf" / "oist_20k_storage.config"
GWDG_CONFIG = ROOT / "conf" / "gwdg.config"
MARMIC_CONFIG = ROOT / "conf" / "marmic.config"
VIPER_CPU_CONFIG = ROOT / "conf" / "viper-cpu.config"
LOCAL_CONFIG = ROOT / "conf" / "local.config"
BASE_CONFIG = ROOT / "conf" / "base.config"


class NextflowConfigContractsTestCase(unittest.TestCase):
    """Protect configuration needed for containerised local runs."""

    def test_nextflow_defaults_define_runtime_params(self) -> None:
        """Keep warning-prone runtime params defined in the root config."""
        config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")

        self.assertIn("ani_score_profile = 'default'", config_text)
        self.assertIn("ani_16s_policy = 'complete'", config_text)
        self.assertIn("busco_db = null", config_text)
        self.assertIn("codetta_db = null", config_text)
        self.assertIn("codetta_db_label = null", config_text)
        self.assertIn("codetta_extra_args = ''", config_text)
        self.assertIn("download_missing_databases = false", config_text)
        self.assertIn("force_runtime_database_rebuild = false", config_text)
        self.assertIn("gcode_rule = 'mean_gene_length_ratio'", config_text)
        self.assertIn("runtime_db_helper_container = 'quay.io/asuq1617/nf-myco_db:0.3'", config_text)
        self.assertIn("runtime_db_scratch_root = null", config_text)
        self.assertIn("task_attempts = 3", config_text)
        self.assertIn("soft_fail_attempts = 3", config_text)
        self.assertIn("taxdump_version = null", config_text)
        self.assertIn("db_download_attempts = 2", config_text)
        self.assertIn("annotation_tools = 'eggnog,cogclassifier,pfam,kofam,padloc'", config_text)
        self.assertIn("singularity_cache_dir = null", config_text)
        self.assertIn("singularity_run_options = ''", config_text)
        self.assertIn("slurm_qos = null", config_text)
        self.assertNotIn("includeConfig 'conf/debug.config'", config_text)
        self.assertIn("includeConfig 'conf/oist.config'", config_text)
        self.assertIn("includeConfig 'conf/gwdg.config'", config_text)
        self.assertIn("includeConfig 'conf/marmic.config'", config_text)
        self.assertIn("includeConfig 'conf/viper-cpu.config'", config_text)
        self.assertNotIn("padloc_db = null", config_text)
        self.assertNotIn("padloc_db_label = null", config_text)
        self.assertNotIn("slurm_account", config_text)
        self.assertNotIn("use_biocontainers", config_text)
        self.assertNotIn("apptainer_cache_dir", config_text)
        self.assertNotIn("apptainer_run_options", config_text)
        self.assertNotIn("max_retries = ", config_text)

    def test_singularity_runtime_surface_uses_renamed_config_contract(self) -> None:
        """Keep the renamed Singularity profile, include, and config file wired in."""
        config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")
        singularity_text = SINGULARITY_CONFIG.read_text(encoding="utf-8")
        docker_text = DOCKER_CONFIG.read_text(encoding="utf-8")
        slurm_text = SLURM_CONFIG.read_text(encoding="utf-8")

        self.assertTrue(SINGULARITY_CONFIG.is_file())
        self.assertFalse((ROOT / "conf" / "apptainer.config").exists())
        self.assertIn("includeConfig 'conf/singularity.config'", config_text)
        self.assertNotIn("includeConfig 'conf/apptainer.config'", config_text)
        self.assertIn("singularity {", config_text)
        self.assertNotIn("apptainer {", config_text)
        self.assertIn("profiles {", singularity_text)
        self.assertIn("singularity.enabled = true", singularity_text)
        self.assertIn("singularity.autoMounts = true", singularity_text)
        self.assertIn("params.singularity_cache_dir", singularity_text)
        self.assertIn("params.singularity_run_options", singularity_text)
        self.assertIn("singularity.enabled = false", docker_text)
        self.assertNotIn("params.slurm_account", slurm_text)
        self.assertNotIn("apptainer.enabled = false", docker_text)

    def test_oist_profile_layers_project_overrides_on_nf_helper(self) -> None:
        """Keep OIST site defaults shared while preserving annotation limits."""
        config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")
        oist_text = OIST_CONFIG.read_text(encoding="utf-8")
        shared_oist = ROOT / "external" / "nf-helper" / "conf" / "sites" / "oist.config"

        self.assertTrue(OIST_CONFIG.is_file())
        self.assertTrue(shared_oist.is_file())
        self.assertIn("includeConfig 'conf/oist.config'", config_text)
        self.assertIn('includeConfig "${projectDir}/external/nf-helper/conf/sites/oist.config"', oist_text)
        self.assertIn("profiles {", oist_text)
        self.assertIn("oist {", oist_text)
        self.assertIn("max_memory = 500.GB", oist_text)
        self.assertIn("max_time = 2.h", oist_text)
        self.assertIn("queue = params.slurm_queue ?: 'short'", oist_text)
        self.assertIn("containerOptions = params.singularity_run_options ?: ''", oist_text)
        self.assertIn("withLabel: process_medium", oist_text)
        self.assertIn("withLabel: process_high", oist_text)
        self.assertNotIn("def buildSlurmClusterOptions", oist_text)
        self.assertNotIn("beforeScript", oist_text)

    def test_gwdg_profile_layers_project_overrides_on_nf_helper(self) -> None:
        """Keep GWDG site defaults shared while preserving annotation selectors."""
        config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")
        gwdg_text = GWDG_CONFIG.read_text(encoding="utf-8")
        shared_gwdg = ROOT / "external" / "nf-helper" / "conf" / "sites" / "gwdg.config"

        self.assertTrue(GWDG_CONFIG.is_file())
        self.assertTrue(shared_gwdg.is_file())
        self.assertIn("includeConfig 'conf/gwdg.config'", config_text)
        self.assertIn('includeConfig "${projectDir}/external/nf-helper/conf/sites/gwdg.config"', gwdg_text)
        self.assertIn("profiles {", gwdg_text)
        self.assertIn("gwdg {", gwdg_text)
        self.assertIn("withLabel: process_medium", gwdg_text)
        self.assertIn("withLabel: process_high", gwdg_text)
        self.assertIn("cpus = { Math.min(16 * task.attempt, params.max_cpus as int) }", gwdg_text)
        self.assertIn("memory = { [50.GB * task.attempt, params.max_memory].min() }", gwdg_text)
        self.assertIn("time = { [6.h * task.attempt, params.max_time].min() }", gwdg_text)
        self.assertIn("cpus = { Math.min(64 * task.attempt, params.max_cpus as int) }", gwdg_text)
        self.assertIn("memory = { [200.GB * task.attempt, params.max_memory].min() }", gwdg_text)
        self.assertIn("time = { [1.d * task.attempt, params.max_time].min() }", gwdg_text)
        self.assertIn("withName: PROKKA", gwdg_text)
        self.assertIn("container = { params.prokka_container }", gwdg_text)
        self.assertIn("withName: CALCULATE_ASSEMBLY_STATS", gwdg_text)
        self.assertIn("container = { params.seqtk_container }", gwdg_text)
        self.assertIn("stageInMode = 'copy'", gwdg_text)
        self.assertNotIn("def buildSlurmClusterOptions", gwdg_text)
        self.assertNotIn("params.slurm_account", gwdg_text)

    def test_marmic_profile_comes_from_nf_helper(self) -> None:
        """Expose Marmic as a reusable shared site profile."""
        config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")
        marmic_text = MARMIC_CONFIG.read_text(encoding="utf-8")
        shared_marmic = ROOT / "external" / "nf-helper" / "conf" / "sites" / "marmic.config"

        self.assertTrue(MARMIC_CONFIG.is_file())
        self.assertTrue(shared_marmic.is_file())
        self.assertIn("includeConfig 'conf/marmic.config'", config_text)
        self.assertIn('includeConfig "${projectDir}/external/nf-helper/conf/sites/marmic.config"', marmic_text)
        self.assertNotIn("assembly", marmic_text)
        self.assertNotIn("trace {", marmic_text)

    def test_viper_cpu_profile_maps_database_downloads_to_the_login_node(self) -> None:
        """Keep Viper site defaults shared and all internet selectors explicit."""
        config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")
        viper_text = VIPER_CPU_CONFIG.read_text(encoding="utf-8")
        shared_viper = ROOT / "external" / "nf-helper" / "conf" / "sites" / "viper-cpu.config"

        self.assertTrue(VIPER_CPU_CONFIG.is_file())
        self.assertTrue(shared_viper.is_file())
        self.assertIn("includeConfig 'conf/viper-cpu.config'", config_text)
        self.assertIn(
            'includeConfig "${projectDir}/external/nf-helper/conf/sites/viper-cpu.config"',
            viper_text,
        )
        self.assertIn("'viper-cpu' {", viper_text)
        for process_name in (
            "DOWNLOAD_BUSCO_DATASET",
            "PREP_TAXDUMP_DATABASE",
            "PREP_CODETTA_DATABASE",
            "DOWNLOAD_CHECKM2_DATABASE",
            "DOWNLOAD_BUSCO_DATABASES",
            "PREP_ANNOTATION_RESOURCE",
        ):
            self.assertIn(process_name, viper_text)
        self.assertIn("executor = 'local'", viper_text)
        self.assertIn("cpus = 1", viper_text)
        self.assertIn("memory = 16.GB", viper_text)

    def test_oist_20k_storage_override_is_opt_in_and_resume_safe(self) -> None:
        """Keep the large-cohort OIST override separate from the default tracked profiles."""
        config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")
        override_text = OIST_20K_STORAGE_CONFIG.read_text(encoding="utf-8")

        self.assertTrue(OIST_20K_STORAGE_CONFIG.is_file())
        self.assertNotIn("includeConfig 'conf/oist_20k_storage.config'", config_text)
        self.assertIn("scratch = '/scratch'", override_text)
        self.assertIn("stageOutMode = 'move'", override_text)
        self.assertIn("withName: BUILD_FASTANI_INPUTS", override_text)
        self.assertIn("scratch = false", override_text)
        self.assertIn("withName: FASTANI", override_text)
        self.assertIn("stageInMode = 'copy'", override_text)
        self.assertNotIn("cleanup = true", override_text)
        self.assertNotIn("workDir =", override_text)

    def test_python_container_uses_shared_repo_owned_helper_image(self) -> None:
        """Use one shared helper image that carries the ANI scientific stack."""
        config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")

        self.assertIn("python_container = 'quay.io/asuq1617/python-scipy:3.12'", config_text)
        self.assertNotIn("python_container = 'python:3.12'", config_text)

    def test_ccfinder_container_points_at_the_clean_runtime_tag(self) -> None:
        """Use the cleaned CRISPRCasFinder image tag by default."""
        config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")

        self.assertIn("ccfinder_container = 'quay.io/asuq1617/ccfinder:4.2.30'", config_text)
        self.assertNotIn("ccfinder_container = 'quay.io/asuq1617/ccfinder:4.3.2'", config_text)

    def test_native_annotation_containers_require_immutable_explicit_references(self) -> None:
        text = NEXTFLOW_CONFIG.read_text()
        for tool in ("eggnog", "cogclassifier", "pfam", "kofam", "padloc"):
            self.assertIn(f"{tool}_container = null", text)

    def test_prokka_container_points_at_the_fixed_runtime_tag(self) -> None:
        """Use the fixed Prokka image tag by default."""
        config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")

        self.assertIn("prokka_container = 'quay.io/asuq1617/prokka:1.15.6'", config_text)
        self.assertNotIn(
            "prokka_container = 'quay.io/biocontainers/prokka:1.15.6--pl5321hdfd78af_0'",
            config_text,
        )

    def test_eggnog_v2_default_is_removed(self) -> None:
        """Reject the previous mapper-v2 default."""
        config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")

        self.assertIn("eggnog_container = null", config_text)
        self.assertNotIn(
            "eggnog_container = 'quay.io/biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_2'",
            config_text,
        )

    def test_codetta_container_points_at_the_fixed_runtime_tag(self) -> None:
        """Use the fixed Codetta image tag by default."""
        config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")

        self.assertIn("codetta_container = 'quay.io/asuq1617/codetta:2.0'", config_text)

    def test_python_helper_processes_are_pinned_explicitly(self) -> None:
        """Keep shared helper processes on the single Python helper image."""
        config_text = BASE_CONFIG.read_text(encoding="utf-8")

        self.assertNotIn("params.standardMaxRetries", config_text)
        self.assertNotIn("def retryThenFinish", config_text)
        self.assertNotIn("soft_fail_attempts", config_text)
        self.assertNotIn("errorStrategy = 'retry'", config_text)
        self.assertIn("errorStrategy = {", config_text)
        self.assertIn(
            "task.attempt <= Math.max((params.task_attempts as int) - 1, 0)",
            config_text,
        )
        self.assertIn(
            "maxRetries = { Math.max((params.task_attempts as int) - 1, 0) }",
            config_text,
        )
        self.assertIn("withName: CLUSTER_ANI {\n        container = { params.python_container }", config_text)
        self.assertIn("withName: SELECT_ANI_REPRESENTATIVES {\n        container = { params.python_container }", config_text)
        self.assertIn("withName: SUMMARISE_16S {\n        container = { params.python_container }", config_text)
        self.assertIn("withName: WRITE_SAMPLE_STATUS {\n        container = { params.python_container }", config_text)

    def test_seqtk_helper_processes_are_pinned_explicitly(self) -> None:
        """Keep seqtk-backed helper processes on the shared seqtk image."""
        config_text = BASE_CONFIG.read_text(encoding="utf-8")

        self.assertIn("withName: STAGE_INPUTS {\n        container = { params.seqtk_container }", config_text)
        self.assertIn("withName: CALCULATE_ASSEMBLY_STATS {\n        container = { params.seqtk_container }", config_text)

    def test_runtime_database_prep_processes_use_the_dedicated_helper_image(self) -> None:
        """Keep the prep entry point on the dedicated runtime DB helper image."""
        root_config_text = NEXTFLOW_CONFIG.read_text(encoding="utf-8")
        config_text = BASE_CONFIG.read_text(encoding="utf-8")

        self.assertNotIn("params.dbDownloadMaxRetries", config_text)
        self.assertIn(
            "task.attempt <= Math.max((params.db_download_attempts as int) - 1, 0)",
            config_text,
        )
        self.assertIn(
            "maxRetries = { Math.max((params.db_download_attempts as int) - 1, 0) }",
            config_text,
        )
        self.assertIn("withName: DOWNLOAD_BUSCO_DATASET {\n        errorStrategy = {", config_text)
        self.assertIn(
            "withName: DOWNLOAD_BUSCO_DATASET {\n        errorStrategy = {\n            task.attempt <= Math.max((params.db_download_attempts as int) - 1, 0)",
            config_text,
        )
        self.assertIn("withLabel: prep_taxdump_database {\n        errorStrategy = {", config_text)
        self.assertIn(
            "withLabel: prep_codetta_database {\n        errorStrategy = {",
            config_text,
        )
        self.assertIn(
            "withLabel: finalise_runtime_database {\n        container = { params.runtime_db_helper_container }",
            config_text,
        )
        self.assertIn(
            "withLabel: merge_runtime_database_reports {\n        container = { params.runtime_db_helper_container }",
            config_text,
        )
        self.assertNotIn("quay.io/asuq1617/nf-myco_db:0.3", config_text)
        self.assertNotIn("quay.io/asuq1617/nf-myco_db:0.2", config_text)
        self.assertNotIn("quay.io/asuq1617/nf-myco_db:0.2", root_config_text)
        self.assertIn(
            "withLabel: download_checkm2_database {\n        errorStrategy = {",
            config_text,
        )
        self.assertIn(
            "withLabel: download_busco_databases {\n        errorStrategy = {",
            config_text,
        )
        self.assertNotIn("withName: PREP_TAXDUMP_DATABASE", config_text)
        self.assertNotIn("withName: DOWNLOAD_CHECKM2_DATABASE", config_text)
        self.assertNotIn("withName: DOWNLOAD_BUSCO_DATABASES", config_text)
        self.assertNotIn("withName: DOWNLOAD_EGGNOG_DATABASE", config_text)
        self.assertNotIn("withName: DOWNLOAD_PADLOC_DATABASE", config_text)
        self.assertNotIn("withName: FINALISE_RUNTIME_DATABASE", config_text)
        self.assertNotIn("withName: MERGE_RUNTIME_DATABASE_REPORTS", config_text)

    def test_main_and_prep_workflows_require_codetta_runtime_path(self) -> None:
        """Keep the Codetta runtime database wired into both workflow entrypoints."""
        main_text = MAIN_WORKFLOW.read_text(encoding="utf-8")
        prep_text = PREP_WORKFLOW.read_text(encoding="utf-8")

        self.assertIn("if (!params.codetta_db) {", main_text)
        self.assertIn("error \"params.codetta_db is required.\"", main_text)
        self.assertIn("destinations.codetta", prep_text)
        self.assertIn("codettaRequest", prep_text)

    def test_docker_profile_forces_amd64_ccfinder_on_arm_hosts(self) -> None:
        """Allow Apple Silicon Docker runs to pull the pinned CCFINDER image."""
        config_text = DOCKER_CONFIG.read_text(encoding="utf-8")

        self.assertIn("withName: CCFINDER", config_text)
        self.assertIn("container = { params.ccfinder_container }", config_text)
        self.assertIn("--platform linux/amd64", config_text)
        self.assertIn("['aarch64', 'arm64']", config_text)

    def test_slurm_profile_inherits_shared_error_strategy(self) -> None:
        """Keep SLURM runtime profiles on the shared retry-then-finish policy."""
        slurm_text = SLURM_CONFIG.read_text(encoding="utf-8")

        self.assertNotIn("def buildSlurmClusterOptions", slurm_text)
        self.assertIn("process.clusterOptions = {", slurm_text)
        self.assertIn("params.slurm_cluster_options instanceof Boolean", slurm_text)
        self.assertIn("Use --slurm_qos 2h or --slurm_cluster_options='--qos=2h'.", slurm_text)
        self.assertIn('params.slurm_qos ? "--qos=${params.slurm_qos}" : null', slurm_text)
        self.assertIn("params.slurm_cluster_options ?: null", slurm_text)
        self.assertNotIn("errorStrategy", slurm_text)

    def test_local_profile_caps_memory_for_workstation_runs(self) -> None:
        """Keep local acceptance runs within a typical workstation memory budget."""
        config_text = LOCAL_CONFIG.read_text(encoding="utf-8")

        self.assertIn("params.max_memory = 16.GB", config_text)
        self.assertNotIn("withName: EGGNOG", config_text)

    def test_reporting_assets_overwrite_on_resumed_runs(self) -> None:
        """Allow resumed local runs to rewrite pipeline info artefacts safely."""
        config_text = BASE_CONFIG.read_text(encoding="utf-8")

        self.assertIn("timeline {\n    enabled = true\n    file = \"${params.outdir}/pipeline_info/timeline.html\"\n    overwrite = true", config_text)
        self.assertIn("report {\n    enabled = true\n    file = \"${params.outdir}/pipeline_info/report.html\"\n    overwrite = true", config_text)
        self.assertIn("trace {\n    enabled = true\n    file = \"${params.outdir}/pipeline_info/trace.tsv\"\n    overwrite = true", config_text)
        self.assertIn("dag {\n    enabled = true\n    file = \"${params.outdir}/pipeline_info/dag.html\"\n    overwrite = true", config_text)

    def test_prepare_workflow_registers_shared_process_names_without_global_suppression(self) -> None:
        """Keep prep-only selector handling out of the main workflow entry point."""
        prep_text = PREP_WORKFLOW.read_text(encoding="utf-8")
        main_text = MAIN_WORKFLOW.read_text(encoding="utf-8")

        self.assertIn(
            "include { INPUT_VALIDATION_AND_STAGING as UNUSED_INPUT_VALIDATION_AND_STAGING }",
            prep_text,
        )
        self.assertIn(
            "include { PER_SAMPLE_ANNOTATION as UNUSED_PER_SAMPLE_ANNOTATION }",
            prep_text,
        )
        self.assertNotIn("nextflow.enable.configProcessNamesValidation = false", prep_text)
        self.assertNotIn("nextflow.enable.configProcessNamesValidation = false", main_text)


if __name__ == "__main__":
    unittest.main()
