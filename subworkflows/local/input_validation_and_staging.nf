include { VALIDATE_INPUTS } from '../../modules/local/validate_inputs'
include { STAGE_INPUTS } from '../../modules/local/stage_inputs'
include { PUBLISHED_RESULTS_IMPORT } from './published_results_import'

/*
 * Validate top-level inputs once, then stage every requested genome to a
 * canonical internal-ID-safe FASTA for downstream per-sample tools.
 */
workflow INPUT_VALIDATION_AND_STAGING {
    take:
    sample_csv
    metadata
    busco_lineages

    main:
    VALIDATE_INPUTS(sample_csv, metadata, busco_lineages)

    if (params.update_from) {
        PUBLISHED_RESULTS_IMPORT(
            VALIDATE_INPUTS.out.validated_samples,
            VALIDATE_INPUTS.out.accession_map,
            VALIDATE_INPUTS.out.sample_status,
            VALIDATE_INPUTS.out.validation_warnings,
            metadata,
            busco_lineages,
        )
        sample_genomes = PUBLISHED_RESULTS_IMPORT.out.new_genomes
        validated = PUBLISHED_RESULTS_IMPORT.out.validated_samples
        accession_mapping = PUBLISHED_RESULTS_IMPORT.out.accession_map
        warnings = PUBLISHED_RESULTS_IMPORT.out.validation_warnings
        initial_status = PUBLISHED_RESULTS_IMPORT.out.sample_status
        reused_genomes = PUBLISHED_RESULTS_IMPORT.out.staged_genomes
        reused_qc = PUBLISHED_RESULTS_IMPORT.out.gcode_qc
        reused_cohort_qc = PUBLISHED_RESULTS_IMPORT.out.gcode_qc_for_cohort_16s
        reused_16s = PUBLISHED_RESULTS_IMPORT.out.sixteen_s_summaries
        reused_busco = PUBLISHED_RESULTS_IMPORT.out.busco_summaries
        reused_codetta = PUBLISHED_RESULTS_IMPORT.out.codetta_summary
        reused_ccfinder = PUBLISHED_RESULTS_IMPORT.out.ccfinder_summary
        reused_prokka = PUBLISHED_RESULTS_IMPORT.out.prokka
        reused_padloc = PUBLISHED_RESULTS_IMPORT.out.padloc
        reused_eggnog = PUBLISHED_RESULTS_IMPORT.out.eggnog
        reused_eggnog_skips = PUBLISHED_RESULTS_IMPORT.out.eggnog_skips
        inherited_version_reports = PUBLISHED_RESULTS_IMPORT.out.inherited_versions
        import_versions = PUBLISHED_RESULTS_IMPORT.out.versions
    } else {
        sample_genomes = VALIDATE_INPUTS.out.validated_samples
            .splitCsv(header: true, sep: '\t')
            .map { row ->
                def meta = row.collectEntries { key, value -> [(key): value] }
                def genomePath = java.nio.file.Path.of(row.genome_fasta).toRealPath()
                tuple(meta, genomePath)
            }
        validated = VALIDATE_INPUTS.out.validated_samples
        accession_mapping = VALIDATE_INPUTS.out.accession_map
        warnings = VALIDATE_INPUTS.out.validation_warnings
        initial_status = VALIDATE_INPUTS.out.sample_status
        reused_genomes = Channel.empty()
        reused_qc = Channel.empty()
        reused_cohort_qc = Channel.empty()
        reused_16s = Channel.empty()
        reused_busco = Channel.empty()
        reused_codetta = Channel.empty()
        reused_ccfinder = Channel.empty()
        reused_prokka = Channel.empty()
        reused_padloc = Channel.empty()
        reused_eggnog = Channel.empty()
        reused_eggnog_skips = Channel.empty()
        inherited_version_reports = Channel.value([])
        import_versions = Channel.empty()
    }

    STAGE_INPUTS(sample_genomes)

    versions = VALIDATE_INPUTS.out.versions.mix(STAGE_INPUTS.out.versions).mix(import_versions)

    emit:
    validated_samples = validated
    accession_map = accession_mapping
    validation_warnings = warnings
    sample_status = initial_status
    staged_genomes = STAGE_INPUTS.out.staged_fasta.mix(reused_genomes)
    new_staged_genomes = STAGE_INPUTS.out.staged_fasta
    reused_gcode_qc = reused_qc
    reused_gcode_qc_for_cohort_16s = reused_cohort_qc
    reused_sixteen_s = reused_16s
    reused_busco_summaries = reused_busco
    reused_codetta_summary = reused_codetta
    reused_ccfinder_summary = reused_ccfinder
    reused_prokka_results = reused_prokka
    reused_padloc_results = reused_padloc
    reused_eggnog_results = reused_eggnog
    reused_eggnog_skip_rows = reused_eggnog_skips
    inherited_versions = inherited_version_reports
    staged_fai = STAGE_INPUTS.out.fai
    versions = versions
}
