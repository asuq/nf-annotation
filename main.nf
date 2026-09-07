#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BUSCO_DATASET_PREP } from './subworkflows/local/busco_dataset_prep'
include { COHORT_16S } from './subworkflows/local/cohort_16s'
include { COHORT_ANI } from './subworkflows/local/cohort_ani'
include { COHORT_TAXONOMY } from './subworkflows/local/cohort_taxonomy'
include { FINAL_OUTPUTS } from './subworkflows/local/final_outputs'
include { INPUT_VALIDATION_AND_STAGING } from './subworkflows/local/input_validation_and_staging'
include { PER_SAMPLE_ANNOTATION } from './subworkflows/local/per_sample_annotation'
include { PER_SAMPLE_QC } from './subworkflows/local/per_sample_qc'

workflow {
    def normaliseBuscoLineages = { rawValue ->
        def rawItems = rawValue instanceof List ? rawValue : [rawValue]
        def lineages = rawItems
            .findAll { it != null }
            .collectMany { it.toString().split(',') as List }
            .collect { it.trim() }
            .findAll { it }
        if (lineages.isEmpty()) {
            error "params.busco_lineages must resolve to one or more lineage names."
        }
        def seen = [] as Set
        def duplicates = lineages.findAll { !seen.add(it) }.unique()
        if (!duplicates.isEmpty()) {
            error "params.busco_lineages must resolve to unique lineage names: ${duplicates.join(', ')}"
        }
        return lineages
    }
    def normaliseBooleanParam = { rawValue, paramName ->
        if (rawValue instanceof Boolean) {
            return rawValue
        }
        if (rawValue == null) {
            return false
        }
        def token = rawValue.toString().trim().toLowerCase()
        if (token in ['true', 't', 'yes', 'y', '1']) {
            return true
        }
        if (token in ['false', 'f', 'no', 'n', '0']) {
            return false
        }
        error "params.${paramName} must be boolean-like: true/false, yes/no, or 1/0."
    }

    if (!params.sample_csv) {
        error "params.sample_csv is required."
    }
    if (!params.metadata) {
        error "params.metadata is required."
    }
    if (!params.taxdump) {
        error "params.taxdump is required."
    }
    if (!params.update_from) {
        if (!params.checkm2_db) { error "params.checkm2_db is required." }
        if (!params.codetta_db) { error "params.codetta_db is required." }
        if (!params.eggnog_db) { error "params.eggnog_db is required." }
    }
    if (params.update_from) {
        def sourceRoot = new File(params.update_from.toString()).canonicalFile.toPath()
        def outputRoot = new File(params.outdir.toString()).canonicalFile.toPath()
        if (!sourceRoot.toFile().isDirectory()) {
            error "params.update_from must point to a published results directory."
        }
        if (outputRoot.startsWith(sourceRoot) || sourceRoot.startsWith(outputRoot)) {
            error "--update_from and --outdir must be separate, non-overlapping directories."
        }
        def existingOutputs = outputRoot.toFile().listFiles()?.findAll { it.name != 'pipeline_info' }
        if (existingOutputs) {
            if (!workflow.resume || !outputRoot.resolve('tables/cohort_update_run.json').toFile().isFile()) {
                error "A cohort update requires a fresh --outdir. Use -resume only for an interrupted update."
            }
        }
    }
    buscoLineagesList = normaliseBuscoLineages.call(params.busco_lineages)
    primaryBuscoColumn = (params.busco_primary_column ?: "BUSCO_${buscoLineagesList[0]}").toString()
    aniAllowIncomplete16s = normaliseBooleanParam.call(
        params.ani_allow_incomplete_16s,
        'ani_allow_incomplete_16s',
    )
    def aniThresholdValue = null
    try {
        aniThresholdValue = params.ani_threshold.toString().toDouble()
    } catch (NumberFormatException ignored) {
        error "params.ani_threshold must be a finite fraction in (0,1)."
    }
    if (!Double.isFinite(aniThresholdValue) || aniThresholdValue <= 0 || aniThresholdValue >= 1) {
        error "params.ani_threshold must be a finite fraction in (0,1)."
    }
    if (!(params.ani_score_profile in ['default', 'isolate', 'mag'])) {
        error "params.ani_score_profile must be default, isolate, or mag."
    }
    if (!(primaryBuscoColumn in buscoLineagesList.collect { lineage -> "BUSCO_${lineage}".toString() })) {
        error "params.busco_primary_column must identify one of the configured BUSCO lineages."
    }

    log.warn 'PADLOC and eggNOG outputs are retained in sample folders but are intentionally excluded from master_table.tsv.'

    sampleCsv = Channel.fromPath(params.sample_csv, checkIfExists: true)
    metadata = Channel.value(file(params.metadata, checkIfExists: true))
    taxdump = Channel.fromPath(params.taxdump, checkIfExists: true)

    INPUT_VALIDATION_AND_STAGING(
        sampleCsv,
        metadata,
        Channel.value(buscoLineagesList),
    )
    // Database paths are required only after preflight identifies an addition.
    hasNewSamples = INPUT_VALIDATION_AND_STAGING.out.new_staged_genomes
        .first()
        .map { item -> true }
    checkm2Db = hasNewSamples.map { present ->
        if (!params.checkm2_db) { error "params.checkm2_db is required for added samples." }
        file(params.checkm2_db, checkIfExists: true)
    }
    codettaDb = hasNewSamples.map { present ->
        if (!params.codetta_db) { error "params.codetta_db is required for added samples." }
        file(params.codetta_db, checkIfExists: true)
    }
    eggnogDb = hasNewSamples.map { present ->
        if (!params.eggnog_db) { error "params.eggnog_db is required for added samples." }
        file(params.eggnog_db, checkIfExists: true)
    }
    buscoLineages = hasNewSamples.flatMap { present -> buscoLineagesList }
    BUSCO_DATASET_PREP(buscoLineages)
    COHORT_TAXONOMY(INPUT_VALIDATION_AND_STAGING.out.validated_samples, metadata, taxdump)
    PER_SAMPLE_QC(
        INPUT_VALIDATION_AND_STAGING.out.new_staged_genomes,
        checkm2Db,
        BUSCO_DATASET_PREP.out.datasets,
    )
    allGcodeQc = PER_SAMPLE_QC.out.gcode_qc.mix(INPUT_VALIDATION_AND_STAGING.out.reused_gcode_qc)
    allSixteenS = PER_SAMPLE_QC.out.sixteen_s_summaries.mix(INPUT_VALIDATION_AND_STAGING.out.reused_sixteen_s)
    allBusco = PER_SAMPLE_QC.out.busco_summaries.mix(INPUT_VALIDATION_AND_STAGING.out.reused_busco_summaries)
    COHORT_16S(
        allSixteenS,
        PER_SAMPLE_QC.out.gcode_qc_for_cohort_16s.mix(INPUT_VALIDATION_AND_STAGING.out.reused_gcode_qc_for_cohort_16s),
        metadata,
    )
    PER_SAMPLE_ANNOTATION(
        INPUT_VALIDATION_AND_STAGING.out.new_staged_genomes,
        PER_SAMPLE_QC.out.gcode_qc,
        codettaDb,
        eggnogDb,
    )
    COHORT_ANI(
        INPUT_VALIDATION_AND_STAGING.out.validated_samples,
        metadata,
        INPUT_VALIDATION_AND_STAGING.out.staged_genomes,
        allGcodeQc,
        allSixteenS,
        allBusco,
        Channel.value(primaryBuscoColumn),
        Channel.value(aniAllowIncomplete16s),
    )
    FINAL_OUTPUTS(
        INPUT_VALIDATION_AND_STAGING.out.validated_samples,
        INPUT_VALIDATION_AND_STAGING.out.sample_status,
        metadata,
        COHORT_TAXONOMY.out.taxonomy,
        allGcodeQc,
        allSixteenS,
        COHORT_ANI.out.parsed_busco,
        PER_SAMPLE_ANNOTATION.out.codetta_summary.mix(INPUT_VALIDATION_AND_STAGING.out.reused_codetta_summary),
        PER_SAMPLE_ANNOTATION.out.ccfinder_summary.mix(INPUT_VALIDATION_AND_STAGING.out.reused_ccfinder_summary),
        PER_SAMPLE_ANNOTATION.out.prokka.mix(INPUT_VALIDATION_AND_STAGING.out.reused_prokka_results),
        PER_SAMPLE_ANNOTATION.out.padloc.mix(INPUT_VALIDATION_AND_STAGING.out.reused_padloc_results),
        PER_SAMPLE_ANNOTATION.out.eggnog.mix(INPUT_VALIDATION_AND_STAGING.out.reused_eggnog_results),
        PER_SAMPLE_ANNOTATION.out.eggnog_skips.mix(INPUT_VALIDATION_AND_STAGING.out.reused_eggnog_skip_rows),
        COHORT_ANI.out.clusters,
        COHORT_ANI.out.ani_metadata,
        COHORT_ANI.out.assembly_stats,
        COHORT_ANI.out.fastani_matrix,
        Channel.value(buscoLineagesList),
        Channel.value(primaryBuscoColumn),
        Channel.value(aniAllowIncomplete16s),
        INPUT_VALIDATION_AND_STAGING.out.versions
            .mix(BUSCO_DATASET_PREP.out.versions)
            .mix(COHORT_TAXONOMY.out.versions)
            .mix(PER_SAMPLE_QC.out.versions)
            .mix(COHORT_16S.out.versions)
            .mix(PER_SAMPLE_ANNOTATION.out.versions)
            .mix(COHORT_ANI.out.versions),
        INPUT_VALIDATION_AND_STAGING.out.inherited_versions,
        workflow.nextflow.version.toString(),
        workflow.manifest.version ?: 'NA',
        workflow.commitId ?: 'NA',
        workflow.containerEngine ?: 'none',
    )
}
