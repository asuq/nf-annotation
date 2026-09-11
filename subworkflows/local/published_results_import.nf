include { PREPARE_COHORT_UPDATE } from '../../modules/local/prepare_cohort_update'
include { IMPORT_PUBLISHED_SAMPLE } from '../../modules/local/import_published_sample'

/* Adapt published samples to the existing per-sample output channel contracts. */
workflow PUBLISHED_RESULTS_IMPORT {
    take:
    validated_samples
    accession_map
    initial_status
    validation_warnings
    metadata
    busco_lineages

    main:
    sourceRoot = file(params.update_from, checkIfExists: true)
    candidateGenomes = validated_samples
        .splitCsv(header: true, sep: '\t')
        .map { row -> java.nio.file.Path.of(row.genome_fasta).toRealPath() }
        .collect()
    previousIdentity = file("${params.outdir}/tables/cohort_update_run.json")
    previousUpdate = Channel.value(previousIdentity.exists() ? [previousIdentity] : [])
    settingsNames = [
        'barrnap_kingdom', 'busco_lineages', 'busco_primary_column', 'gcode_rule',
        'codetta_extra_args', 'ccfinder_extra_args', 'padloc_extra_args', 'eggnog_extra_args',
        'eggnog_only_accessions', 'ani_allow_incomplete_16s', 'ani_threshold', 'ani_score_profile',
        'taxdump', 'taxdump_label', 'checkm2_db', 'checkm2_db_label', 'codetta_db', 'codetta_db_label',
        'busco_db', 'eggnog_db', 'eggnog_db_label', 'prepare_busco_datasets',
        'python_container', 'seqtk_container', 'barrnap_container', 'checkm2_container', 'busco_container',
        'prokka_container', 'codetta_container', 'ccfinder_container', 'padloc_container', 'eggnog_container', 'fastani_container',
    ]
    updateSettings = settingsNames.collectEntries { name -> [(name): params[name]] }

    PREPARE_COHORT_UPDATE(
        validated_samples, accession_map, initial_status, validation_warnings,
        Channel.value(sourceRoot), candidateGenomes, metadata, previousUpdate,
        busco_lineages, Channel.value(updateSettings),
        Channel.value(file(params.outdir).toAbsolutePath().normalize().toString()),
    )

    newGenomes = PREPARE_COHORT_UPDATE.out.new_samples
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row, java.nio.file.Path.of(row.genome_fasta).toRealPath()) }
    reusedSamples = PREPARE_COHORT_UPDATE.out.reused_samples
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = row.findAll { key, value -> !(key in ['source_gcode', 'source_eggnog_status']) }
            tuple(meta, sourceRoot.resolve("samples/${row.accession}"), row.source_gcode, row.source_eggnog_status)
        }
    IMPORT_PUBLISHED_SAMPLE(reusedSamples, busco_lineages)
    imported = IMPORT_PUBLISHED_SAMPLE.out.results
    annotated = imported.filter { item -> item[3] in ['4', '11'] }

    emit:
    validated_samples = PREPARE_COHORT_UPDATE.out.validated_samples
    accession_map = PREPARE_COHORT_UPDATE.out.accession_map
    validation_warnings = PREPARE_COHORT_UPDATE.out.validation_warnings
    sample_status = PREPARE_COHORT_UPDATE.out.sample_status
    new_genomes = newGenomes
    staged_genomes = imported.map { item -> tuple(item[0], item[1].resolve("staged/${item[0].internal_id}.fasta")) }
    gcode_qc = imported.map { item -> tuple(item[0], item[1].resolve('checkm2/checkm2_summary.tsv')) }
    gcode_qc_for_cohort_16s = imported.map { item -> tuple(item[0], item[2].resolve("${item[0].internal_id}_checkm2_summary.tsv")) }
    sixteen_s_summaries = imported.map { item -> tuple(item[0], item[2].resolve("${item[0].internal_id}_best_16S.fna"), item[2].resolve("${item[0].internal_id}_16S_status.tsv")) }
    busco_summaries = imported.flatMap { item -> item[5].collect { lineage -> tuple(item[0], lineage, item[1].resolve("busco/${lineage}/short_summary.json")) } }
    codetta_summary = imported.map { item -> tuple(item[0], item[1].resolve('codetta/codetta_summary.tsv')) }
    ccfinder_summary = annotated.map { item -> tuple(item[0], item[1].resolve('ccfinder/ccfinder_strains.tsv'), item[1].resolve('ccfinder/ccfinder_contigs.tsv'), item[1].resolve('ccfinder/ccfinder_crisprs.tsv')) }
    prokka = annotated.map { item -> tuple(item[0], item[1].resolve('prokka'), item[1].resolve('prokka/prokka.gff'), item[1].resolve('prokka/prokka.faa'), item[1].resolve('prokka/prokka.gbk'), item[1].resolve('prokka/prokka.log')) }
    padloc = annotated.map { item -> tuple(item[0], item[1].resolve('padloc/padloc'), item[1].resolve('padloc/padloc.log')) }
    eggnog = annotated.filter { item -> item[4] != 'skipped' }.map { item -> tuple(item[0], item[1].resolve('eggnog/eggnog'), item[1].resolve('eggnog/eggnog_annotations.tsv'), item[1].resolve('eggnog/eggnog.log')) }
    eggnog_skips = annotated.filter { item -> item[4] == 'skipped' }.map { item -> "${item[0].accession}\tskipped\teggnog_short_circuit\t0\t0\t0" }
    inherited_versions = PREPARE_COHORT_UPDATE.out.inherited_versions.map { report -> [report] }
    versions = PREPARE_COHORT_UPDATE.out.versions.mix(IMPORT_PUBLISHED_SAMPLE.out.versions)
}
