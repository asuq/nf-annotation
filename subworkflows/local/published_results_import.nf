include { COMPLETE_ANI_REFERENCES } from '../../modules/local/complete_ani_references'
include { FINGERPRINT_REUSE_SOURCES; VALIDATE_REUSE_SOURCES } from '../../modules/local/reuse_source_validation'
include { CHECK_COHORT_UPDATE } from '../../modules/local/check_cohort_update'
include { IMPORT_ANI_BATCH } from '../../modules/local/import_ani_batch'
include { PREPARE_COHORT_UPDATE } from '../../modules/local/prepare_cohort_update'
include { IMPORT_PUBLISHED_SAMPLE } from '../../modules/local/import_published_sample'
include { PREPARE_ANNOTATION_BUNDLE } from '../../modules/local/prepare_annotation_bundle'

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
    qcRoot = Channel.value(params.ani_publication_mode == 'reference' ? file(params.qc_publication_root) : [])
    sourceRoots = PublishedSources.paths(params.update_from)
    if (params.ani_publication_mode == 'reference') { ReferenceInputs.requireQcProviders(sourceRoots, params.qc_publication_root) }
    dependencyRoots = ReferenceInputs.paths(sourceRoots).findAll { !(it in sourceRoots) }
    candidateGenomes = validated_samples
        .splitCsv(header: true, sep: '\t')
        .map { row -> java.nio.file.Path.of(row.genome_fasta).toRealPath() }
        .collect()
    previousIdentity = file("${params.outdir}/tables/cohort_update_run.json")
    previousUpdate = Channel.value(previousIdentity.exists() ? [previousIdentity] : [])
    settingsNames = [
        'domain', 'busco_domain_lineages', 'pipeline_stage',
        'busco_lineages', 'busco_primary_column', 'gcode_rule',
        'codetta_extra_args', 'ccfinder_extra_args',
        'ani_16s_policy', 'ani_threshold', 'ani_score_profile',
        'taxdump', 'taxdump_label', 'checkm2_db', 'checkm2_db_label', 'codetta_db', 'codetta_db_label',
        'busco_db', 'prepare_busco_datasets',
        'python_container', 'seqtk_container', 'barrnap_container', 'checkm2_container', 'busco_container',
        'prokka_container', 'codetta_container', 'ccfinder_container', 'fastani_container',
    ]
    if (SamplePolicy.stage(params) == 'ani') { settingsNames.addAll(['ani_import_batches', 'ani_publication_mode']) }
    updateSettings = settingsNames.collectEntries { name -> [(name): params[name]] }
    if (params.ani_publication_mode == 'reference') {
        updateSettings.qc_publication_root = file(params.qc_publication_root).toFile().canonicalPath
    }

    FINGERPRINT_REUSE_SOURCES(validated_samples.map { sourceRoots }, Channel.value(dependencyRoots))
    sourceFingerprints = FINGERPRINT_REUSE_SOURCES.out.fingerprints.map { path ->
        new groovy.json.JsonSlurper().parse(path.toFile()).sort(false) { it.root }
    }.first()
    sourceGroups = sourceFingerprints.flatMap { records ->
        def size = Math.max(1, Math.ceil(records.size() / 8.0) as int)
        records.collate(size).collect { group ->
            def roots = group.collect { file(it.root) }
            tuple(group, roots, ReferenceInputs.paths(roots).findAll { !(it in roots) })
        }
    }
    VALIDATE_REUSE_SOURCES(sourceGroups, Channel.value(SamplePolicy.reuseCodeHash(projectDir)))
    sourceReceipts = VALIDATE_REUSE_SOURCES.out.receipts.toList().map { paths ->
        paths.collectMany { new groovy.json.JsonSlurper().parse(it.toFile()) }.sort(false) { it.source.root }
    }
    checkedSources = sourceFingerprints.map { [it] }.combine(sourceReceipts.map { [it] }).map { item ->
        // Keep nested lists intact when combining value channels.
        [fingerprints: item[0], receipts: item[1]]
    }.first()
    PREPARE_COHORT_UPDATE(
        validated_samples, accession_map, initial_status, validation_warnings,
        Channel.value(sourceRoots), candidateGenomes, metadata,
        busco_lineages, Channel.value(updateSettings),
        Channel.value(file(params.outdir).toAbsolutePath().normalize().toString()),
        Channel.value(SamplePolicy.codeHash(projectDir)),
        checkedSources,
        Channel.value(dependencyRoots),
    )

    CHECK_COHORT_UPDATE(PREPARE_COHORT_UPDATE.out.prepared, previousUpdate)

    newGenomes = CHECK_COHORT_UPDATE.out.new_samples
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row, java.nio.file.Path.of(row.genome_fasta).toRealPath()) }
    reusedSamples = CHECK_COHORT_UPDATE.out.reused_samples
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = row.findAll { key, value -> !(key in ['source_gcode', 'source_index']) }
            tuple(meta, sourceRoots[row.source_index as int].resolve("samples/${row.accession}"), row.source_gcode)
        }
    if (SamplePolicy.stage(params) == 'ani') {
        destination = file(params.outdir).toFile().canonicalFile.toPath()
        destination.toFile().mkdirs()
        batches = reusedSamples.toList().flatMap { samples ->
            if (!samples) { return [] }
            def size = Math.max(1, Math.ceil(samples.size() / (params.ani_import_batches as double)) as int)
            samples.collate(size).withIndex().collect { group, index ->
                tuple(index, group.collect { [meta: it[0], gcode: it[2].toString(), source_index: sourceRoots.indexOf(it[1].parent.parent)] })
            }
        }
        IMPORT_ANI_BATCH(batches, Channel.value(sourceRoots), Channel.value(dependencyRoots), Channel.value(destination), Channel.value(destination.toString()), checkedSources, qcRoot)
        if (params.ani_publication_mode == 'reference') {
            COMPLETE_ANI_REFERENCES(IMPORT_ANI_BATCH.out.audit.toList(), CHECK_COHORT_UPDATE.out.validated_samples,
                Channel.value(destination), Channel.value(ReferenceInputs.paths(sourceRoots)), qcRoot)
            readyImports = IMPORT_ANI_BATCH.out.results.combine(COMPLETE_ANI_REFERENCES.out.complete).map { item -> item.take(3) }
        } else {
            readyImports = IMPORT_ANI_BATCH.out.results
        }
        uncheckedImports = readyImports.flatMap { item ->
            if (!(item instanceof List) || item.size() != 3) { error 'Invalid ANI import batch result tuple' }
            def records = new groovy.json.JsonSlurper().parse(item[0].toFile()).samples
            def channels = item[1]
            def root = item[2]
            records.collect { record ->
                tuple(record.meta, java.nio.file.Paths.get(record.native_root), channels, record.gcode, [], record.bundle_root ?: '')
            }
        }
        expectedImports = CHECK_COHORT_UPDATE.out.reused_samples.splitCsv(header: true, sep: '\t').toList()
        imported = uncheckedImports.toList().map { [it] }.combine(expectedImports.map { [it] }).flatMap { item ->
            def observed = item[0]
            def expected = item[1].collect { it.accession }
            def byAccession = observed.groupBy { it[0].accession }
            if (byAccession.keySet() != (expected as Set) || byAccession.values().any { it.size() != 1 }) {
                error 'ANI imports did not produce exactly one result for every requested sample'
            }
            expected.collect { byAccession[it][0] }
        }
        reusedBundles = imported.filter { it[3] in ['4', '11'] }.map { item ->
            tuple(item[0], java.nio.file.Paths.get(item[5]))
        }
        importVersions = IMPORT_ANI_BATCH.out.versions
        reusedBuscoJson = Channel.empty()
        reusedBuscoTables = imported.combine(busco_lineages.map { values -> [values] }).flatMap { item ->
            SamplePolicy.sampleLineages(item[0], item[6] as List).collect { lineage ->
                tuple(item[0], lineage, item[1].resolve("busco/${lineage}/busco_summary_${lineage}.tsv"))
            }
        }
    } else {
        IMPORT_PUBLISHED_SAMPLE(reusedSamples.map { item ->
            def root = item[1].parent.parent
            tuple(item[0], item[1], item[2], ReferenceInputs.paths([root]).findAll { it != root })
        }, busco_lineages)
        imported = IMPORT_PUBLISHED_SAMPLE.out.results.map { item ->
            tuple(item[0], item[1].first().parent, item[2], item[3], item[4])
        }
        // Revalidate retained native inputs with the current bundle adapter. Search
        // reuse remains keyed to protein/coordinate identities, not producer code.
        bundleInputs = imported.filter { item -> item[3] in ['4', '11'] }.map { item ->
            tuple(
                item[0], item[1].resolve("staged/${item[0].internal_id}.fasta"), item[3],
                item[1].resolve('prokka/prokka.faa'), item[1].resolve('prokka/prokka.gff'),
                item[1].resolve('prokka/prokka.gbk'), item[1].resolve('prokka/prokka.log'),
            )
        }
        PREPARE_ANNOTATION_BUNDLE(bundleInputs)
        reusedBundles = PREPARE_ANNOTATION_BUNDLE.out.bundle
        importVersions = IMPORT_PUBLISHED_SAMPLE.out.versions.mix(PREPARE_ANNOTATION_BUNDLE.out.versions)
        reusedBuscoJson = imported.flatMap { item -> SamplePolicy.sampleLineages(item[0], item[4] as List).collect { lineage -> tuple(item[0], lineage, item[1].resolve("busco/${lineage}/short_summary.json")) } }
        reusedBuscoTables = Channel.empty()
    }
    annotated = imported.filter { item -> item[3] in ['4', '11'] }

    emit:
    validated_samples = CHECK_COHORT_UPDATE.out.validated_samples
    accession_map = CHECK_COHORT_UPDATE.out.accession_map
    validation_warnings = CHECK_COHORT_UPDATE.out.validation_warnings
    sample_status = CHECK_COHORT_UPDATE.out.sample_status
    new_genomes = newGenomes
    staged_genomes = imported.map { item -> tuple(item[0], item[1].resolve("staged/${item[0].internal_id}.fasta")) }
    gcode_qc = imported.map { item -> tuple(item[0], item[1].resolve('checkm2/checkm2_summary.tsv')) }
    gcode_qc_for_cohort_16s = imported.map { item -> tuple(item[0], item[2].resolve("${item[0].internal_id}_checkm2_summary.tsv")) }
    sixteen_s_summaries = imported.map { item -> tuple(item[0], item[2].resolve("${item[0].internal_id}_best_16S.fna"), item[2].resolve("${item[0].internal_id}_16S_status.tsv")) }
    busco_summaries = reusedBuscoJson
    parsed_busco = reusedBuscoTables
    codetta_summary = imported.map { item -> tuple(item[0], item[1].resolve('codetta/codetta_summary.tsv')) }
    ccfinder_summary = annotated.filter { item -> (item[0].qc_scope ?: 'full') == 'full' }.map { item -> tuple(item[0], item[1].resolve('ccfinder/ccfinder_strains.tsv'), item[1].resolve('ccfinder/ccfinder_contigs.tsv'), item[1].resolve('ccfinder/ccfinder_crisprs.tsv')) }
    prokka = annotated.map { item -> tuple(item[0], item[1].resolve('prokka'), item[1].resolve('prokka/prokka.gff'), item[1].resolve('prokka/prokka.faa'), item[1].resolve('prokka/prokka.gbk'), item[1].resolve('prokka/prokka.log')) }
    bundles = reusedBundles
    inherited_versions = CHECK_COHORT_UPDATE.out.inherited_versions.map { report -> [report] }
    versions = PREPARE_COHORT_UPDATE.out.versions.mix(importVersions)
}
