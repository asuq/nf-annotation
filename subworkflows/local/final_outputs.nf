include { REPORT_GENOMES } from '../../modules/local/report_genomes'
include { REFERENCE_ACCEPTANCE } from '../../modules/local/reference_acceptance'
include { BUILD_MASTER_TABLE } from '../../modules/local/build_master_table'
include { COLLECT_VERSIONS } from '../../modules/local/collect_versions'
include { WRITE_SAMPLE_STATUS } from '../../modules/local/write_sample_status'
include { AGGREGATE_ANNOTATIONS; ANNOTATION_ACCEPTANCE } from '../../modules/local/functional_annotation'

/*
 * Build the final cohort tables from the keyed per-sample summaries and gather
 * one authoritative provenance report.
 */
workflow FINAL_OUTPUTS {
    take:
    validated_samples
    initial_status
    metadata
    taxonomy
    checkm2_summaries
    sixteen_s_summaries
    busco_summaries
    codetta_summaries
    ccfinder_summaries
    prokka_results
    annotation_results
    annotation_plan
    annotation_bundles
    annotation_batches
    cohort_evidence
    aaiSummary
    aniSummary
    aniRepresentatives
    ani_metadata
    assembly_stats
    ani_matrix
    busco_lineages
    primary_busco_column
    ani_16s_policy
    version_files
    inherited_versions
    nextflow_version
    pipeline_version
    git_commit
    container_engine

    main:
    checkm2_seed = Channel.value(file("${projectDir}/assets/tables/headers/checkm2_summary.tsv"))
    sixteen_s_seed = Channel.value(file("${projectDir}/assets/tables/headers/16s_status.tsv"))
    codetta_seed = Channel.value(file("${projectDir}/assets/tables/headers/codetta_summary.tsv"))
    ccfinder_seed = Channel.value(file("${projectDir}/assets/tables/headers/ccfinder_strains.tsv"))
    finalOutputsCollectDir = file("${workflow.workDir}/collect/${workflow.sessionId}/final_outputs")

    extractExitCode = { logFile ->
        def file = logFile.toFile()
        if (!file.exists()) {
            return 'NA'
        }
        def matches = file.readLines().findAll { line -> line.startsWith('exit_code=') }
        return matches ? matches[-1].split('=', 2)[1].trim() : 'NA'
    }

    unpackTuple = { item, channelName, expectedSize ->
        if (!(item instanceof List)) {
            def actualType = item == null ? 'null' : item.getClass().getName()
            throw new IllegalArgumentException(
                "${channelName} expected a ${expectedSize}-value tuple, received ${actualType}."
            )
        }
        if (item.size() != expectedSize) {
            def preview = item.take(Math.min(item.size(), 3))
            throw new IllegalArgumentException(
                "${channelName} expected a ${expectedSize}-value tuple, " +
                    "received ${item.size()} value(s): ${preview}"
            )
        }
        return item
    }

    combined_checkm2 = checkm2_seed
        .mix(checkm2_summaries.map { item ->
            def values = unpackTuple.call(item, 'checkm2_summaries', 2)
            values[1]
        })
        .collectFile(
            name: 'checkm2_summaries.tsv',
            keepHeader: true,
            skip: 1,
            newLine: true,
        )

    combined_16s = sixteen_s_seed
        .mix(sixteen_s_summaries.map { item ->
            def values = unpackTuple.call(item, 'sixteen_s_summaries', 3)
            values[2]
        })
        .collectFile(
            name: '16s_statuses.tsv',
            keepHeader: true,
            skip: 1,
            newLine: true,
        )

    combined_ccfinder = ccfinder_seed
        .mix(ccfinder_summaries.map { item ->
            def values = unpackTuple.call(item, 'ccfinder_summaries', 4)
            values[1]
        })
        .collectFile(
            name: 'ccfinder_strains.tsv',
            keepHeader: true,
            skip: 1,
            newLine: true,
        )

    combined_codetta = codetta_seed
        .mix(codetta_summaries.map { item ->
            def values = unpackTuple.call(item, 'codetta_summaries', 2)
            values[1]
        })
        .collectFile(
            name: 'codetta_summary.tsv',
            keepHeader: true,
            skip: 1,
            newLine: true,
        )

    collected_busco = busco_summaries

    prokkaManifest = Channel
        .of('accession\texit_code\tgff_size\tfaa_size')
        .concat(
            prokka_results.map { item ->
                def values = unpackTuple.call(item, 'prokka_results', 6)
                def meta = values[0]
                def gff = values[2]
                def faa = values[3]
                def log = values[5]
                "${meta.accession}\t${extractExitCode.call(log)}\t${gff.toFile().length()}\t${faa.toFile().length()}"
            }
        )
        .collectFile(
            name: 'prokka_manifest.tsv',
            newLine: true,
            sort: false,
            storeDir: finalOutputsCollectDir,
        )

    BUILD_MASTER_TABLE(
        validated_samples,
        metadata,
        busco_lineages,
        taxonomy,
        combined_checkm2,
        combined_16s,
        collected_busco,
        combined_codetta,
        combined_ccfinder,
        aniSummary,
        aaiSummary,
        assembly_stats,
        Channel.value(SamplePolicy.codeHash(projectDir)),
    )

    WRITE_SAMPLE_STATUS(
        validated_samples,
        initial_status,
        busco_lineages,
        metadata,
        taxonomy,
        combined_checkm2,
        combined_16s,
        collected_busco,
        combined_codetta,
        combined_ccfinder,
        prokkaManifest,
        aniSummary,
        aaiSummary,
        assembly_stats,
        primary_busco_column,
        ani_16s_policy,
        Channel.value(SamplePolicy.codeHash(projectDir)),
    )

    def isAni = SamplePolicy.stage(params) == 'clustering'
    publicationDependencies = isAni && params.publication_mode == 'reference' ? Channel.value(ReferenceInputs.paths(PublishedSources.paths(params.update_from))) : Channel.value([])
    if (!SamplePolicy.annotationTools(params)) {
        manifests = annotation_bundles.map { items -> items.collect { it.resolve('bundle.json') } }
        references = annotation_bundles.map { params.publication_mode == 'reference' ? file("${params.outdir}/cohort_references.json") : [] }
        REPORT_GENOMES(annotation_plan, BUILD_MASTER_TABLE.out.master_table, WRITE_SAMPLE_STATUS.out.sample_status, manifests, references, cohort_evidence)
        reportTables = REPORT_GENOMES.out.tables
        reportManifest = REPORT_GENOMES.out.manifest
        reportAcceptance = REPORT_GENOMES.out.acceptance
    } else {
        AGGREGATE_ANNOTATIONS(annotation_plan, BUILD_MASTER_TABLE.out.master_table, WRITE_SAMPLE_STATUS.out.sample_status,
            annotation_bundles, annotation_results.toList().map { rows -> tuple(rows.collect { [accession: it[0].accession, tool: it[0].tool] }, rows.collect { it[1] }) }, annotation_batches.toList(),
            Channel.value([]), Channel.value([]), cohort_evidence)
        reportTables = AGGREGATE_ANNOTATIONS.out.tables
        reportManifest = AGGREGATE_ANNOTATIONS.out.manifest
        reportAcceptance = AGGREGATE_ANNOTATIONS.out.acceptance
    }

    final_versions = BUILD_MASTER_TABLE.out.versions
        .mix(WRITE_SAMPLE_STATUS.out.versions)

    collected_versions = version_files
        .mix(final_versions)
        .collect()
        .map { files ->
            files.toList().sort { versionFile -> versionFile.toString() }
        }

    COLLECT_VERSIONS(
        collected_versions,
        inherited_versions,
        busco_lineages,
        nextflow_version,
        pipeline_version,
        git_commit,
        container_engine,
    )

    if (isAni && params.publication_mode == 'reference') {
        REFERENCE_ACCEPTANCE(reportAcceptance, COLLECT_VERSIONS.out.versions_table,
            Channel.value(file(params.outdir)), publicationDependencies, Channel.value([]))
    } else {
        ANNOTATION_ACCEPTANCE(reportAcceptance, COLLECT_VERSIONS.out.versions_table)
    }

    emit:
    master_table = reportTables.map { files -> files.find { it.name == 'master_table.tsv' } }
    sample_status = reportTables.map { files -> files.find { it.name == 'sample_status.tsv' } }
    annotation_manifest = reportManifest
    ani_representatives = aniRepresentatives
    versions_table = COLLECT_VERSIONS.out.versions_table
    versions = final_versions
}
