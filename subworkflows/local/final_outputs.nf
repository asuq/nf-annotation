include { BUILD_MASTER_TABLE } from '../../modules/local/build_master_table'
include { COLLECT_VERSIONS } from '../../modules/local/collect_versions'
include { SELECT_ANI_REPRESENTATIVES } from '../../modules/local/select_ani_representatives'
include { WRITE_SAMPLE_STATUS } from '../../modules/local/write_sample_status'

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
    padloc_results
    eggnog_results
    eggnog_skips
    ani_clusters
    ani_metadata
    assembly_stats
    ani_matrix
    busco_lineages
    primary_busco_column
    ani_allow_incomplete_16s
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

    countTopLevelFiles = { directoryPath ->
        def files = directoryPath.toFile().listFiles()
        return files == null ? 0 : files.count { file -> file.isFile() }
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

    padlocManifest = Channel
        .of('accession\texit_code\tresult_file_count')
        .concat(
            padloc_results.map { item ->
                def values = unpackTuple.call(item, 'padloc_results', 3)
                def meta = values[0]
                def padlocDir = values[1]
                def log = values[2]
                "${meta.accession}\t${extractExitCode.call(log)}\t${countTopLevelFiles.call(padlocDir)}"
            }
        )
        .collectFile(
            name: 'padloc_manifest.tsv',
            newLine: true,
            sort: false,
            storeDir: finalOutputsCollectDir,
        )

    eggnogManifest = Channel
        .of('accession\tstatus\twarnings\texit_code\tannotations_size\tresult_file_count')
        .concat(
            eggnog_results.map { item ->
                def values = unpackTuple.call(item, 'eggnog_results', 4)
                def meta = values[0]
                def eggnogDir = values[1]
                def annotations = values[2]
                def log = values[3]
                def exitCode = extractExitCode.call(log)
                def annotationsSize = annotations.toFile().length()
                def resultFileCount = countTopLevelFiles.call(eggnogDir)
                def status = (exitCode == '0' && annotationsSize > 0) ? 'done' : 'failed'
                def warnings = ''
                if (status == 'failed') {
                    warnings = exitCode == '0' ? 'missing_eggnog_outputs' : 'eggnog_failed'
                }
                "${meta.accession}\t${status}\t${warnings}\t${exitCode}\t${annotationsSize}\t${resultFileCount}"
            }
        )
        .concat(eggnog_skips)
        .collectFile(
            name: 'eggnog_manifest.tsv',
            newLine: true,
            sort: false,
            storeDir: finalOutputsCollectDir,
        )

    SELECT_ANI_REPRESENTATIVES(
        ani_clusters,
        ani_metadata,
        ani_matrix,
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
        SELECT_ANI_REPRESENTATIVES.out.ani_summary,
        assembly_stats,
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
        padlocManifest,
        eggnogManifest,
        SELECT_ANI_REPRESENTATIVES.out.ani_summary,
        assembly_stats,
        primary_busco_column,
        ani_allow_incomplete_16s,
    )

    final_versions = SELECT_ANI_REPRESENTATIVES.out.versions
        .mix(BUILD_MASTER_TABLE.out.versions)
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

    emit:
    master_table = BUILD_MASTER_TABLE.out.master_table
    sample_status = WRITE_SAMPLE_STATUS.out.sample_status
    ani_representatives = SELECT_ANI_REPRESENTATIVES.out.ani_representatives
    versions_table = COLLECT_VERSIONS.out.versions_table
    versions = final_versions
}
