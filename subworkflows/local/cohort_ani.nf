include { SUMMARISE_BUSCO } from '../../modules/local/summarise_busco'
include { CALCULATE_ASSEMBLY_STATS } from '../../modules/local/calculate_assembly_stats'
include { BUILD_FASTANI_INPUTS } from '../../modules/local/build_fastani_inputs'
include { FASTANI } from '../../modules/local/fastani'
include { CLUSTER_ANI } from '../../modules/local/cluster_ani'

/*
 * Build the ANI-eligible cohort, run all-vs-all FastANI, and emit stable ANI
 * cluster memberships for downstream representative selection.
 */
workflow COHORT_ANI {
    take:
    validated_samples
    metadata
    staged_genomes
    gcode_qc
    sixteen_s_summaries
    busco_summaries
    primary_busco_column
    ani_16s_policy

    main:
    SUMMARISE_BUSCO(busco_summaries)

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

    staged_manifest = Channel
        .of('accession\tinternal_id\tstaged_filename')
        .concat(
            staged_genomes.map { item ->
                def values = unpackTuple.call(item, 'staged_genomes', 2)
                def meta = values[0]
                def genome = values[1]
                "${meta.accession}\t${meta.internal_id}\t${genome.getName()}"
            }
        )
        .collectFile(name: 'staged_genomes.tsv', newLine: true, sort: false)

    staged_fasta_files = staged_genomes
        .map { item ->
            def values = unpackTuple.call(item, 'staged_genomes', 2)
            values[1]
        }
        .collect()

    combined_checkm2 = gcode_qc
        .map { item ->
            def values = unpackTuple.call(item, 'gcode_qc', 2)
            values[1]
        }
        .collectFile(
            name: 'checkm2_summaries.tsv',
            keepHeader: true,
            skip: 1,
            newLine: true,
        )

    combined_16s = sixteen_s_summaries
        .map { item ->
            def values = unpackTuple.call(item, 'sixteen_s_summaries', 3)
            values[2]
        }
        .collectFile(
            name: '16s_statuses.tsv',
            keepHeader: true,
            skip: 1,
            newLine: true,
        )

    busco_tables = SUMMARISE_BUSCO.out.summary
        .map { item ->
            def values = unpackTuple.call(item, 'busco_summary', 3)
            def lineage = values[1]
            def summary = values[2]
            tuple("busco_summary_${lineage}.tsv", summary)
        }
        .collectFile(
            keepHeader: true,
            skip: 1,
        ) { item ->
            def values = unpackTuple.call(item, 'busco_table_file', 2)
            def outputName = values[0]
            def summary = values[1]
            [outputName, summary]
        }
        .collect()

    CALCULATE_ASSEMBLY_STATS(staged_manifest, staged_fasta_files)

    BUILD_FASTANI_INPUTS(
        validated_samples,
        metadata,
        staged_manifest,
        staged_fasta_files,
        combined_checkm2,
        combined_16s,
        busco_tables,
        CALCULATE_ASSEMBLY_STATS.out.stats,
        primary_busco_column,
        ani_16s_policy,
    )
    FASTANI(BUILD_FASTANI_INPUTS.out.fastani_inputs, BUILD_FASTANI_INPUTS.out.paths)
    CLUSTER_ANI(FASTANI.out.matrix, BUILD_FASTANI_INPUTS.out.metadata)

    versions = SUMMARISE_BUSCO.out.versions
        .mix(CALCULATE_ASSEMBLY_STATS.out.versions)
        .mix(BUILD_FASTANI_INPUTS.out.versions)
        .mix(FASTANI.out.versions)
        .mix(CLUSTER_ANI.out.versions)

    emit:
    parsed_busco = busco_tables
    fastani_inputs = BUILD_FASTANI_INPUTS.out.fastani_inputs
    fastani_paths = BUILD_FASTANI_INPUTS.out.paths
    ani_metadata = BUILD_FASTANI_INPUTS.out.metadata
    ani_exclusions = BUILD_FASTANI_INPUTS.out.exclusions
    assembly_stats = CALCULATE_ASSEMBLY_STATS.out.stats
    fastani_matrix = FASTANI.out.matrix
    fastani_raw = FASTANI.out.raw_output
    fastani_log = FASTANI.out.log
    clusters = CLUSTER_ANI.out.clusters
    versions = versions
}
