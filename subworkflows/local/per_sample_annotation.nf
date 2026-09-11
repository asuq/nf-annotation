include { PROKKA } from '../../modules/local/prokka'
include { PREPARE_ANNOTATION_BUNDLE } from '../../modules/local/prepare_annotation_bundle'
include { CCFINDER } from '../../modules/local/ccfinder'
include { CODETTA } from '../../modules/local/codetta'
include { SUMMARISE_CODETTA } from '../../modules/local/summarise_codetta'
include { SUMMARISE_CCFINDER } from '../../modules/local/summarise_ccfinder'

/*
 * Run the gcode-gated annotation branch by joining staged genomes to the
 * per-sample CheckM2 summary that assigned the translation table.
 */
workflow PER_SAMPLE_ANNOTATION {
    take:
    sample_genomes
    gcode_summaries
    codetta_db

    main:
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

    sample_by_accession = sample_genomes.map { item ->
        def values = unpackTuple.call(item, 'sample_genomes', 2)
        def meta = values[0]
        def genome = values[1]
        tuple(meta.accession, meta, genome)
    }

    gcode_by_accession = gcode_summaries
        .map { item ->
            def values = unpackTuple.call(item, 'gcode_summaries', 2)
            values[1]
        }
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(row.accession.toString(), row.Gcode.toString())
        }

    annotation_candidates = sample_by_accession
        .join(gcode_by_accession)
        .filter { item ->
            def values = unpackTuple.call(item, 'annotation_candidates_join', 4)
            def gcode = values[3]
            ['4', '11'].contains(gcode.toString())
        }
        .map { item ->
            def values = unpackTuple.call(item, 'annotation_candidates', 4)
            def meta = values[1]
            def genome = values[2]
            def gcode = values[3]
            tuple(meta, genome, gcode.toString())
        }

    CODETTA(sample_genomes.combine(codetta_db))
    SUMMARISE_CODETTA(CODETTA.out.summary_input)

    PROKKA(annotation_candidates)
    PREPARE_ANNOTATION_BUNDLE(PROKKA.out.bundle_inputs)
    CCFINDER(annotation_candidates)
    SUMMARISE_CCFINDER(CCFINDER.out.result_json)
    versions = PROKKA.out.versions
        .mix(PREPARE_ANNOTATION_BUNDLE.out.versions)
        .mix(CODETTA.out.versions)
        .mix(SUMMARISE_CODETTA.out.versions)
        .mix(CCFINDER.out.versions)
        .mix(SUMMARISE_CCFINDER.out.versions)

    emit:
    codetta = CODETTA.out.results
    codetta_summary = SUMMARISE_CODETTA.out.summary
    prokka = PROKKA.out.results
    bundles = PREPARE_ANNOTATION_BUNDLE.out.bundle
    ccfinder = CCFINDER.out.results
    ccfinder_summary = SUMMARISE_CCFINDER.out.summaries
    versions = versions
}
