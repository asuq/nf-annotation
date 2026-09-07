include { PROKKA } from '../../modules/local/prokka'
include { CCFINDER } from '../../modules/local/ccfinder'
include { CODETTA } from '../../modules/local/codetta'
include { SUMMARISE_CODETTA } from '../../modules/local/summarise_codetta'
include { SUMMARISE_CCFINDER } from '../../modules/local/summarise_ccfinder'
include { PADLOC } from '../../modules/local/padloc'
include { EGGNOG } from '../../modules/local/eggnog'

/*
 * Run the gcode-gated annotation branch by joining staged genomes to the
 * per-sample CheckM2 summary that assigned the translation table.
 */
workflow PER_SAMPLE_ANNOTATION {
    take:
    sample_genomes
    gcode_summaries
    codetta_db
    eggnog_db

    main:
    def eggnogOnlyAccessions = null
    if (params.eggnog_only_accessions != null) {
        def rawAccessions = params.eggnog_only_accessions instanceof Collection
            ? (params.eggnog_only_accessions as Collection)
            : params.eggnog_only_accessions.toString().split(',')
        eggnogOnlyAccessions = rawAccessions
            .collect { it.toString().trim() }
            .findAll { !it.isEmpty() } as Set
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
    CCFINDER(annotation_candidates)
    SUMMARISE_CCFINDER(CCFINDER.out.result_json)
    PADLOC(PROKKA.out.padloc_inputs)

    eggnog_inputs = PROKKA.out.eggnog_inputs
    if (eggnogOnlyAccessions != null && !eggnogOnlyAccessions.isEmpty()) {
        eggnog_inputs = eggnog_inputs.filter { item ->
            def values = unpackTuple.call(item, 'eggnog_inputs', 2)
            def meta = values[0]
            eggnogOnlyAccessions.contains(meta.accession.toString())
        }
    }

    EGGNOG(eggnog_inputs.combine(eggnog_db))

    eggnogSkippedRows = gcode_by_accession
        .filter { item ->
            item[1] in ['4', '11'] && eggnogOnlyAccessions != null &&
                !eggnogOnlyAccessions.isEmpty() && !eggnogOnlyAccessions.contains(item[0])
        }
        .map { item -> "${item[0]}\tskipped\teggnog_short_circuit\t0\t0\t0" }

    versions = PROKKA.out.versions
        .mix(CODETTA.out.versions)
        .mix(SUMMARISE_CODETTA.out.versions)
        .mix(CCFINDER.out.versions)
        .mix(SUMMARISE_CCFINDER.out.versions)
        .mix(PADLOC.out.versions)
        .mix(EGGNOG.out.versions)

    emit:
    codetta = CODETTA.out.results
    codetta_summary = SUMMARISE_CODETTA.out.summary
    prokka = PROKKA.out.results
    ccfinder = CCFINDER.out.results
    ccfinder_summary = SUMMARISE_CCFINDER.out.summaries
    padloc = PADLOC.out.results
    eggnog = EGGNOG.out.results
    eggnog_skips = eggnogSkippedRows
    versions = versions
}
