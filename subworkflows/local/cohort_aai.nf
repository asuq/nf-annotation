include { PREPARE_AAI_BATCH; FASTAAI_COMPARE; CLUSTER_AAI } from '../../modules/local/fastaai'

workflow COHORT_AAI {
    take:
    samples
    metadata
    checkm2
    busco
    stats
    bundles
    primary_busco_column
    main:
    if (SamplePolicy.stage(params) == 'qc') {
        summaries = Channel.value(file("${projectDir}/assets/tables/headers/aai_not_requested.tsv"))
        representatives = Channel.empty()
        versions = Channel.empty()
        evidence = Channel.empty()
    } else {
        requested = samples.splitCsv(header: true, sep: '\t').toList()
        batches = requested.map { [it] }.combine(bundles.map { [it] }).flatMap { item ->
            def rows = item[0].sort(false) { a, b -> a.accession <=> b.accession }
            def byAccession = item[1].collectEntries { path ->
                def manifest = new groovy.json.JsonSlurper().parse(path.resolve('bundle.json').toFile())
                [(manifest.accession): path]
            }
            def size = Math.max(1, Math.ceil(rows.size() / 8.0) as int)
            rows.collate(size).withIndex().collect { group, index ->
                tuple(index, group, group.collect { byAccession[it.accession] }.findAll { it != null })
            }
        }
        destination = file(params.outdir)
        destination.toFile().mkdirs()
        sources = params.update_from ? PublishedSources.paths(params.update_from) : []
        PREPARE_AAI_BATCH(batches, metadata, checkm2.first(), busco, stats.first(),
            Channel.value(sources), Channel.value(ReferenceInputs.paths(sources)), Channel.value(destination), primary_busco_column)
        domainRecords = PREPARE_AAI_BATCH.out.batch.toList().flatMap { paths ->
            def rows = paths.collectMany { new groovy.json.JsonSlurper().parse(it.toFile()).samples }
            rows.groupBy { it.domain }.collect { domain, members ->
                def ordered = members.sort(false) { a, b -> a.accession <=> b.accession }
                tuple(domain, ordered, ordered.findAll { it.status == 'ready' }.collect { file(it.profile) })
            }
        }
        FASTAAI_COMPARE(domainRecords)
        CLUSTER_AAI(FASTAAI_COMPARE.out.comparison.toSortedList { a, b -> a[0] <=> b[0] }.map { rows -> rows.collect { it[1] } }, samples)
        evidence = FASTAAI_COMPARE.out.comparison.map { domain, directory -> tuple("cohort/fastaai/${domain}", directory) }
            .mix(CLUSTER_AAI.out.evidence.map { tuple('cohort/aai_clusters', it) })
            .mix(PREPARE_AAI_BATCH.out.batch.toList().map { tuple('cohort/fastaai/profiles', file("${params.outdir}/cohort/fastaai/profiles")) })
        summaries = CLUSTER_AAI.out.summary
        representatives = CLUSTER_AAI.out.representatives
        versions = PREPARE_AAI_BATCH.out.versions.mix(FASTAAI_COMPARE.out.versions).mix(CLUSTER_AAI.out.versions)
    }
    emit:
    evidence = evidence
    summary = summaries
    representatives = representatives
    versions = versions
}
