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
        groupIndex = samples.splitCsv(header: true, sep: '\t').toList().map { rows ->
            def ordered = rows.sort(false) { a, b -> a.accession <=> b.accession }
            def size = Math.max(1, Math.ceil(ordered.size() / 8.0) as int)
            ordered.collate(size).withIndex().collectMany { group, index ->
                group.collect { row -> [row.accession, [sample: row, group: [id: index, size: group.size()]]] }
            }.collectEntries()
        }.first()
        qcRows = checkm2.map { it[1] }.splitCsv(header: true, sep: '\t')
        paired = qcRows.filter { it.Gcode in ['4', '11'] }.map { tuple(it.accession, it) }
            .join(bundles.map { tuple(it[0].accession, it[1]) }, remainder: true, failOnDuplicate: true)
            .map { accession, qc, bundle ->
                if (qc == null) { error "Protein bundle lacks resolved QC: ${accession}" }
                tuple(accession, qc, bundle ?: [])
            }
        ready = paired.mix(qcRows.filter { !(it.Gcode in ['4', '11']) }.map { tuple(it.accession, it, []) })
        batches = ready.combine(groupIndex).map { accession, qc, bundle, index ->
            def item = index[accession]
            if (item == null) { error "Foreign AAI QC sample: ${accession}" }
            tuple(groupKey(item.group, item.group.size as int), item.sample + [qc: qc], bundle)
        }.groupTuple().map { key, records, bundlePaths ->
            tuple(key.getGroupTarget().id, records.sort(false) { a, b -> a.accession <=> b.accession },
                bundlePaths.findAll { it }.sort(false) { it.toString() })
        }
        destination = file(params.outdir)
        destination.toFile().mkdirs()
        sources = params.update_from ? PublishedSources.paths(params.update_from) : []
        PREPARE_AAI_BATCH(batches, metadata,
            Channel.value(sources), Channel.value(ReferenceInputs.paths(sources)), Channel.value(destination), primary_busco_column)
        domainRecords = PREPARE_AAI_BATCH.out.batch.toList().flatMap { paths ->
            def rows = paths.collectMany { new groovy.json.JsonSlurper().parse(it.toFile()).samples }
            rows.groupBy { it.domain }.collect { domain, members ->
                def ordered = members.sort(false) { a, b -> a.accession <=> b.accession }
                tuple(domain, ordered, ordered.findAll { it.status == 'ready' }.collect { file(it.profile) })
            }
        }
        FASTAAI_COMPARE(domainRecords)
        CLUSTER_AAI(FASTAAI_COMPARE.out.comparison.toSortedList { a, b -> a[0] <=> b[0] }.map { rows -> rows.collect { it[1] } }, samples, metadata, stats.first(), busco)
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
