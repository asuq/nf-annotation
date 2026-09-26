include { ANNOTATION_PREFLIGHT } from '../../modules/local/functional_annotation'
include { PLAN_ANNOTATION_GROUPS; PLAN_ANNOTATION_PART; MERGE_ANNOTATION_PLANS } from '../../modules/local/functional_annotation'
include { ANNOTATION_SEARCH } from '../../modules/local/functional_annotation'
include { NORMALIZE_ANNOTATION } from '../../modules/local/functional_annotation'
include { REUSE_ANNOTATION } from '../../modules/local/functional_annotation'
include { COMPLETE_EGGNOG_BATCH } from '../../modules/local/functional_annotation'

/* Resolve all requested methods before starting expensive per-sample searches. */
workflow ANNOTATION_RESOURCES {
    main:
    ['eggnog_only_accessions', 'eggnog_extra_args', 'padloc_extra_args'].each { removed ->
        if (params.containsKey(removed)) {
            error "--${removed} has been removed; use the shared --annotation_tools workflow and its fixed scientific methods"
        }
    }
    def supported = ['eggnog', 'cogclassifier', 'pfam', 'kofam', 'padloc']
    if (!params.eggnog_group_samples.toString().isInteger() || (params.eggnog_group_samples as int) < 1) {
        error '--eggnog_group_samples must be a positive integer'
    }
    def selected = params.annotation_tools ? params.annotation_tools.toString().split(',') as List : []
    if (selected.unique(false).size() != selected.size() || selected.any { !(it in supported) }) {
        error "--annotation_tools must be a unique comma-separated subset of ${supported.join(',')}"
    }
    if (selected) {
        AnnotationExecution.requireEngine(workflow.containerEngine)
    }
    def cpus = Math.min(params.annotation_cpus as int, params.max_cpus as int)
    def memory = Math.min((params.annotation_memory as nextflow.util.MemoryUnit).toBytes(), (params.max_memory as nextflow.util.MemoryUnit).toBytes()) / (1024.0 * 1024 * 1024)
    def config = [annotation_tools: selected.join(','), helper_container: params.python_container, tools: [:]]
    def resourcePaths = []
    def runtimePaths = []
    selected.each { tool ->
        def container = params["${tool}_container"]
        if (!container) { error "--${tool}_container requires an immutable container reference or SIF path" }
        def resource = tool == 'padloc' ? null : params["${tool}_db"]
        if (tool != 'padloc' && !resource) { error "--${tool}_db is required for an enabled annotation tool" }
        if (resource) { resourcePaths.add(file(resource, checkIfExists: true)) }
        if (container.toString().endsWith('.sif')) { runtimePaths.add(file(container, checkIfExists: true)) }
        config.tools[tool] = [container: container.toString(), resource: resource ? file(resource).toAbsolutePath().toString() : null, cpus: cpus, memory_gib: memory]
    }
    if (params.python_container.toString().endsWith('.sif')) { runtimePaths.add(file(params.python_container, checkIfExists: true)) }
    ANNOTATION_PREFLIGHT(Channel.value(config), Channel.value(resourcePaths), Channel.value(runtimePaths.unique()))

    emit:
    receipt = ANNOTATION_PREFLIGHT.out.receipt
}

/* Keep the full sample/tool plan while executing each shared batch once. */
workflow FUNCTIONAL_ANNOTATION {
    take:
    samples
    bundles
    preflight
    source_results

    main:
    // This list is available from the manifest before sample jobs run.
    declaredAccessions = samples.splitCsv(header: true, sep: '\t').toList().map { rows ->
        rows.collect { it.accession.toString() } as Set
    }
    seenBundles = [] as Set
    checkedBundles = bundles.combine(declaredAccessions).map { item ->
        def accession = item[0].accession.toString()
        if (!(accession in item[2]) || !seenBundles.add(accession)) {
            error "Duplicate or foreign protein bundle: ${accession}"
        }
        tuple(item[0], item[1])
    }
    // Collection is now only for final reporting, never for launching searches.
    collectedBundles = checkedBundles.toList().map { items ->
        items.sort(false) { left, right -> left[0].accession <=> right[0].accession }.collect { it[1] }
    }
    sharedPreflight = preflight.first()
    resourceManifests = sharedPreflight.map { receipt ->
        def record = new groovy.json.JsonSlurper().parse(receipt.toFile())
        record.tools.values().findAll { it.resource_path }.collect {
            file("${it.resource_path}/annotation_resource.json", checkIfExists: true)
        }
    }
    individualInputs = checkedBundles.map { item ->
        def accession = item[0].accession.toString()
        tuple([key: accession, scope: 'individual', members: [accessions: [accession]]], [item[1]])
    }
    PLAN_ANNOTATION_GROUPS(samples, sharedPreflight, source_results, resourceManifests)
    groupIndex = PLAN_ANNOTATION_GROUPS.out.groups.first().combine(sharedPreflight).map { item ->
        def record = new groovy.json.JsonSlurper().parse(item[0].toFile())
        def checked = new groovy.json.JsonSlurper().parse(item[1].toFile())
        def published = file("${params.outdir}/annotation_batches")
        if (published.exists() || java.nio.file.Files.isSymbolicLink(published)) {
            if (java.nio.file.Files.isSymbolicLink(published) || !published.toFile().isDirectory()) {
                error 'Published annotation_batches must be an ordinary directory; use a fresh --outdir'
            }
            published.toFile().listFiles().each { directory ->
                if (!directory.isDirectory() || java.nio.file.Files.isSymbolicLink(directory.toPath())) {
                    error 'Existing annotation_batches are incompatible with the current plan; use a fresh --outdir'
                }
                def packet = new groovy.json.JsonSlurper().parse(new File(directory, 'batch_result.json'))
                if (packet.batch.search_method != checked.tools.eggnog?.search_method ||
                    packet.batch.members.any { !(it.accession in record.accessions) }) {
                    error 'Existing annotation_batches are incompatible with the current plan; use a fresh --outdir'
                }
            }
        }
        record.groups.collectMany { group ->
            group.accessions.collect { accession -> [accession, group] }
        }.collectEntries()
    }
    groupedInputs = checkedBundles.combine(groupIndex)
        .filter { item -> item[2].containsKey(item[0].accession.toString()) }
        .map { item ->
            def group = item[2][item[0].accession.toString()]
            tuple(groupKey(group, group.accessions.size()), item[1])
        }
        .groupTuple(remainder: true)
        .map { item ->
            def group = item[0].getGroupTarget()
            tuple([key: group.key, scope: 'eggnog', members: group], item[1])
        }
    partInputs = individualInputs.mix(groupedInputs)
        .combine(source_results.map { roots -> [PublishedSources.items(roots)] }).map { item ->
        def meta = item[0]
        def sources = PublishedSources.items(item[2])
        def previous = []
        sources.each { source ->
            meta.members.accessions.each { accession ->
                def tools = meta.scope == 'eggnog' ? ['eggnog'] : ['cogclassifier', 'pfam', 'kofam', 'padloc']
                tools.each { tool ->
                    def path = source.resolve("samples/${accession}/annotation/${tool}")
                    if (path.exists()) { previous.add(path) }
                }
            }
            (meta.members.previous_batch_ids ?: []).each { key ->
                def path = source.resolve("annotation_batches/${key}")
                if (path.exists()) { previous.add(path) }
            }
        }
        tuple(meta, item[1], previous, sources.collect { it.toString() })
    }
    PLAN_ANNOTATION_PART(partInputs, sharedPreflight, resourceManifests)
    fragments = PLAN_ANNOTATION_PART.out.plan.toList().map { paths -> paths.sort(false) { it.toString() } }
    MERGE_ANNOTATION_PLANS(samples, sharedPreflight, fragments, source_results, resourceManifests)
    // Batch members share one task directory. Collapse only inside a ready
    // fragment, so an unrelated unfinished fragment cannot hold up execution.
    seenExecutionKeys = [] as Set
    plannedTasks = PLAN_ANNOTATION_PART.out.tasks.flatMap { item ->
        def directory = item[1]
        def record = new groovy.json.JsonSlurper().parse(directory.resolve('annotation_plan.json').toFile())
        if (item[0].scope == 'eggnog') {
            def expected = record.tasks.findAll { it.batch_id }.collect { it.batch_id } as Set
            def published = file("${params.outdir}/annotation_batches")
            if (published.exists()) {
                published.toFile().listFiles().each { previous ->
                    def packet = new groovy.json.JsonSlurper().parse(new File(previous, 'batch_result.json'))
                    if (packet.batch.members.any { it.accession in record.accessions } && !(previous.name in expected)) {
                        error 'Existing annotation_batches are incompatible with the current plan; use a fresh --outdir'
                    }
                }
            }
        }
        def contexts = record.tasks.findAll { it.action != 'skip' }.collect { original ->
            def row = original.collectEntries { key, value -> [(key): value == null ? 'NA' : value.toString()] }
            def resource = row.resource == 'NA' ? [] : file(row.resource, checkIfExists: true)
            def bundlePath = java.nio.file.Paths.get(row.bundle)
            def stagedBundle = bundlePath.isAbsolute() ? bundlePath : directory.resolve(bundlePath).normalize()
            if (!bundlePath.isAbsolute() && !stagedBundle.startsWith(directory.toAbsolutePath().normalize())) {
                error 'Planned bundle path escapes its staged fragment'
            }
            tuple(row, directory.resolve(row.task_directory), file(stagedBundle, checkIfExists: true), resource)
        }
        contexts.groupBy { it[0].task_directory }.collect { key, items ->
            AnnotationExecution.validateExecutionGroup(key, items)
            items[0]
        }
    }.map { item ->
        def key = item[0].task_directory
        if (!seenExecutionKeys.add(key)) { error "Duplicate functional execution key: ${key}" }
        item
    }
    runTasks = plannedTasks.filter { item -> item[0].action == 'run' }
    ANNOTATION_SEARCH(runTasks)
    plannedExecutionContexts = runTasks.map { item -> tuple(item[0].task_directory, item) }
    nativeExecutionResults = ANNOTATION_SEARCH.out.raw_results
        .map { item -> tuple(item[0].task_directory, tuple(item[0], item[1], item[2])) }
    successfulNativeResults = plannedExecutionContexts.join(nativeExecutionResults, failOnDuplicate: true, remainder: true)
        .filter { joined ->
            if (!joined[1]) { error "Unknown native functional execution key: ${joined[0]}" }
            joined[2] != null
        }
        .map { joined ->
            def key = joined[0]
            def planned = joined[1]
            def nativeResult = joined[2]
            def plannedMeta = planned[0]
            def nativeMeta = nativeResult[0]
            def mismatchedField = ['task_directory', 'tool', 'accession', 'batch_id'].find { field ->
                (plannedMeta[field] ?: 'NA').toString() != (nativeMeta[field] ?: 'NA').toString()
            }
            if (mismatchedField) {
                error "Native functional result identity mismatch for ${key}: ${mismatchedField}"
            }
            tuple(plannedMeta, planned[1], planned[2], planned[3], nativeResult[1])
        }
    // Failed native tasks retain a real nonzero Nextflow exit and are therefore
    // rerun on resume. Collect their saved diagnostics after all searches finish.
    failedRawResults = runTasks
        .combine(ANNOTATION_SEARCH.out.raw_results.count())
        .map { item ->
            def raw = AnnotationExecution.failedNativeRaw(params.outdir, workflow.start, item[0])
            raw ? tuple(item[0], item[1], item[2], item[3], raw) : null
        }
    nativeResults = successfulNativeResults.mix(failedRawResults)
    individualTasks = plannedTasks.filter { item -> item[0].batch_id == 'NA' }
    normalizationInputs = nativeResults
        .filter { item -> item[0].batch_id == 'NA' }.mix(
        individualTasks.filter { item -> item[0].action == 'renormalize' }
            .map { item -> tuple(item[0], item[1], item[2], item[3], item[1].resolve('previous/raw')) }
    )
    NORMALIZE_ANNOTATION(normalizationInputs)
    REUSE_ANNOTATION(individualTasks.filter { item -> item[0].action == 'reuse' }.map { item -> tuple(item[0], item[1]) })
    batchInputs = nativeResults
        .filter { item -> item[0].batch_id != 'NA' }.mix(
        plannedTasks.filter { item -> item[0].batch_id != 'NA' && item[0].action in ['reuse', 'renormalize'] }
            .map { item -> tuple(item[0], item[1], item[2], item[3], item[1].resolve('previous_batch/raw')) }
    )
    COMPLETE_EGGNOG_BATCH(batchInputs)
    batchResults = COMPLETE_EGGNOG_BATCH.out.members.flatMap { members ->
        (members instanceof Collection ? members : [members]).collect { member ->
            tuple([accession: member.parent.parent.name, tool: 'eggnog'], member)
        }
    }

    emit:
    results = NORMALIZE_ANNOTATION.out.result.mix(REUSE_ANNOTATION.out.result, batchResults)
    native_batches = COMPLETE_EGGNOG_BATCH.out.native_batch
    plan = MERGE_ANNOTATION_PLANS.out.plan
    bundle_files = collectedBundles
}
