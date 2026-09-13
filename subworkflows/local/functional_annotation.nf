include { ANNOTATION_PREFLIGHT } from '../../modules/local/functional_annotation'
include { PLAN_ANNOTATIONS } from '../../modules/local/functional_annotation'
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
    collectedBundles = bundles.toList().map { items ->
        items.sort(false) { left, right -> left[0].accession <=> right[0].accession }.collect { it[1] }
    }
    // The planner rechecks these small files at their recorded absolute paths.
    // Staging them makes their parent directories visible inside its container.
    resourceManifests = preflight.map { receipt ->
        def record = new groovy.json.JsonSlurper().parse(receipt.toFile())
        record.tools.values().findAll { it.resource_path }.collect {
            file("${it.resource_path}/annotation_resource.json", checkIfExists: true)
        }
    }
    PLAN_ANNOTATIONS(samples, preflight, collectedBundles, source_results, resourceManifests)
    checkedTaskDirectories = PLAN_ANNOTATIONS.out.tasks.combine(PLAN_ANNOTATIONS.out.plan)
        .map { item ->
            def directory = item[0]
            def plan = item[1]
            def record = new groovy.json.JsonSlurper().parse(plan.toFile())
            def expected = record.tasks.findAll { it.batch_id }.collect { it.batch_id } as Set
            def published = file("${params.outdir}/annotation_batches")
            if (published.exists() || java.nio.file.Files.isSymbolicLink(published)) {
                if (java.nio.file.Files.isSymbolicLink(published) || !published.toFile().isDirectory()) {
                    error 'Published annotation_batches must be an ordinary directory; use a fresh --outdir'
                }
                def incompatible = published.toFile().listFiles().findAll {
                    !(it.name in expected) || !it.isDirectory() || java.nio.file.Files.isSymbolicLink(it.toPath())
                }
                if (incompatible) {
                    error 'Existing annotation_batches are incompatible with the current plan; use a fresh --outdir'
                }
            }
            directory
        }
    plannedTasks = PLAN_ANNOTATIONS.out.table
        .splitCsv(header: true, sep: '\t')
        .filter { row -> row.action != 'skip' }
        .unique { row -> row.task_directory }
        .combine(checkedTaskDirectories)
        .map { item ->
            def row = item[0]
            def taskDir = item[1].resolve(row.task_directory)
            def resource = row.resource == 'NA' ? [] : file(row.resource, checkIfExists: true)
            tuple(row, taskDir, file(row.bundle, checkIfExists: true), resource)
        }
    ANNOTATION_SEARCH(plannedTasks.filter { item -> item[0].action == 'run' })
    // Failed native tasks retain a real nonzero Nextflow exit and are therefore
    // rerun on resume. Collect their saved diagnostics after all searches finish.
    failedRawResults = plannedTasks.filter { item -> item[0].action == 'run' }
        .combine(ANNOTATION_SEARCH.out.raw_results.count())
        .map { item ->
            def raw = AnnotationExecution.failedNativeRaw(params.outdir, workflow.start, item[0])
            raw ? tuple(item[0], item[1], item[2], item[3], raw) : null
        }
    nativeResults = ANNOTATION_SEARCH.out.raw_results.mix(failedRawResults)
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
    plan = PLAN_ANNOTATIONS.out.plan
    bundle_files = collectedBundles
}
