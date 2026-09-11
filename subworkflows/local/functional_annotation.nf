include { ANNOTATION_PREFLIGHT } from '../../modules/local/functional_annotation'
include { PLAN_ANNOTATIONS } from '../../modules/local/functional_annotation'
include { ANNOTATION_SEARCH } from '../../modules/local/functional_annotation'
include { NORMALIZE_ANNOTATION } from '../../modules/local/functional_annotation'
include { REUSE_ANNOTATION } from '../../modules/local/functional_annotation'

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

/* Run, reuse or renormalize each task from an explicit full-cohort plan. */
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
    PLAN_ANNOTATIONS(samples, preflight, collectedBundles, source_results)
    plannedTasks = PLAN_ANNOTATIONS.out.table
        .splitCsv(header: true, sep: '\t')
        .filter { row -> row.action != 'skip' }
        .combine(PLAN_ANNOTATIONS.out.tasks)
        .map { item ->
            def row = item[0]
            def taskDir = item[1].resolve(row.task_directory)
            def resource = row.resource == 'NA' ? [] : file(row.resource, checkIfExists: true)
            tuple(row, taskDir, file(row.bundle, checkIfExists: true), resource)
        }
    ANNOTATION_SEARCH(plannedTasks.filter { item -> item[0].action == 'run' })
    normalizationInputs = ANNOTATION_SEARCH.out.raw_results.mix(
        plannedTasks.filter { item -> item[0].action == 'renormalize' }
            .map { item -> tuple(item[0], item[1], item[2], item[3], item[1].resolve('previous/raw')) }
    )
    NORMALIZE_ANNOTATION(normalizationInputs)
    REUSE_ANNOTATION(plannedTasks.filter { item -> item[0].action == 'reuse' }.map { item -> tuple(item[0], item[1]) })

    emit:
    results = NORMALIZE_ANNOTATION.out.result.mix(REUSE_ANNOTATION.out.result)
    plan = PLAN_ANNOTATIONS.out.plan
    bundle_files = collectedBundles
}
