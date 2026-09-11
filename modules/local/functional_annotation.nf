/* Shared native execution, validation and publication for every core tool. */
process ANNOTATION_PREFLIGHT {
    label 'process_single'
    container params.python_container
    cache false
    errorStrategy 'finish'
    maxRetries 0
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: true

    input:
    val configuration
    path resource_inputs, stageAs: 'resource_inputs/resource??'
    path runtime_inputs, stageAs: 'runtime_inputs/runtime??.sif'

    output:
    path 'annotation_preflight.json', emit: receipt

    script:
    def encoded = groovy.json.JsonOutput.toJson(configuration).bytes.encodeBase64().toString()
    """
    python3 -c 'import base64; from pathlib import Path; Path("configuration.json").write_bytes(base64.b64decode("${encoded}"))'
    python3 "\$(command -v prepare_annotation_tasks.py)" preflight \
        --config configuration.json --output annotation_preflight.json
    """
}

process PLAN_ANNOTATIONS {
    label 'process_single'
    container params.python_container
    // The preflight receipt is recreated after full resource checks. Hash its
    // content and the published inputs so an identical plan keeps stable paths.
    cache 'deep'
    errorStrategy 'finish'
    maxRetries 0
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: true,
        saveAs: { name -> name in ['annotation_plan.tsv', 'annotation_plan.json'] ? name : null }

    input:
    path samples
    path receipt
    path bundles, stageAs: 'bundles/bundle??'
    path source_results, name: 'source_results'

    output:
    path 'planned', emit: tasks
    path 'annotation_plan.json', emit: plan
    path 'annotation_plan.tsv', emit: table

    script:
    def bundleArgs = (bundles instanceof Collection ? bundles : [bundles]).collect { "--bundle '${it}'" }.join(' ')
    def sourceArgs = source_results ? "--source '${source_results}'" : ''
    """
    python3 "\$(command -v prepare_annotation_tasks.py)" plan \
        --samples '${samples}' --preflight '${receipt}' ${bundleArgs} \
        ${sourceArgs} --output planned
    cp planned/annotation_plan.json annotation_plan.json
    cp planned/annotation_plan.tsv annotation_plan.tsv
    """
}

process ANNOTATION_SEARCH {
    tag "${meta.accession}:${meta.tool}"
    cpus { meta.cpus as int }
    memory { "${meta.memory_gib} GB" }
    time { params.max_time }
    maxForks params.annotation_max_forks
    container { meta.container }
    errorStrategy 'finish'
    maxRetries 0

    input:
    tuple val(meta), path(task_directory, name: 'task'), path(bundle, name: 'bundle'), path(resource, name: 'resource')

    output:
    tuple val(meta), path(task_directory), path(bundle), path(resource), path('raw'), emit: raw_results

    script:
    """
    test '${task.cpus}' -eq '${meta.cpus}'
    test '${task.memory.toBytes()}' -eq '${Math.round((meta.memory_gib as double) * 1024 * 1024 * 1024)}'
    bash task/run.sh
    """
}

process NORMALIZE_ANNOTATION {
    tag "${meta.accession}:${meta.tool}"
    label 'process_single'
    container params.python_container
    errorStrategy 'finish'
    maxRetries 0
    publishDir { "${params.outdir}/samples/${meta.accession}/annotation" }, mode: 'copy', overwrite: true,
        saveAs: { name -> name == 'result' ? meta.tool : null }

    input:
    tuple val(meta), path(task_directory, name: 'task'), path(bundle, name: 'bundle'), path(resource, name: 'resource'), path(raw, name: 'raw')

    output:
    tuple val(meta), path('result'), emit: result

    script:
    def resourceArg = resource ? "--resource '${resource}'" : ''
    """
    python3 "\$(command -v prepare_annotation_tasks.py)" normalize \
        --task task --raw raw --bundle bundle ${resourceArg} --output result
    """
}

process REUSE_ANNOTATION {
    tag "${meta.accession}:${meta.tool}"
    label 'process_single'
    container params.python_container
    errorStrategy 'finish'
    maxRetries 0
    publishDir { "${params.outdir}/samples/${meta.accession}/annotation" }, mode: 'copy', overwrite: true,
        saveAs: { name -> name == 'result' ? meta.tool : null }

    input:
    tuple val(meta), path(task_directory, name: 'task')

    output:
    tuple val(meta), path('result'), emit: result

    script:
    """
    python3 "\$(command -v prepare_annotation_tasks.py)" reuse --task task --output result
    """
}

process AGGREGATE_ANNOTATIONS {
    label 'process_single'
    container params.python_container
    cache false
    errorStrategy 'finish'
    maxRetries 0
    publishDir params.outdir, mode: 'copy', overwrite: true,
        saveAs: { name -> name.startsWith('report/') ? name.substring(7) : null }

    input:
    path plan
    path master, name: 'upstream_master.tsv'
    path sample_status, name: 'upstream_sample_status.tsv'
    path bundles, stageAs: 'bundles/bundle??'
    path results, stageAs: 'results/result???'

    output:
    path 'report/tables/*', emit: tables
    path 'report/annotation_results.json', emit: manifest
    path 'report/annotation_acceptance.json', emit: acceptance

    script:
    def bundleArgs = (bundles instanceof Collection ? bundles : [bundles]).collect { "--bundle '${it}'" }.join(' ')
    def resultArgs = (results instanceof Collection ? results : [results]).collect { "--result '${it}'" }.join(' ')
    """
    python3 "\$(command -v aggregate_annotations.py)" --plan '${plan}' \
        --master upstream_master.tsv --sample-status upstream_sample_status.tsv \
        ${bundleArgs} ${resultArgs} --output report
    """
}

process ANNOTATION_ACCEPTANCE {
    label 'process_single'
    container params.python_container
    errorStrategy 'finish'
    maxRetries 0

    input:
    path acceptance
    path reporting_finished

    output:
    path 'annotation_complete.txt', emit: complete

    script:
    """
    python3 - '${acceptance}' <<'PY'
    import json
    import sys
    from pathlib import Path
    record = json.loads(Path(sys.argv[1]).read_text())
    if not record['complete']:
        for failure in record['failures']:
            print(f"{failure['accession']} / {failure['tool']}: {failure['reason']}", file=sys.stderr)
        sys.exit('Requested annotation analyses did not all succeed; published statuses describe the incomplete run')
    Path('annotation_complete.txt').write_text('All requested annotation analyses passed validation.\\n')
    PY
    """
}
