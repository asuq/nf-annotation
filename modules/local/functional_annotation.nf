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
    if (params.annotation_tools) {
        AnnotationExecution.requireRuntime(workflow.containerEngine, task.container, params.python_container)
    }
    def encoded = groovy.json.JsonOutput.toJson(configuration).bytes.encodeBase64().toString()
    """
    python3 -c 'import base64; from pathlib import Path; Path("configuration.json").write_bytes(base64.b64decode("${encoded}"))'
    python3 "\$(command -v prepare_annotation_tasks.py)" preflight \
        --config configuration.json --output annotation_preflight.json
    """
}


process PLAN_ANNOTATION_GROUPS {
    label 'process_single'
    container params.python_container
    cache 'deep'
    errorStrategy 'finish'
    maxRetries 0
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: true

    input:
    path samples
    path receipt
    path source_results, name: 'source_results'
    path resource_manifests, stageAs: 'resource_manifests/manifest??.json'

    output:
    path 'annotation_groups.json', emit: groups

    script:
    if (params.annotation_tools) {
        AnnotationExecution.requireRuntime(workflow.containerEngine, task.container, params.python_container)
    }
    def sourceArgs = source_results ? "--source '${source_results}'" : ''
    """
    python3 "\$(command -v prepare_annotation_tasks.py)" groups \
        --samples '${samples}' --preflight '${receipt}' ${sourceArgs} \
        --size '${params.eggnog_group_samples}' --output annotation_groups.json
    """
}

process PLAN_ANNOTATION_PART {
    tag "${meta.key}"
    label 'process_single'
    container params.python_container
    cache 'deep'
    errorStrategy 'finish'
    maxRetries 0

    input:
    tuple val(meta), path(bundles, stageAs: 'bundles/bundle??'), path(previous_inputs, stageAs: 'previous/entry??'), val(source_root)
    path receipt
    path resource_manifests, stageAs: 'resource_manifests/manifest??.json'

    output:
    tuple val(meta), path('planned'), emit: tasks
    path 'planned/annotation_plan.json', emit: plan

    script:
    if (params.annotation_tools) {
        AnnotationExecution.requireRuntime(workflow.containerEngine, task.container, params.python_container)
    }
    def bundleList = groovy.json.JsonOutput.toJson((bundles instanceof Collection ? bundles : [bundles]).collect { it.toString() })
    def members = groovy.json.JsonOutput.toJson(meta.members)
    def sourceArgs = source_root ? "--source '${source_root}'" : ''
    """
    cat > members.json <<'ANNOTATION_MEMBERS'
${members}
ANNOTATION_MEMBERS
    cat > bundle_list.json <<'ANNOTATION_BUNDLES'
${bundleList}
ANNOTATION_BUNDLES
    python3 "\$(command -v prepare_annotation_tasks.py)" plan-part \
        --members members.json --bundle-list bundle_list.json --preflight '${receipt}' \
        --scope '${meta.scope}' ${sourceArgs} --output planned
    """
}

process MERGE_ANNOTATION_PLANS {
    label 'process_single'
    container params.python_container
    cache 'deep'
    errorStrategy 'finish'
    maxRetries 0
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', overwrite: true

    input:
    path samples
    path receipt
    path fragments, stageAs: 'fragments/plan??.json'
    path source_results, name: 'source_results'
    path resource_manifests, stageAs: 'resource_manifests/manifest??.json'

    output:
    path 'annotation_plan.json', emit: plan
    path 'annotation_plan.tsv', emit: table

    script:
    if (params.annotation_tools) {
        AnnotationExecution.requireRuntime(workflow.containerEngine, task.container, params.python_container)
    }
    def fragmentList = groovy.json.JsonOutput.toJson((fragments instanceof Collection ? fragments : [fragments]).collect { it.toString() })
    def sourceArgs = source_results ? "--source '${source_results}'" : ''
    """
    cat > fragment_list.json <<'ANNOTATION_FRAGMENTS'
${fragmentList}
ANNOTATION_FRAGMENTS
    python3 "\$(command -v prepare_annotation_tasks.py)" merge-plans \
        --samples '${samples}' --preflight '${receipt}' --fragment-list fragment_list.json \
        ${sourceArgs} --output merged
    cp merged/annotation_plan.json annotation_plan.json
    cp merged/annotation_plan.tsv annotation_plan.tsv
    """
}

process ANNOTATION_SEARCH {
    tag { meta.batch_id && meta.batch_id != 'NA' ? "eggnog-batch:${meta.batch_id.substring(0, 12)}" : "${meta.accession}:${meta.tool}" }
    cpus { meta.cpus as int }
    memory { "${meta.memory_gib} GB" }
    time { params.max_time }
    maxForks params.annotation_max_forks
    container { meta.container }
    publishDir "${params.outdir}/pipeline_info/annotation_storage", mode: 'copy', overwrite: true,
        saveAs: { name -> name == 'annotation_storage.json' ? "${meta.task_directory}.json" : null }
    // An abruptly terminated task may not emit raw files. The complete plan lets
    // aggregation record that missing result as failed and reject the final gate
    // after independent samples finish. Nextflow retains the task error/trace.
    errorStrategy {
        AnnotationExecution.preserveFailure(task, meta, params.outdir, workflow.start)
        'ignore'
    }
    maxRetries 0

    input:
    tuple val(meta), path(task_directory, name: 'task'), path(bundle, name: 'bundle'), path(resource, name: 'resource')

    output:
    tuple val(meta), path('raw'), path('annotation_storage.json'), emit: raw_results

    script:
    AnnotationExecution.requireRuntime(workflow.containerEngine, task.container, meta.container)
    def cleanupPolicy = (params.annotation_cleanup ?: 'success').toString()
    """
    # Nextflow invokes the durable script by absolute path even under scratch.
    durable_work="\$(dirname "\${BASH_SOURCE[0]}")"
    test -f "\${durable_work}/.command.run"
    test '${task.cpus}' -eq '${meta.cpus}'
    test '${task.memory.toBytes()}' -eq '${Math.round((meta.memory_gib as double) * 1024 * 1024 * 1024)}'
    finalize_storage() {
      native_exit=\$?
      trap - EXIT TERM INT
      set +e
      bash "\$(command -v annotation_storage.sh)" finalize \
        --execution-root "\$PWD" \
        --durable-work "\${durable_work}" \
        --task-directory '${meta.task_directory}' \
        --tool '${meta.tool}' \
        --batch-id '${meta.batch_id ?: 'NA'}' \
        --native-exit "\${native_exit}" \
        --policy '${cleanupPolicy}'
      storage_exit=\$?
      if [[ "\${native_exit}" -ne 0 ]]; then
          exit "\${native_exit}"
      fi
      exit "\${storage_exit}"
    }
    trap finalize_storage EXIT
    trap 'exit 143' TERM
    trap 'exit 130' INT
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
    if (params.annotation_tools) {
        AnnotationExecution.requireRuntime(workflow.containerEngine, task.container, params.python_container)
    }
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
    if (params.annotation_tools) {
        AnnotationExecution.requireRuntime(workflow.containerEngine, task.container, params.python_container)
    }
    """
    python3 "\$(command -v prepare_annotation_tasks.py)" reuse --task task --output result
    """
}

process COMPLETE_EGGNOG_BATCH {
    tag "${meta.batch_id.substring(0, 12)}:${meta.action}"
    label 'process_single'
    container params.python_container
    errorStrategy 'finish'
    maxRetries 0
    publishDir params.outdir, mode: 'copy', overwrite: true,
        saveAs: { name -> name.startsWith('published/') ? name.substring(10) : null }

    input:
    tuple val(meta), path(task_directory, name: 'task'), path(batch, name: 'batch'), path(resource, name: 'resource'), path(raw, name: 'raw')

    output:
    path 'published/annotation_batches/*', emit: native_batch
    path 'published/samples/*/annotation/eggnog', emit: members

    script:
    if (params.annotation_tools) {
        AnnotationExecution.requireRuntime(workflow.containerEngine, task.container, params.python_container)
    }
    """
    python3 "\$(command -v prepare_annotation_tasks.py)" complete-batch \
        --task task --raw raw --batch batch --resource resource --output published
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
    path batches, stageAs: 'batches/batch???'

    output:
    path 'report/tables/*', emit: tables
    path 'report/annotation_results.json', emit: manifest
    path 'report/annotation_acceptance.json', emit: acceptance

    script:
    if (params.annotation_tools) {
        AnnotationExecution.requireRuntime(workflow.containerEngine, task.container, params.python_container)
    }
    def bundleList = groovy.json.JsonOutput.toJson((bundles instanceof Collection ? bundles : [bundles]).collect { it.toString() })
    def resultList = groovy.json.JsonOutput.toJson((results instanceof Collection ? results : [results]).collect { it.toString() })
    def batchList = groovy.json.JsonOutput.toJson((batches instanceof Collection ? batches : [batches]).collect { it.toString() })
    """
    cat > bundle_list.json <<'ANNOTATION_BUNDLE_LIST'
${bundleList}
ANNOTATION_BUNDLE_LIST
    cat > result_list.json <<'ANNOTATION_RESULT_LIST'
${resultList}
ANNOTATION_RESULT_LIST
    cat > batch_list.json <<'ANNOTATION_BATCH_LIST'
${batchList}
ANNOTATION_BATCH_LIST
    python3 "\$(command -v aggregate_annotations.py)" --plan '${plan}' \
        --master upstream_master.tsv --sample-status upstream_sample_status.tsv \
        --bundle-list bundle_list.json --result-list result_list.json \
        --batch-list batch_list.json --output report
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
    if (params.annotation_tools) {
        AnnotationExecution.requireRuntime(workflow.containerEngine, task.container, params.python_container)
    }
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
