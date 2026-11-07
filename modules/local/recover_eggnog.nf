/* Recover only unfinished annotations, using the unchanged pinned native runtime. */
process RECOVER_EGGNOG_BATCH {
    tag "${meta.batch_id.substring(0, 12)}"
    cpus { meta.cpus as int }
    memory { "${meta.memory_gib} GB" }
    time { params.max_time }
    maxForks params.annotation_max_forks
    container { meta.container }
    scratch true
    stageOutMode 'copy'
    publishDir "${params.outdir}/pipeline_info/annotation_storage", mode: 'copy', overwrite: true,
        saveAs: { name -> name == 'annotation_storage.json' ? "${meta.task_directory}.json" : null }
    errorStrategy {
        AnnotationExecution.preserveFailure(task, meta, params.outdir, workflow.start)
        'ignore'
    }
    maxRetries 0

    input:
    tuple val(meta), path(task_directory, name: 'task'), path(bundle, name: 'bundle'), path(resource, name: 'resource'), path(previous, name: 'previous'), path(members, stageAs: 'members/member????')

    output:
    tuple val(meta), path('raw'), path('annotation_storage.json'), emit: raw_results

    script:
    AnnotationExecution.requireRuntime(workflow.containerEngine, task.container, meta.container)
    def cleanupPolicy = (params.annotation_cleanup ?: 'success').toString()
    """
    durable_work="\$(dirname "\${BASH_SOURCE[0]}")"
    test -f "\${durable_work}/.command.run"
    test '${task.cpus}' -eq '${meta.cpus}'
    test '${task.memory.toBytes()}' -eq '${Math.round((meta.memory_gib as double) * 1024 * 1024 * 1024)}'
    finalize_storage() {
      native_exit=\$?
      trap - EXIT TERM INT
      set +e
      bash "\$(command -v annotation_storage.sh)" finalize \
        --execution-root "\$PWD" --durable-work "\${durable_work}" \
        --task-directory '${meta.task_directory}' --tool '${meta.tool}' \
        --batch-id '${meta.batch_id}' --native-exit "\${native_exit}" \
        --policy '${cleanupPolicy}'
      storage_exit=\$?
      if [[ "\${native_exit}" -ne 0 ]]; then exit "\${native_exit}"; fi
      exit "\${storage_exit}"
    }
    trap finalize_storage EXIT
    trap 'exit 143' TERM
    trap 'exit 130' INT
    # Nextflow scratch must reside on a worker-local filesystem for this remedy.
    case "\$(stat -f -c %T .)" in
      ext2/ext3|ext4|xfs|btrfs|zfs|tmpfs|overlayfs) ;;
      *) echo 'eggNOG recovery requires worker-local scratch' >&2; exit 1 ;;
    esac
    python3 "\$(command -v recover_eggnog_batch.py)" \
      --task task --batch bundle --resource resource --previous previous --members members/member*
    """
}
