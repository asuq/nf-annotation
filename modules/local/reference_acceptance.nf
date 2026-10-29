process REFERENCE_ACCEPTANCE {
    label 'process_cohort'
    container params.python_container
    stageInMode 'symlink'
    cache false
    publishDir "${params.outdir}/pipeline_info/reference_acceptance", mode: 'copy', overwrite: true, pattern: 'reference_progress.json'
    errorStrategy 'finish'
    maxRetries 0

    input:
    path acceptance
    path reporting_finished
    path publication, name: 'publication'
    path dependencies, stageAs: 'dependencies/root????'
    path snapshot, name: 'source_snapshot.json'

    output:
    path 'annotation_complete.txt', emit: complete
    path 'reference_progress.json', emit: progress

    script:
    def workers = Math.min(Math.min(task.cpus as int, params.source_fingerprint_workers as int), 8)
    def snapshotArg = snapshot ? '--snapshot source_snapshot.json' : '--source-validation publication/tables/source_validation.json'
    """
    python3 "\$(command -v verify_reference_acceptance.py)" --acceptance '${acceptance}' \
        --publication publication ${snapshotArg} --output annotation_complete.txt --workers ${workers} --progress reference_progress.json
    """
}
