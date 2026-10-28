process REFERENCE_ACCEPTANCE {
    label 'process_cohort'
    container params.python_container
    stageInMode 'symlink'
    cache false
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

    script:
    def snapshotArg = snapshot ? '--snapshot source_snapshot.json' : ''
    """
    python3 "\$(command -v verify_reference_acceptance.py)" --acceptance '${acceptance}' \
        --publication publication ${snapshotArg} --output annotation_complete.txt
    """
}
