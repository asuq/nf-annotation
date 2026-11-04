/* Directory cache keys alone do not detect edits to nested publication files. */
process FINGERPRINT_REUSE_SOURCES {
    label 'process_cohort'
    stageInMode 'symlink'
    container params.python_container
    cache false
    publishDir "${params.outdir}/pipeline_info/source_fingerprinting", mode: 'copy', overwrite: true, pattern: 'fingerprint_*.json'
    errorStrategy 'finish'
    maxRetries 0

    input:
    path sources, stageAs: 'sources/source????'
    path dependencies, stageAs: 'dependencies/root????'
    path previous, name: 'previous_validation.json'

    output:
    path 'source_fingerprints.json', emit: fingerprints
    path 'validated_sources.json', emit: receipts
    path 'fingerprint_metrics.json', emit: metrics
    path 'fingerprint_progress.json', emit: progress

    script:
    def previousArg = previous ? '--previous previous_validation.json' : ''
    def workers = Math.min(Math.min(task.cpus as int, params.source_fingerprint_workers as int), 8)
    def roots = groovy.json.JsonOutput.toJson((sources instanceof Collection ? sources : [sources]).collect { it.toString() })
    """
    cat > sources.json <<'SOURCE_ROOTS'
${roots}
SOURCE_ROOTS
    python3 "\$(command -v fingerprint_sources)" \
        --sources sources.json --output source_fingerprints.json --metrics fingerprint_metrics.json \
        --workers ${workers} --progress fingerprint_progress.json --receipts validated_sources.json ${previousArg}
    """
}
