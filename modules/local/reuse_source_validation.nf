/* Directory cache keys alone do not detect edits to nested publication files. */
process FINGERPRINT_REUSE_SOURCES {
    label 'process_cohort'
    stageInMode 'symlink'
    container params.python_container
    cache false
    errorStrategy 'finish'
    maxRetries 0

    input:
    path sources, stageAs: 'sources/source????'

    output:
    path 'source_fingerprints.json', emit: fingerprints
    path 'fingerprint_metrics.json', emit: metrics

    script:
    def roots = groovy.json.JsonOutput.toJson((sources instanceof Collection ? sources : [sources]).collect { it.toString() })
    """
    cat > sources.json <<'SOURCE_ROOTS'
${roots}
SOURCE_ROOTS
    python3 "\$(command -v reuse_source_validation.py)" fingerprint \
        --sources sources.json --output source_fingerprints.json --metrics fingerprint_metrics.json
    """
}

process VALIDATE_REUSE_SOURCES {
    tag "${records.size()} source(s)"
    label 'process_cohort'
    stageInMode 'symlink'
    container params.python_container
    cache 'lenient'
    maxForks params.ani_import_max_forks
    errorStrategy 'finish'
    maxRetries 0

    input:
    tuple val(records), path(sources, stageAs: 'sources/source????')
    val validator_code

    output:
    path 'validated_sources.json', emit: receipts

    script:
    def roots = groovy.json.JsonOutput.toJson((sources instanceof Collection ? sources : [sources]).collect { it.toString() })
    def fingerprints = groovy.json.JsonOutput.toJson(records)
    """
    # validator_code=${validator_code}
    cat > sources.json <<'SOURCE_ROOTS'
${roots}
SOURCE_ROOTS
    cat > fingerprints.json <<'SOURCE_FINGERPRINTS'
${fingerprints}
SOURCE_FINGERPRINTS
    python3 "\$(command -v reuse_source_validation.py)" validate \
        --sources sources.json --fingerprints fingerprints.json --output validated_sources.json
    """
}
