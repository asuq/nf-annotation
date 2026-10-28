process COMPLETE_ANI_REFERENCES {
    label 'process_cohort'
    container params.python_container
    stageInMode 'symlink'
    cache false
    errorStrategy 'finish'
    maxRetries 0

    input:
    path audits, stageAs: 'audits/audit????.json'
    path samples
    path destination, name: 'destination'
    path dependencies, stageAs: 'dependencies/root????'

    output:
    path 'reference_complete.json', emit: complete

    script:
    def paths = groovy.json.JsonOutput.toJson((audits instanceof Collection ? audits : [audits]).collect { it.toString() })
    """
    cat > audits.json <<'REFERENCE_AUDITS'
${paths}
REFERENCE_AUDITS
    python3 "\$(command -v complete_ani_references.py)" --audits audits.json \
        --samples '${samples}' --destination destination --receipt reference_complete.json
    """
}
