/* Bounded bookkeeping allocations; publish each sample once on shared storage. */
process IMPORT_ANI_BATCH {
    tag "batch ${batch_id} (${records.size()} samples)"
    label 'process_cohort'
    stageInMode 'symlink'
    container params.python_container
    cache false
    maxForks params.ani_import_max_forks
    errorStrategy 'finish'
    maxRetries 0

    input:
    tuple val(batch_id), val(records)
    path source_samples, stageAs: 'source_roots/source????'
    path dependencies, stageAs: 'dependencies/root????'
    path destination, name: 'destination'
    val destination_root
    val source_validation

    output:
    tuple path('import_audit.json'), path('channels'), val(destination_root), emit: results
    path 'import_audit.json', emit: audit
    path 'versions.yml', emit: versions

    script:
    def validationJson = groovy.json.JsonOutput.toJson(source_validation)
    def recordJson = groovy.json.JsonOutput.toJson(records)
    def sourceJson = groovy.json.JsonOutput.toJson((source_samples instanceof Collection ? source_samples : [source_samples]).collect { it.toString() })
    """
    cat > records.json <<'ANI_IMPORT_RECORDS'
${recordJson}
ANI_IMPORT_RECORDS
    cat > sources.json <<'ANI_IMPORT_SOURCES'
${sourceJson}
ANI_IMPORT_SOURCES
    if [[ ! -L destination ]]; then
        echo 'ANI publication roots require link staging on shared storage; do not use stageInMode copy.' >&2
        exit 1
    fi
    cat > source_validation.json <<'SOURCE_VALIDATION'
${validationJson}
SOURCE_VALIDATION
    import_status=0
    python3 "\$(command -v import_ani_batch.py)" \
        --records records.json --sources sources.json --destination destination \
        --mode '${params.ani_publication_mode}' --source-validation source_validation.json \
        --channels channels --audit import_audit.json || import_status=\$?
    mkdir -p destination/pipeline_info/ani_import
    if [[ -f import_audit.json ]]; then
        cp import_audit.json destination/pipeline_info/ani_import/batch_${batch_id}.json
    fi
    if [[ "\${import_status}" != 0 ]]; then exit "\${import_status}"; fi
    cat <<'VERSIONS' > versions.yml
"${task.process}":
  transform: "batched_validated_ani_reuse"
VERSIONS
    """
}
