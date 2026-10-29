/* Validate the entire reusable cohort before starting any new analyses. */
process PREPARE_COHORT_UPDATE {
    tag 'cohort_update'
    label 'process_cohort'
    stageInMode 'symlink'
    cache 'deep'
    publishDir "${params.outdir}/pipeline_info/cohort_validation", mode: 'copy', overwrite: true, pattern: 'cohort_validation_progress.json'

    input:
    path validated_samples, name: 'requested_samples.tsv'
    path accession_map, name: 'requested_accession_map.tsv'
    path initial_status, name: 'initial_status.tsv'
    path validation_warnings, name: 'requested_validation_warnings.tsv'
    path source_results, stageAs: 'source_results/source????'
    path genome_inputs, stageAs: 'candidate_genomes/genome??????'
    path metadata, name: 'metadata_input.tsv'
    val busco_lineages
    val update_settings
    val destination
    val policy_code
    val source_validation
    path dependencies, stageAs: 'dependencies/root????'

    output:
    path 'cohort_validation_progress.json', emit: progress
    path 'prepared', emit: prepared
    path 'validated_samples.tsv', emit: validated_samples
    path 'accession_map.tsv', emit: accession_map
    path 'validation_warnings.tsv', emit: validation_warnings
    path 'sample_status.tsv', emit: sample_status
    path 'new_samples.tsv', emit: new_samples
    path 'reused_samples.tsv', emit: reused_samples
    path 'cohort_update.tsv', emit: audit
    path 'cohort_update_run.json', emit: identity
    path 'inherited_versions.tsv', emit: inherited_versions
    path 'versions.yml', emit: versions

    script:
    def workers = Math.min(Math.min(task.cpus as int, params.cohort_validation_workers as int), 8)
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def lineageArgs = busco_lineages.collect { "--busco-lineage ${quote.call(it)}" }.join(' ')
    def sourceValidationJson = groovy.json.JsonOutput.toJson(source_validation)
    def settingsJson = groovy.json.JsonOutput.toJson(update_settings)
    def sourceArgs = PublishedSources.arguments(source_results, '--source-results')
    """
    # sample_policy_sha256=${policy_code}

    cat <<'SOURCE_VALIDATION' > source_validation.json
${sourceValidationJson}
SOURCE_VALIDATION
    cat <<'UPDATE_SETTINGS' > update_settings.json
${settingsJson}
UPDATE_SETTINGS
    python3 "\$(command -v prepare_cohort_update.py)" \
        --validated-samples ${quote.call(validated_samples)} \
        --accession-map ${quote.call(accession_map)} \
        --initial-status initial_status.tsv \
        --validation-warnings ${quote.call(validation_warnings)} \
        ${sourceArgs} \
        --genome-inputs candidate_genomes \
        --metadata ${quote.call(metadata)} \
        --settings update_settings.json --source-validation source_validation.json \
        --gcode-rule ${quote.call(params.gcode_rule)} \
        --destination ${quote.call(destination)} \
        ${lineageArgs} --outdir prepared --workers ${workers} --progress cohort_validation_progress.json
    cp source_validation.json prepared/source_validation.json
    cp prepared/* .
    printf '\"%s\":\\n  python: \"%s\"\\n  script: \"bin/prepare_cohort_update.py\"\\n' \
        '${task.process}' "\$(python3 --version 2>&1 | sed 's/^Python //')" > versions.yml
    """
}
