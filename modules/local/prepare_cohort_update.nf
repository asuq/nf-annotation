/* Validate the entire reusable cohort before starting any new analyses. */
process PREPARE_COHORT_UPDATE {
    tag 'cohort_update'
    label 'process_single'
    cache 'deep'
    publishDir(
        "${params.outdir}/tables",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            filename in ['validated_samples.tsv', 'accession_map.tsv', 'validation_warnings.tsv',
                         'sample_status.tsv', 'cohort_update.tsv', 'cohort_update_run.json',
                         'inherited_versions.tsv'] ? filename : null
        },
    )

    input:
    path validated_samples, name: 'requested_samples.tsv'
    path accession_map, name: 'requested_accession_map.tsv'
    path initial_status, name: 'initial_status.tsv'
    path validation_warnings, name: 'requested_validation_warnings.tsv'
    path source_results, name: 'source_results'
    path genome_inputs, stageAs: 'candidate_genomes/genome??????'
    path metadata, name: 'metadata_input.tsv'
    path previous_update, name: 'previous_update.json'
    val busco_lineages
    val update_settings
    val destination

    output:
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
    def quote = { value -> "'" + value.toString().replace("'", "'\"'\"'") + "'" }
    def lineageArgs = busco_lineages.collect { "--busco-lineage ${quote.call(it)}" }.join(' ')
    def previousArg = previous_update ? '--previous-update previous_update.json' : ''
    def settingsJson = groovy.json.JsonOutput.toJson(update_settings)
    """
    cat <<'UPDATE_SETTINGS' > update_settings.json
    ${settingsJson}
    UPDATE_SETTINGS
    python3 "\$(command -v prepare_cohort_update.py)" \
        --validated-samples ${quote.call(validated_samples)} \
        --accession-map ${quote.call(accession_map)} \
        --initial-status initial_status.tsv \
        --validation-warnings ${quote.call(validation_warnings)} \
        --source-results source_results \
        --genome-inputs candidate_genomes \
        --metadata ${quote.call(metadata)} \
        --settings update_settings.json \
        --destination ${quote.call(destination)} \
        ${lineageArgs} ${previousArg} --outdir prepared
    cp prepared/* .
    printf '\"%s\":\\n  python: \"%s\"\\n  script: \"bin/prepare_cohort_update.py\"\\n' \
        '${task.process}' "\$(python3 --version 2>&1 | sed 's/^Python //')" > versions.yml
    """
}
