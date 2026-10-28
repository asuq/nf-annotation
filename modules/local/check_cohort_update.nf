/* Always check resume identity, without invalidating cached source validation. */
process CHECK_COHORT_UPDATE {
    label 'process_single'
    container params.python_container
    cache false
    errorStrategy 'finish'
    maxRetries 0
    publishDir(
        "${params.outdir}/tables",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            filename in ['validated_samples.tsv', 'accession_map.tsv', 'validation_warnings.tsv',
                         'cohort_update.tsv', 'cohort_update_run.json',
                         'inherited_versions.tsv', 'source_validation.json'] ? filename : null
        },
    )

    input:
    path prepared
    path previous_update, name: 'previous_update.json'

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
    path 'source_validation.json', emit: source_validation

    script:
    def previousArg = previous_update ? '--previous-update previous_update.json' : ''
    """
    python3 "\$(command -v check_cohort_update.py)" \
        --prepared '${prepared}' ${previousArg} --output .
    """
}
