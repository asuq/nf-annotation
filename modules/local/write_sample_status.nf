/*
 * Build the authoritative sample-status table by overlaying downstream module
 * outcomes onto the validation-time status seed.
 */
process WRITE_SAMPLE_STATUS {
    tag "sample_status"
    label 'process_single'
    cache 'deep'
    errorStrategy 'finish'
    maxRetries 0

    input:
    path validated_samples
    path initial_status, name: 'initial_status.tsv'
    val busco_lineages
    path metadata
    path taxonomy
    path checkm2
    path sixteen_s_status
    path busco_tables, name: 'busco_tables/busco_table??.tsv'
    path codetta_summary
    path ccfinder_strains
    path prokka_manifest
    path ani_summary
    path assembly_stats
    val primary_busco_column
    val ani_allow_incomplete_16s

    output:
    path 'sample_status.tsv', emit: sample_status
    path 'versions.yml', emit: versions

    script:
    def buscoTableList = busco_tables instanceof Collection ? busco_tables : [busco_tables]
    def buscoArgs = buscoTableList.collect { "--busco \"${it}\"" }.join(' \\\n        ')
    def lineageArgs = (busco_lineages as List).collect {
        "--busco-lineage \"${it}\""
    }.join(' \\\n        ')
    def allowIncomplete16sArg = ani_allow_incomplete_16s ? '--ani-allow-incomplete-16s' : ''
    """
    script_path="\$(command -v build_sample_status.py)"
    python3 "\${script_path}" \
        --validated-samples "${validated_samples}" \
        --initial-status "${initial_status}" \
        ${lineageArgs} \
        --metadata "${metadata}" \
        --taxonomy "${taxonomy}" \
        --checkm2 "${checkm2}" \
        --16s-status "${sixteen_s_status}" \
        ${buscoArgs} \
        --codetta-summary "${codetta_summary}" \
        --ccfinder-strains "${ccfinder_strains}" \
        --prokka-manifest "${prokka_manifest}" \
        --ani "${ani_summary}" \
        --assembly-stats "${assembly_stats}" \
        --primary-busco-column "${primary_busco_column}" \
        ${allowIncomplete16sArg} \
        --output sample_status.tsv

    cat <<EOF > versions.yml
    "${task.process}":
      python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
      script: "bin/build_sample_status.py"
    EOF
    """.stripIndent()

}
