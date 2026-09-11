/*
 * Assemble the final master table from the validated manifest, metadata block,
 * and the stable derived-summary tables collected earlier in the workflow.
 */
process BUILD_MASTER_TABLE {
    tag "master_table"
    label 'process_single'
    cache 'deep'
    errorStrategy 'finish'
    maxRetries 0

    input:
    path validated_samples
    path metadata
    val busco_lineages
    path taxonomy
    path checkm2
    path sixteen_s_status
    path busco_tables, name: 'busco_tables/busco_table??.tsv'
    path codetta_summary
    path ccfinder_strains
    path ani_summary
    path assembly_stats

    output:
    path 'master_table.tsv', emit: master_table
    path 'versions.yml', emit: versions

    script:
    def buscoTableList = busco_tables instanceof Collection ? busco_tables : [busco_tables]
    def buscoArgs = buscoTableList.collect { "--busco \"${it}\"" }.join(' \\\n        ')
    def lineageArgs = (busco_lineages as List).collect {
        "--busco-lineage \"${it}\""
    }.join(' \\\n        ')
    """
    build_master_table.py \
        --validated-samples "${validated_samples}" \
        --metadata "${metadata}" \
        ${lineageArgs} \
        --taxonomy "${taxonomy}" \
        --checkm2 "${checkm2}" \
        --16s-status "${sixteen_s_status}" \
        ${buscoArgs} \
        --codetta-summary "${codetta_summary}" \
        --ccfinder-strains "${ccfinder_strains}" \
        --ani "${ani_summary}" \
        --assembly-stats "${assembly_stats}" \
        --output master_table.tsv

    cat <<EOF > versions.yml
    "${task.process}":
      python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
      script: "bin/build_master_table.py"
    EOF
    """.stripIndent()

}
