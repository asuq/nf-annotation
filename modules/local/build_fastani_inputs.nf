/*
 * Decide ANI eligibility, build canonical FastANI input files, and emit the
 * ANI-ready metadata table required by cluster_ani.py.
 */
process BUILD_FASTANI_INPUTS {
    tag "fastani_prep"
    label 'process_single'
    cache 'deep'
    publishDir(
        { "${params.outdir}/cohort/fastani" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            filename in ['ani_metadata.tsv', 'ani_exclusions.tsv', 'fastani_paths.txt']
                ? filename
                : null
        },
    )

    input:
    path validated_samples
    path metadata
    path staged_manifest
    path staged_fastas
    path checkm2_summaries
    path sixteen_s_statuses
    path busco_tables, name: 'busco_tables/busco_table??.tsv'
    path assembly_stats
    val primary_busco_column
    val ani_16s_policy

    output:
    path 'fastani_inputs', emit: fastani_inputs
    path 'fastani_paths.txt', emit: paths
    path 'ani_metadata.tsv', emit: metadata
    path 'ani_exclusions.tsv', emit: exclusions
    path 'versions.yml', emit: versions

    script:
    def buscoTableList = busco_tables instanceof Collection ? busco_tables : [busco_tables]
    def buscoArgs = buscoTableList.collect { "--busco \"${it}\"" }.join(' \\\n        ')
    """build_fastani_inputs.py \
    --validated-samples "${validated_samples}" \
    --metadata "${metadata}" \
    --staged-manifest "${staged_manifest}" \
    --checkm2 "${checkm2_summaries}" \
    --16s-status "${sixteen_s_statuses}" \
    ${buscoArgs} \
    --primary-busco-column "${primary_busco_column}" \
    --assembly-stats "${assembly_stats}" \
    --ani-16s-policy "${ani_16s_policy}" \
    --outdir .

    cat <<EOF > versions.yml
"${task.process}":
  python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
  script: "bin/build_fastani_inputs.py"
EOF
"""

}
