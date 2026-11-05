/*
 * Run ANI complete-linkage clustering and emit stable cluster memberships.
 */
process CLUSTER_ANI {
    tag "ani_cluster"
    label 'process_medium'
    publishDir(
        { "${params.outdir}/cohort/ani_clusters" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
    )

    input:
    path ani_matrix
    path ani_metadata

    output:
    path 'cluster.tsv', emit: clusters
    path 'ani_summary.tsv', emit: summary
    path 'ani_representatives.tsv', emit: representatives
    path 'versions.yml', emit: versions

    script:
    def aniThreshold = params.ani_threshold
    """
    python3 "\$(command -v analyse_ani.py)" \
        --ani-matrix "${ani_matrix}" \
        --ani-metadata "${ani_metadata}" \
        --matrix-name-column matrix_name \
        --threshold ${aniThreshold} \
        --ani-score-profile "${params.score_profile}" \
        --threads ${task.cpus} \
        --outdir .

    cat <<EOF > versions.yml
    "${task.process}":
      python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
      script: "bin/analyse_ani.py"
    EOF
    """

}
