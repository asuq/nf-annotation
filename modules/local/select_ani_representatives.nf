/*
 * Select one representative per ANI cluster and emit both the published
 * representative table and the accession-keyed ANI summary for final joins.
 */
process SELECT_ANI_REPRESENTATIVES {
    tag "ani_reps"
    label 'process_medium'
    publishDir(
        "${params.outdir}/cohort/ani_clusters",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            filename in ['ani_representatives.tsv', 'ani_summary.tsv']
                ? filename
                : null
        },
    )

    input:
    path ani_clusters
    path ani_metadata
    path ani_matrix

    output:
    path 'ani_summary.tsv', emit: ani_summary
    path 'ani_representatives.tsv', emit: ani_representatives
    path 'versions.yml', emit: versions

    script:
    def aniScoreProfile = params.ani_score_profile
    """
    select_ani_representatives.py \
        --ani-clusters "${ani_clusters}" \
        --ani-metadata "${ani_metadata}" \
        --ani-matrix "${ani_matrix}" \
        --ani-score-profile "${aniScoreProfile}" \
        --ani-summary-output ani_summary.tsv \
        --ani-representatives-output ani_representatives.tsv

    cat <<EOF > versions.yml
    "${task.process}":
      python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
      script: "bin/select_ani_representatives.py"
    EOF
    """

}
