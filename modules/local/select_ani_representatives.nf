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

    stub:
    '''
    cat <<'EOF' > ani_summary.tsv
    Accession	Cluster_ID	Is_Representative	ANI_to_Representative	Score
    sample_a	C000001	yes	100.0000	7.500000
    EOF
    cat <<'EOF' > ani_representatives.tsv
    Cluster_ID	Representative_Accession	Organism_Name	CheckM2_Completeness	CheckM2_Contamination	BUSCO	Assembly_Level	N50	Cluster_Size
    C000001	sample_a	Sample_A	95.00	1.00	C:98.0%[S:98.0%,D:0.0%],F:1.0%,M:1.0%,n:200	Scaffold	50000	1
    EOF
    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      script: "bin/select_ani_representatives.py"
    EOF
    '''
}
