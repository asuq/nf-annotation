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
    path 'versions.yml', emit: versions

    script:
    def aniThreshold = params.ani_threshold
    """
    cluster_ani.py \
        --ani-matrix "${ani_matrix}" \
        --ani-metadata "${ani_metadata}" \
        --matrix-name-column matrix_name \
        --threshold ${aniThreshold} \
        --threads ${task.cpus} \
        --outdir .

    cat <<EOF > versions.yml
    "${task.process}":
      python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
      script: "bin/cluster_ani.py"
    EOF
    """

    stub:
    '''
    cat <<'EOF' > cluster.tsv
    Accession	Cluster_ID	Matrix_Name
    sample_a	C000001	fastani_inputs/sample_a.fasta
    EOF
    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      script: "bin/cluster_ani.py"
    EOF
    '''
}
