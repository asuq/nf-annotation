/*
 * Parse Barrnap raw outputs in the shared Python helper image and emit the
 * per-sample 16S summary artefacts used by downstream cohort and reporting
 * steps.
 */
process SUMMARISE_16S {
    tag "${meta.accession}"
    label 'process_single'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/16s" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename ->
            if (filename == "${meta.internal_id}_best_16S.fna") {
                return 'best_16S.fna'
            }
            if (filename == "${meta.internal_id}_16S_status.tsv") {
                return '16S_status.tsv'
            }
            return null
        },
    )

    input:
    tuple val(meta), path(rrna_gff), path(rrna_fasta), path(barrnap_log)

    output:
    tuple val(meta), path("${meta.internal_id}_best_16S.fna"), path("${meta.internal_id}_16S_status.tsv"), emit: summaries
    path 'versions.yml', emit: versions

    script:
    """
    summarise_16s.py \
        --accession "${meta.accession}" \
        --rrna-gff "${rrna_gff}" \
        --rrna-fasta "${rrna_fasta}" \
        --outdir .

    cp best_16S.fna "${meta.internal_id}_best_16S.fna"
    cp 16S_status.tsv "${meta.internal_id}_16S_status.tsv"

    cat <<EOF > versions.yml
    "${task.process}":
      python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
      summarise_16s: "bin/summarise_16s.py"
    EOF
    """

}
