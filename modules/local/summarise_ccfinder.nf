/*
 * Summarise one CRISPRCasFinder JSON output into strain, contig, and CRISPR
 * TSVs. The Python CLI emits a failed sample summary when the JSON is invalid.
 */
process SUMMARISE_CCFINDER {
    tag "${meta.accession}"
    label 'process_single'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/ccfinder" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
    )

    input:
    tuple val(meta), path(result_json)

    output:
    tuple val(meta), path('ccfinder_strains.tsv'), path('ccfinder_contigs.tsv'), path('ccfinder_crisprs.tsv'), emit: summaries
    tuple val(meta), path('ccfinder_strains.tsv'), emit: strains
    path 'versions.yml', emit: versions

    script:
    """
    summarise_ccfinder.py \
        --accession "${meta.accession}" \
        --result-json "${result_json}" \
        --outdir .

    cat <<EOF > versions.yml
    "${task.process}":
      python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
      script: "bin/summarise_ccfinder.py"
    EOF
    """

}
