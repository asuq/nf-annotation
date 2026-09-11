/*
 * Summarise one Codetta output directory into the stable reporting row.
 */
process SUMMARISE_CODETTA {
    tag "${meta.accession}"
    label 'process_single'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/codetta" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
    )

    input:
    tuple val(meta), path(codetta_dir)

    output:
    tuple val(meta), path('codetta_summary.tsv'), emit: summary
    path 'versions.yml', emit: versions

    script:
    """
    script_path="\$(command -v summarise_codetta.py)"
    python3 "\${script_path}" \
        --accession "${meta.accession}" \
        --codetta-dir "${codetta_dir}" \
        --output codetta_summary.tsv

    python_version="\$(python3 --version 2>&1 | sed 's/^Python //')"
    printf '%s\n' \
        '"${task.process}":' \
        "  python: \"\${python_version}\"" \
        '  script: "bin/summarise_codetta.py"' \
        > versions.yml
    """

}
