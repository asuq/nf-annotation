/* Preserve a validated, published input for annotation and later reannotation. */
process PREPARE_ANNOTATION_BUNDLE {
    tag "${meta.accession}"
    label 'process_single'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/annotation" },
        mode: 'copy', overwrite: true,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
    )

    input:
    tuple val(meta), path(genome), val(gcode), path(faa), path(gff), path(gbk), path(prokka_log)

    output:
    tuple val(meta), path('bundle'), emit: bundle
    path 'versions.yml', emit: versions

    script:
    """
    python3 "\$(command -v prepare_annotation_bundle.py)" \
        --accession '${meta.accession}' --gcode '${gcode}' \
        --genome '${genome}' --faa '${faa}' --gff '${gff}' --gbk '${gbk}' \
        --prokka-log '${prokka_log}' --outdir bundle
    python_version="\$(python3 --version 2>&1 | sed 's/^Python //')"
    printf '"%s":\\n  python: "%s"\\n  script: "bin/prepare_annotation_bundle.py"\\n' \
        '${task.process}' "\${python_version}" > versions.yml
    """
}
