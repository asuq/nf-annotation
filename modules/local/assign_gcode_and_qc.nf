/*
 * Merge paired CheckM2 reports into a single per-sample QC summary.
 */
process ASSIGN_GCODE_AND_QC {
    tag "${meta.accession}"
    label 'process_single'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/checkm2" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'checkm2_summary.tsv' ? filename : null },
    )

    input:
    tuple val(meta), path(checkm2_gcode4_report, name: 'checkm2_gcode4_report.tsv'), path(checkm2_gcode11_report, name: 'checkm2_gcode11_report.tsv')

    output:
    tuple val(meta), path('checkm2_summary.tsv'), emit: summary
    tuple val(meta), path("${meta.internal_id}_checkm2_summary.tsv"), emit: cohort_summary
    path 'versions.yml', emit: versions

    script:
    """
    summarise_checkm2.py \
        --accession "${meta.accession}" \
        --gcode4-report "${checkm2_gcode4_report}" \
        --gcode11-report "${checkm2_gcode11_report}" \
        --output checkm2_summary.tsv

    cp checkm2_summary.tsv "${meta.internal_id}_checkm2_summary.tsv"

    cat <<EOF > versions.yml
    "${task.process}":
      python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
      script: "bin/summarise_checkm2.py"
    EOF
    """

    stub:
    """
    summarise_checkm2.py \
        --accession "${meta.accession}" \
        --gcode4-report "${checkm2_gcode4_report}" \
        --gcode11-report "${checkm2_gcode11_report}" \
        --output checkm2_summary.tsv
    cp checkm2_summary.tsv "${meta.internal_id}_checkm2_summary.tsv"
    cat <<'EOF' > versions.yml
    "${task.process}":
      python: "stub"
      script: "bin/summarise_checkm2.py"
    EOF
    """
}
