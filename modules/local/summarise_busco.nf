/*
 * Parse one BUSCO machine-readable summary into a stable single-row TSV.
 */
process SUMMARISE_BUSCO {
    tag "${meta.accession} / ${lineage}"
    label 'process_single'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/busco/${lineage}" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == "busco_summary_${lineage}.tsv" ? filename : null },
    )

    input:
    tuple val(meta), val(lineage), path(busco_summary_json)

    output:
    tuple val(meta), val(lineage), path("busco_summary_${lineage}.tsv"), emit: summary
    path 'versions.yml', emit: versions

    script:
    """summarise_busco.py \
    --accession "${meta.accession}" \
    --summary "${busco_summary_json}" \
    --lineage "${lineage}" \
    --output "busco_summary_${lineage}.tsv"

cat <<EOF > versions.yml
"${task.process}":
  python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
  script: "bin/summarise_busco.py"
EOF
"""

}
