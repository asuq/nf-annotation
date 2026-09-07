/* Copy published artefacts and restore unique execution filenames for joins. */
process IMPORT_PUBLISHED_SAMPLE {
    tag "${meta.accession}"
    label 'process_single'
    cache 'deep'
    publishDir(
        "${params.outdir}/samples",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'sample' ? meta.accession : null },
    )

    input:
    tuple val(meta), path(source_sample, name: 'source_sample'), val(gcode), val(eggnog_status)
    val busco_lineages

    output:
    tuple val(meta), path('sample'), path('channels'), val(gcode), val(eggnog_status), val(busco_lineages), emit: results
    path 'versions.yml', emit: versions

    script:
    """
    mkdir sample channels
    cp -R source_sample/. sample/
    cp sample/16s/best_16S.fna channels/${meta.internal_id}_best_16S.fna
    cp sample/16s/16S_status.tsv channels/${meta.internal_id}_16S_status.tsv
    cp sample/checkm2/checkm2_summary.tsv channels/${meta.internal_id}_checkm2_summary.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      transform: "copy_validated_published_sample"
    EOF
    """
}
