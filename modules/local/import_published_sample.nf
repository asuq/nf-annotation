/* Copy published artefacts and restore unique execution filenames for joins. */
process IMPORT_PUBLISHED_SAMPLE {
    tag "${meta.accession}"
    label 'process_single'
    cache 'deep'
    publishDir(
        "${params.outdir}/samples",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename.startsWith('sample/') ? "${meta.accession}/${filename.substring(7)}" : null },
    )

    input:
    tuple val(meta), path(source_sample, name: 'source_sample'), val(gcode)
    val busco_lineages

    output:
    // Each retained upstream entry owns a separate publication destination;
    // none can replace the parent of a rebuilt annotation bundle or result.
    tuple val(meta), path('sample/*', arity: '1..*', hidden: true), path('channels'), val(gcode), val(busco_lineages), emit: results
    path 'versions.yml', emit: versions

    script:
    """
    mkdir sample channels
    python3 - <<'PY'
    import shutil
    from pathlib import Path
    for source in Path('source_sample').iterdir():
        if source.name != 'annotation':
            destination = Path('sample') / source.name
            if source.is_dir():
                shutil.copytree(source, destination)
            else:
                shutil.copyfile(source, destination)
    PY
    cp sample/16s/best_16S.fna channels/${meta.internal_id}_best_16S.fna
    cp sample/16s/16S_status.tsv channels/${meta.internal_id}_16S_status.tsv
    cp sample/checkm2/checkm2_summary.tsv channels/${meta.internal_id}_checkm2_summary.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      transform: "copy_validated_published_sample"
    EOF
    """
}
