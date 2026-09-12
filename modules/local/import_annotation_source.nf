/* Copy a validated source independently of old work directories and caches. */
process IMPORT_ANNOTATION_SOURCE {
    label 'process_single'
    container params.python_container
    cache false
    errorStrategy 'finish'
    maxRetries 0
    publishDir params.outdir, mode: 'copy', overwrite: true,
        saveAs: { name -> name.startsWith('imported/inherited/samples/') || name.startsWith('imported/inherited/tables/') ? name.substring(19) : null }

    input:
    path source_results, name: 'source_results'
    path preflight

    output:
    path 'imported/upstream_master.tsv', emit: master
    path 'imported/upstream_sample_status.tsv', emit: sample_status
    path 'imported/source_samples.tsv', emit: samples
    // Only publish files: copying a parent directory could replace concurrently
    // published functional results or regenerated cohort tables on resume.
    path 'imported/inherited/**', type: 'file', hidden: true, emit: published_files

    script:
    """
    python3 "\$(command -v prepare_annotation_tasks.py)" import \
        --source source_results --output imported
    """
}
