/* Copy a validated source independently of old work directories and caches. */
process IMPORT_ANNOTATION_SOURCE {
    label 'process_single'
    container params.python_container
    cache false
    errorStrategy 'finish'
    maxRetries 0
    publishDir params.outdir, mode: 'copy', overwrite: true,
        saveAs: { name -> name.startsWith('imported/inherited/') ? name.substring(19) : null }

    input:
    path source_results, name: 'source_results'
    path preflight

    output:
    path 'imported/upstream_master.tsv', emit: master
    path 'imported/upstream_sample_status.tsv', emit: sample_status
    path 'imported/source_samples.tsv', emit: samples
    path 'imported/inherited/samples', emit: sample_tree
    path 'imported/inherited/tables', emit: inherited_tables

    script:
    """
    python3 "\$(command -v prepare_annotation_tasks.py)" import \
        --source source_results --output imported
    """
}
