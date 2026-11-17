process ASSEMBLY_ELIGIBILITY {
    label 'process_single'
    container params.python_container

    input:
    path samples
    path stats

    output:
    path 'post_qc_samples.tsv', emit: samples

    script:
    """
    assembly_eligibility.py "${samples}" "${stats}" post_qc_samples.tsv
    """
}
