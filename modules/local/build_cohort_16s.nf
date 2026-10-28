/*
 * Build and publish the cohort-level 16S FASTA and manifest outputs from the
 * per-sample 16S summaries generated during the QC branch.
 */
process BUILD_COHORT_16S {
    tag "cohort_16s"
    label 'process_cohort'
    cache 'deep'
    publishDir(
        "${params.outdir}/cohort/16s",
        mode: 'copy',
        overwrite: true,
    )

    input:
    path best_fastas, stageAs: 'summaries/*'
    path status_tables, stageAs: 'summaries/*'
    path qc_tables, stageAs: 'qc/*'
    path metadata

    output:
    path 'all_best_16S.fna', emit: best_fasta
    path 'all_best_16S_manifest.tsv', emit: best_manifest
    path 'all_partial_16S.fna', emit: partial_fasta
    path 'all_partial_16S_manifest.tsv', emit: partial_manifest
    path 'low_quality_best_16S.fna', emit: low_quality_best_fasta
    path 'low_quality_best_16S_manifest.tsv', emit: low_quality_best_manifest
    path 'low_quality_partial_16S.fna', emit: low_quality_partial_fasta
    path 'low_quality_partial_16S_manifest.tsv', emit: low_quality_partial_manifest
    path 'versions.yml', emit: versions

    script:
    """
    shopt -s nullglob
    printf '%s\\0' summaries/*_16S_status.tsv > status_paths.list0
    python3 "\$(command -v build_all_cohort_16s.py)" \
        --statuses status_paths.list0 --qc-dir qc --metadata "${metadata}" --output .

    python_version="\$(python3 --version 2>&1 | sed 's/^Python //')"
    printf '"%s":\n  python: "%s"\n  script: "%s"\n' \
        "${task.process}" \
        "\${python_version}" \
        "bin/build_all_cohort_16s.py" \
        > versions.yml
    """

}
