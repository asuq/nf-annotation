/*
 * Build and publish the cohort-level 16S FASTA and manifest outputs from the
 * per-sample 16S summaries generated during the QC branch.
 */
process BUILD_COHORT_16S {
    tag "cohort_16s"
    label 'process_single'
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
    script_path="\$(command -v concat_best_16s.py)"
    printf 'accession\tinternal_id\tstatus_tsv\tbest_16s_fasta\tqc_tsv\n' > cohort_inputs.tsv

    shopt -s nullglob
    for status_table in summaries/*_16S_status.tsv; do
        stem="\$(basename "\${status_table}" _16S_status.tsv)"
        internal_id="\${stem}"
        best_fasta="summaries/\${stem}_best_16S.fna"
        qc_table="qc/\${stem}_checkm2_summary.tsv"
        if [[ ! -f "\${best_fasta}" ]]; then
            printf 'Missing matching best_16S FASTA for %s\n' "\${status_table}" >&2
            exit 1
        fi
        if [[ ! -f "\${qc_table}" ]]; then
            printf 'Missing matching CheckM2 QC summary for %s\n' "\${status_table}" >&2
            exit 1
        fi
        accession="\$(awk -F '\t' 'NR == 2 { print \$1 }' "\${status_table}")"
        if [[ -z "\${accession}" ]]; then
            printf 'Missing accession in %s\n' "\${status_table}" >&2
            exit 1
        fi
        if [[ -z "\${internal_id}" ]]; then
            printf 'Missing internal ID in %s\n' "\${status_table}" >&2
            exit 1
        fi
        printf '%s\t%s\t%s\t%s\t%s\n' \
            "\${accession}" \
            "\${internal_id}" \
            "\$(pwd)/\${status_table}" \
            "\$(pwd)/\${best_fasta}" \
            "\$(pwd)/\${qc_table}" \
            >> cohort_inputs.tsv
    done

    python3 "\${script_path}" \
        --inputs cohort_inputs.tsv \
        --metadata "${metadata}" \
        --cohort-kind intact \
        --output-fasta all_best_16S.fna \
        --output-manifest all_best_16S_manifest.tsv

    python3 "\${script_path}" \
        --inputs cohort_inputs.tsv \
        --metadata "${metadata}" \
        --cohort-kind partial \
        --output-fasta all_partial_16S.fna \
        --output-manifest all_partial_16S_manifest.tsv

    python3 "\${script_path}" \
        --inputs cohort_inputs.tsv \
        --metadata "${metadata}" \
        --cohort-kind low_quality_intact \
        --output-fasta low_quality_best_16S.fna \
        --output-manifest low_quality_best_16S_manifest.tsv

    python3 "\${script_path}" \
        --inputs cohort_inputs.tsv \
        --metadata "${metadata}" \
        --cohort-kind low_quality_partial \
        --output-fasta low_quality_partial_16S.fna \
        --output-manifest low_quality_partial_16S_manifest.tsv

    python_version="\$(python3 --version 2>&1 | sed 's/^Python //')"
    printf '"%s":\n  python: "%s"\n  script: "%s"\n' \
        "${task.process}" \
        "\${python_version}" \
        "bin/concat_best_16s.py" \
        > versions.yml
    """

}
