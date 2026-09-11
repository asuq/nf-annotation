/*
 * Validate the sample manifest and metadata table before any per-sample work.
 */
process VALIDATE_INPUTS {
    tag "${sample_csv.baseName}"
    label 'process_single'
    publishDir(
        "${params.outdir}/tables",
        mode: 'copy',
        overwrite: true,
        enabled: !params.update_from,
        saveAs: { filename ->
            filename in ['validated_samples.tsv', 'accession_map.tsv', 'validation_warnings.tsv']
                ? filename
                : null
        },
    )

    input:
    path sample_csv
    path metadata
    val busco_lineages

    output:
    path 'validated_samples.tsv', emit: validated_samples
    path 'accession_map.tsv', emit: accession_map
    path 'validation_warnings.tsv', emit: validation_warnings
    path 'sample_status.tsv', emit: sample_status
    path 'versions.yml', emit: versions

    script:
    def lineageArgs = (busco_lineages as List).collect {
        "--busco-lineage \"${it}\""
    }.join(' \\\n    ')
    """validate_inputs.py \
    --sample-csv "${sample_csv}" \
    --metadata "${metadata}" \
    ${lineageArgs} \
    --defer-genome-fasta-check \
    --genome-base-dir "${workflow.launchDir}" \
    --outdir .

cat <<EOF > versions.yml
"${task.process}":
  python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
  script: "bin/validate_inputs.py"
EOF
"""

}
