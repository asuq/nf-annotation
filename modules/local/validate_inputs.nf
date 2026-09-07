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
            filename in ['validated_samples.tsv', 'accession_map.tsv', 'validation_warnings.tsv', 'sample_status.tsv']
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
    --outdir .

cat <<EOF > versions.yml
"${task.process}":
  python: "\$(python3 --version 2>&1 | sed 's/^Python //')"
  script: "bin/validate_inputs.py"
EOF
"""

    stub:
    def stubGenome = file("${projectDir}/assets/fixtures/stub/genomes/TEST_ACC.fasta").toString()
    def sampleStatusColumns = [
        'accession',
        'internal_id',
        'is_new',
        'validation_status',
        'taxonomy_status',
        'barrnap_status',
        'checkm2_gcode4_status',
        'checkm2_gcode11_status',
        'gcode_status',
        'gcode',
        'low_quality',
    ] + (busco_lineages as List).collect { "busco_${it}_status" } + [
        'codetta_status',
        'prokka_status',
        'ccfinder_status',
        'padloc_status',
        'eggnog_status',
        'ani_included',
        'ani_exclusion_reason',
        'warnings',
        'notes',
    ]
    def sampleStatusStubValues = [
        accession: 'TEST_ACC',
        internal_id: 'TEST_ACC',
        is_new: 'false',
        validation_status: 'done',
        warnings: 'stub_warning',
        notes: 'stub warning',
        gcode: 'NA',
        low_quality: 'NA',
        ani_included: 'na',
    ]
    def sampleStatusRow = sampleStatusColumns.collect { column ->
        if (sampleStatusStubValues.containsKey(column)) {
            return sampleStatusStubValues[column]
        }
        return column.endsWith('_status') ? 'na' : ''
    }.join('\t')
    """cat <<'EOF' > validated_samples.tsv
accession	is_new	assembly_level	genome_fasta	internal_id
TEST_ACC	false	NA	${stubGenome}	TEST_ACC
EOF
cat <<'EOF' > accession_map.tsv
accession	internal_id	is_new	assembly_level	genome_fasta	metadata_present
TEST_ACC	TEST_ACC	false	NA	${stubGenome}	true
EOF
cat <<'EOF' > validation_warnings.tsv
accession	warning_code	message
EOF
cat <<'EOF' > sample_status.tsv
${sampleStatusColumns.join('\t')}
${sampleStatusRow}
EOF
cat <<'EOF' > versions.yml
"${task.process}":
  python: "stub"
  script: "bin/validate_inputs.py"
EOF
"""
}
