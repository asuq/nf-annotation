/*
 * Stage each genome to an internal-ID-safe FASTA name for downstream tools.
 */
process STAGE_INPUTS {
    tag "${meta.accession}"
    label 'process_single'
    cache 'deep'
    publishDir(
        { "${params.outdir}/samples/${meta.accession}/staged" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
    )

    input:
    tuple val(meta), path(genome, name: 'input_genome')

    output:
    tuple val(meta), path("${meta.internal_id}.fasta"), emit: staged_fasta
    tuple val(meta), path("${meta.internal_id}.fasta.fai"), emit: fai, optional: true
    path 'versions.yml', emit: versions

    script:
    def internalId = meta.internal_id ?: ''
    """
    if [[ -z "${internalId}" ]]; then
        echo "meta.internal_id is required for STAGE_INPUTS." >&2
        exit 1
    fi

    seqtk seq -A "${genome}" > "${internalId}.fasta"

    if command -v samtools >/dev/null 2>&1; then
        samtools faidx "${internalId}.fasta"
    fi

    seqtk_version="\$(command -v seqtk >/dev/null 2>&1 && seqtk 2>&1 | awk '/^Version:/ { print \$2; exit }' || true)"
    seqtk_version="\${seqtk_version:-NA}"
    samtools_version="\$(command -v samtools >/dev/null 2>&1 && samtools --version 2>&1 | head -n 1 || echo NA)"
    printf '"%s":\n  seqtk: "%s"\n  samtools: "%s"\n' \
      "${task.process}" \
      "\${seqtk_version}" \
      "\${samtools_version}" \
      > versions.yml
    """

    stub:
    """
    cat <<'EOF' > "${meta.internal_id}.fasta"
    >${meta.internal_id}
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    EOF
    : > "${meta.internal_id}.fasta.fai"
    cat <<'EOF' > versions.yml
    "${task.process}":
      seqtk: "stub"
      samtools: "stub"
    EOF
    """
}
