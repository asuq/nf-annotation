/*
 * Calculate in-house assembly statistics from staged FASTA files so ANI and
 * final reporting can reuse one authoritative TSV.
 */
process CALCULATE_ASSEMBLY_STATS {
    tag "assembly_stats"
    label 'process_medium'
    cache 'deep'
    publishDir(
        { "${params.outdir}/cohort/assembly_stats" },
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'versions.yml' ? null : filename },
    )

    input:
    path staged_manifest
    path staged_fastas

    output:
    path 'assembly_stats.tsv', emit: stats
    path 'versions.yml', emit: versions

    script:
    """
    calculate_assembly_stats.sh \
        --staged-manifest "${staged_manifest}" \
        --jobs "${task.cpus}" \
        --output assembly_stats.tsv

    seqtk_version="\$(command -v seqtk >/dev/null 2>&1 && seqtk 2>&1 | awk '/^Version:/ { print \$2; exit }' || true)"
    seqtk_version="\${seqtk_version:-NA}"
    printf '"%s":\n  seqtk: "%s"\n  script: "%s"\n' \
      "${task.process}" \
      "\${seqtk_version}" \
      'bin/calculate_assembly_stats.sh' \
      > versions.yml
    """

    stub:
    """
    printf 'accession\\tn50\\tscaffolds\\tgenome_size\\tgc_content\\n' > assembly_stats.tsv
    awk -F '\\t' 'NR > 1 { print \$1 "\\t80\\t1\\t80\\t50" }' '${staged_manifest}' >> assembly_stats.tsv
    cat <<'EOF' > versions.yml
    "${task.process}":
      seqtk: "stub"
      script: "bin/calculate_assembly_stats.sh"
    EOF
    """
}
