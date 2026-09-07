/*
 * Merge per-process versions files with runtime/container/resource context into
 * one final TSV report.
 */
process COLLECT_VERSIONS {
    tag "versions"
    label 'process_single'
    cache 'deep'
    publishDir(
        "${params.outdir}/tables",
        mode: 'copy',
        overwrite: true,
        saveAs: { filename -> filename == 'tool_and_db_versions.tsv' ? filename : null },
    )

    input:
    path version_files, stageAs: 'version_files/versions??.yml'
    path inherited_versions, stageAs: 'inherited_versions/report??.tsv'
    val busco_lineages
    val nextflow_version
    val pipeline_version
    val git_commit
    val container_engine

    output:
    path 'tool_and_db_versions.tsv', emit: versions_table

    script:
    def lineageArgs = (busco_lineages as List).collect {
        "--busco-lineage \"${it}\""
    }.join(' \\\n        ')
    def containerRefs = [
        'python'  : params.python_container ?: 'NA',
        'seqtk'   : params.seqtk_container ?: 'NA',
        'barrnap' : params.barrnap_container ?: 'NA',
        'checkm2' : params.checkm2_container ?: 'NA',
        'busco'   : params.busco_container ?: 'NA',
        'prokka'  : params.prokka_container ?: 'NA',
        'ccfinder': params.ccfinder_container ?: 'NA',
        'codetta' : params.codetta_container ?: 'NA',
        'padloc'  : params.padloc_container ?: 'NA',
        'eggnog'  : params.eggnog_container ?: 'NA',
        'fastani' : params.fastani_container ?: 'NA',
    ]
    def containerArgs = containerRefs.collect { name, value ->
        "--container-ref \"${name}=${value}\""
    }.join(' \\\n        ')
    def inheritedArgs = inherited_versions ? '--inherited-version-dir inherited_versions' : ''
    """
    script_path="\$(command -v collect_versions.py)"
    python3 "\${script_path}" \
        --version-dir version_files \
        ${inheritedArgs} \
        --nextflow-version "${nextflow_version}" \
        --pipeline-version "${pipeline_version}" \
        --git-commit "${git_commit}" \
        --container-engine "${container_engine}" \
        --checkm2-db "${params.checkm2_db ?: 'NA'}" \
        --checkm2-db-label "${params.checkm2_db_label ?: 'NA'}" \
        --taxdump "${params.taxdump ?: 'NA'}" \
        --taxdump-label "${params.taxdump_label ?: 'NA'}" \
        --codetta-db "${params.codetta_db ?: 'NA'}" \
        --codetta-db-label "${params.codetta_db_label ?: 'NA'}" \
        ${lineageArgs} \
        --busco-db "${params.busco_db ?: 'NA'}" \
        --eggnog-db "${params.eggnog_db ?: 'NA'}" \
        --eggnog-db-label "${params.eggnog_db_label ?: 'NA'}" \
        ${containerArgs} \
        --output tool_and_db_versions.tsv
    """

    stub:
    """
    cat <<'EOF' > tool_and_db_versions.tsv
    component	kind	version	image_or_path	notes
    nextflow	runtime	stub	NA	workflow
    python	runtime	stub	NA	reported by VALIDATE_INPUTS
    codetta	tool	v2.0	NA	reported by CODETTA
    codetta_source_commit	tool	863359ed326276602d44e48227b6003ac6ffd266	NA	reported by CODETTA
    codetta	container	NA	${params.codetta_container ?: 'NA'}	params container reference
    codetta_db	database	${params.codetta_db_label ?: 'NA'}	${params.codetta_db ?: 'NA'}	Codetta profile database
    EOF
    """
}
