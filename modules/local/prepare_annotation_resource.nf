/* Prepare one immutable annotation resource in its matching tool runtime. */
process PREP_ANNOTATION_RESOURCE {
    tag "${component}"
    label 'process_medium'
    label 'prep_annotation_resource'
    container { params["${component}_container"] }
    stageInMode 'symlink'

    input:
    tuple val(component), val(destination), path(destination_parent), val(destination_name), val(download_enabled), val(version), path(source_manifest)

    output:
    path 'annotation_resource_report.tsv', emit: report
    path 'versions.yml', emit: versions

    script:
    """
    python "\$(command -v prepare_runtime_databases.py)" \
        --${component}-dest "${destination_parent}/${destination_name}" \
        ${version ? "--${component}-version '${version}'" : ''} \
        --remote-source-manifest "${source_manifest}" \
        ${download_enabled ? '--download' : ''} \
        --report annotation_resource_report.tsv
    printf '%s\\n' '"${task.process}":' '  helper: "annotation resource schema 1"' > versions.yml
    """

    stub:
    """
    printf 'component\\tstatus\\tsource\\tdestination\\tdetails\\n' > annotation_resource_report.tsv
    printf '%s\\tstub\\tstub\\t%s\\tqualification=pending\\n' '${component}' '${destination}' >> annotation_resource_report.tsv
    printf '%s\\n' '"${task.process}":' '  helper: "stub"' > versions.yml
    """
}
