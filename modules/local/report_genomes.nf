/* Genome reporting when functional searches are disabled, keyed to the exact table and bundle manifests. */
process REPORT_GENOMES {
    label 'process_cohort'
    container params.python_container
    cache 'deep'
    publishDir params.outdir, mode: 'copy', overwrite: true,
        saveAs: { name -> name.startsWith('report/') ? name.substring(7) : null }
    input:
    path plan
    path master, name: 'master.tsv'
    path status, name: 'status.tsv'
    path bundles, stageAs: 'bundles/bundle??????.json'
    path reference, name: 'cohort_references.json'
    tuple val(cohort_prefixes), path(cohort_payloads, stageAs: 'cohort_inputs/input????')

    output:
    path 'report/tables/*', emit: tables
    path 'report/annotation_results.json', emit: manifest
    path 'report/annotation_acceptance.json', emit: acceptance
    script:
    def cohortEntries = groovy.json.JsonOutput.toJson((cohort_payloads instanceof Collection ? cohort_payloads : [cohort_payloads]).withIndex().collect { value, index -> [prefix: cohort_prefixes[index], path: value.toString()] })
    def bundleList = groovy.json.JsonOutput.toJson((bundles instanceof Collection ? bundles : [bundles]).collect { it.toString() })
    def referenceArg = reference ? '--reference cohort_references.json' : ''
    """
    cat > cohort_entries.json <<'COHORT_EVIDENCE'
${cohortEntries}
COHORT_EVIDENCE
    cat > bundles.json <<'BUNDLES'
${bundleList}
BUNDLES
    python3 "\$(command -v report_genomes.py)" --plan '${plan}' --master master.tsv --status status.tsv \
        --bundles bundles.json ${referenceArg} --cohort-entries cohort_entries.json --output report
    """
}
