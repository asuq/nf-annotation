process PREPARE_AAI_BATCH {
    tag "batch ${batch_id}"
    label 'process_medium'
    container params.fastaai_container
    stageInMode 'symlink'
    cache false
    maxForks 2
    errorStrategy 'finish'
    maxRetries 0
    publishDir "${params.outdir}/pipeline_info/fastaai", mode: 'copy', overwrite: true,
        saveAs: { name -> name == 'preprocessing_metrics.json' ? "batch_${batch_id}.json" : null }
    input:
    tuple val(batch_id), val(records), path(bundles, stageAs: 'bundles/bundle????')
    path metadata, name: 'metadata.tsv'
    path checkm2, name: 'checkm2.tsv'
    path busco, stageAs: 'busco/lineage????.tsv'
    path stats, name: 'assembly_stats.tsv'
    path sources, stageAs: 'sources/source????'
    path dependencies, stageAs: 'dependencies/root????'
    path destination, name: 'destination'
    val primary_busco_column
    output:
    path 'batch.json', emit: batch
    path 'preprocessing_metrics.json', emit: metrics
    path 'versions.yml', emit: versions
    script:
    def values = { input -> groovy.json.JsonOutput.toJson((input instanceof Collection ? input : [input]).collect { it.toString() }) }
    """
    cat > records.json <<'RECORDS'
${groovy.json.JsonOutput.toJson(records)}
RECORDS
    cat > bundles.json <<'BUNDLES'
${values.call(bundles)}
BUNDLES
    cat > busco.json <<'BUSCO'
${values.call(busco)}
BUSCO
    cat > sources.json <<'SOURCES'
${values.call(sources)}
SOURCES
    prepare_aai.py --records records.json --metadata metadata.tsv --checkm2 checkm2.tsv \
        --busco busco.json --stats assembly_stats.tsv --bundles bundles.json --sources sources.json \
        --destination destination --domain '${params.domain}' --primary-busco-column '${primary_busco_column}' --threads ${task.cpus} --output batch.json
    python3 - <<'PY'
import json
from pathlib import Path
runtime=json.loads(Path('/opt/nf-annotation/provenance/fastaai.json').read_text())
Path('versions.yml').write_text('"${task.process}":\\n  fastaai: "'+runtime['native_version']+'"\\n  fastaai_source: "'+runtime['upstream_commit']+'"\\n  fastaai_markers: "'+runtime['hmm_sha256']+'"\\n')
PY
    """
    stub:
    def paths = { value -> (value instanceof Collection ? value : [value]).collect { it.toString() } }
    def request = [records:records, bundles:paths.call(bundles), busco:paths.call(busco), sources:paths.call(sources),
        metadata:metadata.toString(), checkm2:checkm2.toString(), stats:stats.toString(), destination:destination.toString(),
        primary_busco_column:primary_busco_column, domain:params.domain, threads:task.cpus, output:'batch.json']
    """
    cat > request.json <<'REQUEST'
${groovy.json.JsonOutput.toJson(request)}
REQUEST
    python3 '${projectDir}/tests/fastaai_workflow_fixture.py' prepare request.json
    printf '\"${task.process}\":\\n  fastaai: \"stub\"\\n' > versions.yml
    """

}

process FASTAAI_COMPARE {
    tag "${domain}"
    label 'process_high'
    container params.fastaai_container
    scratch true
    cache 'lenient'
    errorStrategy 'finish'
    maxRetries 0
    publishDir "${params.outdir}/cohort/fastaai", mode: 'copy', overwrite: true,
        saveAs: { name -> name == 'comparison' ? domain : null }
    input:
    tuple val(domain), val(records), path(profiles, stageAs: 'profiles/profile??????')
    output:
    tuple val(domain), path('comparison'), emit: comparison
    path 'versions.yml', emit: versions
    script:
    def paths = groovy.json.JsonOutput.toJson((profiles instanceof Collection ? profiles : [profiles]).collect { it.toString() })
    """
    cat > records.json <<'RECORDS'
${groovy.json.JsonOutput.toJson(records)}
RECORDS
    cat > profiles.json <<'PROFILES'
${paths}
PROFILES
    compare_aai.py --records records.json --profiles profiles.json --output comparison --threads ${task.cpus}
    python3 - <<'PY'
import json
from pathlib import Path
rows=json.loads(Path('records.json').read_text())
for row in rows:
    if row['profile']: row['profile']=Path(row['profile']).name
Path('comparison/records.json').write_text(json.dumps(dict(domain='${domain}',samples=rows),sort_keys=True)+'\\n')
Path('versions.yml').write_text('"${task.process}":\\n  script: "bin/compare_aai.py"\\n')
PY
    """
    stub:
    """
    cat > records.json <<'RECORDS'
${groovy.json.JsonOutput.toJson(records)}
RECORDS
    python3 '${projectDir}/tests/fastaai_workflow_fixture.py' compare records.json --domain '${domain}'
    printf '\"${task.process}\":\\n  fastaai: \"stub\"\\n' > versions.yml
    """

}

process CLUSTER_AAI {
    label 'process_high'
    container params.python_container
    cache 'deep'
    publishDir "${params.outdir}/cohort", mode: 'copy', overwrite: true,
        saveAs: { name -> name == 'aai' ? 'aai_clusters' : null }
    input:
    path domains, stageAs: 'domains/domain??'
    path samples
    output:
    path 'aai/aai_summary.tsv', emit: summary
    path 'aai/aai_representatives.tsv', emit: representatives
    path 'aai/cluster.tsv', emit: clusters
    path 'aai/method.json', emit: method
    path 'aai/scores.json', emit: scores
    path 'aai/metrics.json', emit: metrics
    path 'aai', emit: evidence
    path 'versions.yml', emit: versions
    script:
    def paths = groovy.json.JsonOutput.toJson((domains instanceof Collection ? domains : [domains]).collect { it.toString() })
    """
    cat > domains.json <<'DOMAINS'
${paths}
DOMAINS
    cluster_aai.py --domains domains.json --samples '${samples}' --threshold ${params.aai_threshold} \
        --score-profile '${params.score_profile}' --output aai
    cat > versions.yml <<'VERSIONS'
"${task.process}":
  script: "bin/cluster_aai.py"
VERSIONS
    """
}
