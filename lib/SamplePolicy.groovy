import groovy.json.JsonOutput

/** Resolve declared policy configuration; Python validates individual records. */
class SamplePolicy {
    static String codeHash(Object projectDir) {
        def source = new File(projectDir.toString(), 'bin/sample_policy.py')
        java.security.MessageDigest.getInstance('SHA-256').digest(source.bytes).encodeHex().toString()
    }

    static String reuseCodeHash(Object projectDir) {
        def digest = java.security.MessageDigest.getInstance('SHA-256')
        new File(projectDir.toString(), 'bin').listFiles().findAll { it.name.endsWith('.py') }.sort { it.name }.each { source ->
            digest.update((source.name + '\u0000').getBytes('UTF-8'))
            digest.update(source.bytes)
            digest.update((byte) 0)
        }
        digest.digest().encodeHex().toString()
    }

    static String stage(params) {
        def value = params.pipeline_stage
        if (!(value in ['full', 'qc', 'ani'])) { throw new IllegalArgumentException('pipeline_stage must be full, qc or ani') }
        value.toString()
    }

    static String annotationTools(params) {
        stage(params) == 'full' ? (params.annotation_tools ?: '').toString() : ''
    }

    static void rejectRemovedOptions(params) {
        ['domain_column', 'prokka_kingdom', 'barrnap_kingdom'].each { name ->
            if (params.containsKey(name)) {
                throw new IllegalArgumentException("--${name} is unsupported; use --domain bacteria, archaea or metadata")
            }
        }
    }

    static String requireDomain(params) {
        rejectRemovedOptions(params)
        if (!(params.domain in ['bacteria', 'archaea', 'metadata'])) {
            throw new IllegalArgumentException('--domain is required and must be bacteria, archaea or metadata')
        }
        params.domain.toString()
    }

    static Map configuration(params) {
        [domain: requireDomain(params), busco_domain_lineages: params.busco_domain_lineages,
         busco_primary_column: params.busco_primary_column]
    }

    static List lineages(params, List legacy) {
        def profiles = params.busco_domain_lineages
        if (!(profiles instanceof Map)) {
            throw new IllegalArgumentException('busco_domain_lineages must be a map')
        }
        if (!profiles) { return legacy }
        if (profiles.values().any { !(it instanceof List) || it.isEmpty() }) {
            throw new IllegalArgumentException('BUSCO domain profiles require nonempty lineage lists')
        }
        profiles.values().collectMany { it }.unique(false)
    }

    static String encoded(params) {
        JsonOutput.toJson(configuration(params)).bytes.encodeBase64().toString()
    }

    static List sampleLineages(Map meta, List legacy) {
        meta.busco_lineages ? meta.busco_lineages.toString().split(',').toList() : legacy
    }
}
