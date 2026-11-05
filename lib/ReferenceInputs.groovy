import groovy.json.JsonSlurper

/** Discover mounts for the single v0.5 publication contract. */
class ReferenceInputs {
    static Map manifest(java.nio.file.Path root) {
        def value = new JsonSlurper().parse(root.resolve('annotation_results.json').toFile())
        if (!(value.pipeline_contract in ['nf-annotation-v0.5', 'nf-annotation-reference-v0.5'])) {
            throw new IllegalArgumentException('Reuse requires v0.5 publications; v0.4 inputs and caches are unsupported')
        }
        value
    }

    static void requireQcProviders(List roots, Object declaredRoot) {
        def boundary = new File(declaredRoot.toString()).canonicalFile.toPath()
        roots.each { raw ->
            def root = new File(raw.toString()).canonicalFile.toPath()
            def record = manifest(root)
            def providers = [root]
            if (record.pipeline_contract == 'nf-annotation-reference-v0.5') {
                def refs = new JsonSlurper().parse(root.resolve('cohort_references.json').toFile())
                def anchor = root.resolve(refs.qc_publication_root).normalize()
                providers = refs.sources.collect { anchor.resolve(it.path).toFile().canonicalFile.toPath() }
            }
            providers.each { provider ->
                if (!provider.startsWith(boundary) || manifest(provider).pipeline_contract != 'nf-annotation-v0.5'
                    || provider.resolve('cohort/fastani').toFile().exists()
                    || provider.resolve('tables/cohort_update_run.json').toFile().exists()) {
                    throw new IllegalArgumentException('Reference providers must be original v0.5 QC publications inside --qc_publication_root')
                }
            }
        }
    }

    static List paths(List roots) {
        def seen = new LinkedHashSet()
        def visit
        visit = { raw, List active ->
            def root = new File(raw.toString()).canonicalFile.toPath()
            if (root in active || active.size() >= 64) { throw new IllegalArgumentException('Circular publication reference') }
            if (root in seen) { return }
            seen.add(root)
            def publication = manifest(root)
            if (publication.pipeline_contract == 'nf-annotation-reference-v0.5') {
                def reference = root.resolve('cohort_references.json')
                if (java.nio.file.Files.isSymbolicLink(reference)) { throw new IllegalArgumentException('Reference manifests must be ordinary files') }
                def record = new JsonSlurper().parse(reference.toFile())
                if (record.schema_version != 1 || !(record.sources instanceof List)
                    || !(record.qc_publication_root instanceof String) || new File(record.qc_publication_root).isAbsolute()) {
                    throw new IllegalArgumentException('Invalid v0.5 relative reference manifest')
                }
                def boundary = root.resolve(record.qc_publication_root).normalize().toFile().canonicalFile.toPath()
                if (root.startsWith(boundary) || boundary.startsWith(root)) { throw new IllegalArgumentException('QC root overlaps reference output') }
                seen.add(boundary)
                record.sources.each { source ->
                    if (!(source.path instanceof String) || new File(source.path).isAbsolute() || source.path.split('/').contains('..')) {
                        throw new IllegalArgumentException('Invalid relative QC provider path')
                    }
                    def provider = boundary.resolve(source.path).toFile().canonicalFile.toPath()
                    if (!provider.startsWith(boundary)) { throw new IllegalArgumentException('QC provider escapes declared root') }
                    visit.call(provider, active + [root])
                }
            }
        }
        roots.each { visit.call(it, []) }
        seen.toList()
    }
}
