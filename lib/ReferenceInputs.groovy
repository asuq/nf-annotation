import groovy.json.JsonSlurper

/** Discover dependency mounts; Python performs authoritative content validation. */
class ReferenceInputs {
    static void requireQcProviders(List roots, Object declaredRoot) {
        def boundary = new File(declaredRoot.toString()).canonicalFile.toPath()
        roots.each { raw ->
            def root = new File(raw.toString()).canonicalFile.toPath()
            def manifest = new JsonSlurper().parse(root.resolve('annotation_results.json').toFile())
            def providers = [root]
            if (manifest.pipeline_contract == 'nf-annotation-ani-reference-v1') {
                throw new IllegalArgumentException('New reference outputs require original QC or v2 inputs; legacy absolute views are not migrated')
            }
            if (manifest.pipeline_contract == 'nf-annotation-ani-reference-v2') {
                def record = new JsonSlurper().parse(root.resolve('ani_references.json').toFile())
                if (!(record.qc_publication_root instanceof String) || new File(record.qc_publication_root).isAbsolute()
                    || !(record.sources instanceof List)) { throw new IllegalArgumentException('Invalid relative QC declaration') }
                def anchor = root.resolve(record.qc_publication_root).normalize()
                providers = record.sources.collect { source ->
                    if (!(source.path instanceof String) || new File(source.path).isAbsolute() || source.path.split('/').contains('..')) {
                        throw new IllegalArgumentException('Invalid relative QC provider path')
                    }
                    anchor.resolve(source.path).toFile().canonicalFile.toPath()
                }
            }
            providers.each { provider ->
                if (!provider.startsWith(boundary) || provider.resolve('ani_references.json').toFile().exists() || java.nio.file.Files.isSymbolicLink(provider.resolve('ani_references.json'))
                    || provider.resolve('cohort/fastani').toFile().exists()
                    || provider.resolve('tables/cohort_update_run.json').toFile().exists()) {
                    throw new IllegalArgumentException('Reference providers must be original QC publications inside --qc_publication_root')
                }
            }
        }
    }

    static List paths(List roots) {
        def seen = new LinkedHashSet()
        def visit
        visit = { raw, List active ->
            def root = new File(raw.toString()).canonicalFile.toPath()
            if (root in active || active.size() >= 64) { throw new IllegalArgumentException('Circular or excessively deep publication reference') }
            if (root in seen) { return }
            seen.add(root)
            def reference = root.resolve('ani_references.json').toFile()
            if (reference.exists() || java.nio.file.Files.isSymbolicLink(reference.toPath())) {
                if (java.nio.file.Files.isSymbolicLink(reference.toPath())) {
                    throw new IllegalArgumentException('Reference manifests must be ordinary files')
                }
                def record = new JsonSlurper().parse(reference)
                if (!(record.sources instanceof List)) { throw new IllegalArgumentException('Invalid reference source list') }
                def boundary = root
                if (record.schema_version == 2) {
                    if (!(record.qc_publication_root instanceof String) || new File(record.qc_publication_root).isAbsolute()) {
                        throw new IllegalArgumentException('QC publication root must be relative')
                    }
                    boundary = root.resolve(record.qc_publication_root).normalize().toFile().canonicalFile.toPath()
                    if (root.startsWith(boundary) || boundary.startsWith(root)) { throw new IllegalArgumentException('QC root overlaps reference output') }
                    seen.add(boundary)
                } else if (record.schema_version != 1) {
                    throw new IllegalArgumentException('Unknown reference schema')
                }
                record.sources.each { source ->
                    if (record.schema_version == 2) {
                        if (!(source.path instanceof String) || new File(source.path).isAbsolute()
                            || source.path.split('/').contains('..')) { throw new IllegalArgumentException('Invalid relative QC provider path') }
                        def provider = boundary.resolve(source.path).toFile().canonicalFile.toPath()
                        if (!provider.startsWith(boundary)) { throw new IllegalArgumentException('QC provider escapes declared root') }
                        visit.call(provider, active + [root])
                    } else {
                        if (!(source.root instanceof String) || !new File(source.root).isAbsolute()) {
                            throw new IllegalArgumentException('Legacy reference source paths must be absolute')
                        }
                        visit.call(source.root, active + [root])
                    }
                }
            }
        }
        roots.each { visit.call(it, []) }
        seen.toList()
    }
}
