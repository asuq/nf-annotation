import groovy.json.JsonSlurper

/** Discover dependency mounts; Python performs authoritative content validation. */
class ReferenceInputs {
    static List paths(List roots) {
        def seen = new LinkedHashSet()
        def visit
        visit = { raw, List active ->
            def root = new File(raw.toString()).canonicalFile.toPath()
            if (root in active || active.size() >= 64) { throw new IllegalArgumentException('Circular or excessively deep publication reference') }
            if (root in seen) { return }
            seen.add(root)
            def reference = root.resolve('ani_references.json').toFile()
            if (reference.exists()) {
                if (java.nio.file.Files.isSymbolicLink(reference.toPath())) {
                    throw new IllegalArgumentException('Reference manifests must be ordinary files')
                }
                def record = new JsonSlurper().parse(reference)
                if (!(record.sources instanceof List)) { throw new IllegalArgumentException('Invalid reference source list') }
                record.sources.each { source ->
                    if (!(source.root instanceof String) || !new File(source.root).isAbsolute()) {
                        throw new IllegalArgumentException('Reference source paths must be absolute')
                    }
                    visit.call(source.root, active + [root])
                }
            }
        }
        roots.each { visit.call(it, []) }
        seen.toList()
    }
}
