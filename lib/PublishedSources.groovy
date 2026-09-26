/** One source remains supported; lists and comma-separated paths add sources. */
class PublishedSources {
    static List items(raw) {
        if (raw == null || raw == '') { return [] }
        return raw instanceof Collection ? raw as List : [raw]
    }

    static List paths(raw) {
        def values = raw instanceof Collection ? raw : (raw == null ? [] : raw.toString().split(',', -1).toList())
        if (values.any { it == null || !it.toString().trim() }) {
            throw new IllegalArgumentException('--update_from contains an empty source path')
        }
        def paths = values.collect { new File(it.toString().trim()).canonicalFile.toPath() }
        if (paths.unique(false).size() != paths.size()) {
            throw new IllegalArgumentException('--update_from contains duplicate source directories')
        }
        return paths
    }

    static String arguments(raw, String option = '--source') {
        return items(raw).collect { option + " '" + it.toString().replace("'", "'\"'\"'") + "'" }.join(' ')
    }
}
