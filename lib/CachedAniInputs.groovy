import java.nio.file.Path

/** Preserve equivalent collector bytes so execution-only changes retain ANI cache keys. */
class CachedAniInputs {
    static List<Path> find(Object launch, Object work, boolean enabled, boolean hasMatrix) {
        if (!enabled || !hasMatrix) { return [] }
        def root = new File(work.toString()).canonicalFile.toPath()
        def candidates = new LinkedHashSet<Path>()
        // Nextflow rotates the launch log before executing the new workflow.
        // Logs locate candidates only; fresh table contents authorise reuse below.
        (1..9).each { number ->
            def log = new File(launch.toString(), ".nextflow.log.${number}")
            if (log.isFile()) {
                log.eachLine { line ->
                    if (line.contains('COHORT_ANI:BUILD_FASTANI_INPUTS')) {
                        def completed = line =~ /Task completed > .*workDir: (.+)\]/
                        if (completed.find()) {
                            def path = new File(completed.group(1)).canonicalFile.toPath()
                            if (path.startsWith(root) && path.resolve('.exitcode').toFile().isFile()
                                && path.resolve('.exitcode').toFile().text.trim() == '0') {
                                candidates.add(path)
                            }
                        }
                        def cached = line =~ /\[([0-9a-f]{2})\/([0-9a-f]+)\] Cached process >/
                        if (cached.find()) {
                            def parent = root.resolve(cached.group(1)).toFile()
                            def prefix = cached.group(2)
                            parent.listFiles()?.findAll { it.isDirectory() && it.name.startsWith(prefix) }?.each { entry ->
                                if (new File(entry, '.exitcode').isFile() && new File(entry, '.exitcode').text.trim() == '0') {
                                    def path = entry.canonicalFile.toPath()
                                    if (path.startsWith(root)) { candidates.add(path) }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (hasMatrix && candidates.isEmpty()) {
            throw new IllegalArgumentException('Cannot retain completed ANI cache: keep the original launch logs and work directory when resuming this publication')
        }
        candidates.toList()
    }

    static boolean sameRows(Path current, Path previous) {
        if (!previous.toFile().isFile()) { return false }
        def left = current.toFile().readLines('UTF-8')
        def right = previous.toFile().readLines('UTF-8')
        if (left.isEmpty() || right.isEmpty() || left[0] != right[0]) { return false }
        if (left == right) { return true }
        // Do not treat physical lines as rows for quoted/multiline TSV fields.
        // These require unchanged bytes/order, rather than a guessed equivalence.
        if ((left + right).any { it.contains('"') }) { return false }
        def width = left[0].split('\t', -1).length
        if ((left.drop(1) + right.drop(1)).any { it && it.split('\t', -1).length != width }) { return false }
        left.drop(1).sort(false) == right.drop(1).sort(false)
    }

    static Path table(Path current, List<Path> candidates, String relative) {
        for (def root : candidates) {
            def previous = root.resolve(relative)
            if (sameRows(current, previous)) { return previous }
        }
        if (!candidates.isEmpty()) {
            throw new IllegalArgumentException('Cannot retain completed ANI cache: fresh summary rows differ from retained task inputs; use a new output directory for changed analysis')
        }
        current
    }

    static List<Path> busco(List<Path> current, List<Path> candidates) {
        for (def root : candidates) {
            def directory = root.resolve('busco_tables').toFile()
            def previous = directory.listFiles()?.findAll { it.isFile() }?.sort { it.name }?.collect { it.toPath() }
            if (previous == null || previous.size() != current.size()) { continue }
            def remaining = new ArrayList(current)
            for (def old : previous) {
                def matched = remaining.find { sameRows(it, old) }
                if (matched == null) { break }
                remaining.remove(matched)
            }
            if (remaining.isEmpty()) { return previous }
        }
        if (!candidates.isEmpty()) {
            throw new IllegalArgumentException('Cannot retain completed ANI cache: fresh summary rows differ from retained task inputs; use a new output directory for changed analysis')
        }
        current
    }
}
