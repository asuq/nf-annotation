/** Estimate serial clustering workspace from its actual cohort representation. */
class ClusterResources {
    static final long GIB = 1024L * 1024L * 1024L

    static String bounded(double bytes, Object maximum) {
        def cap = maximum as nextflow.util.MemoryUnit
        long requested = Math.min(Math.ceil(bytes / GIB) * GIB, cap.toBytes()) as long
        return "${requested} B"
    }

    static String ani(Object matrix, Object maximum) {
        def header = matrix.toRealPath().toFile().withReader { reader ->
            def line = reader.readLine()
            while (line != null && !line.trim()) { line = reader.readLine() }
            line?.trim()
        }
        long count = header?.isInteger() ? header.toLong() : 0L
        // Dense ANI/distance arrays, condensed linkage storage and fixed overhead.
        return bounded(2d * GIB + 32d * count * count, maximum)
    }

    static String aai(Object inputs, Object maximum) {
        def paths = inputs instanceof Collection ? inputs : [inputs]
        long largest = paths.collect { it.toRealPath().resolve('jaccard.npy').toFile().length() }.max() ?: 0L
        // Domains run sequentially; retain one condensed domain plus linkage workspace.
        return bounded(2d * GIB + 4d * largest, maximum)
    }
}
