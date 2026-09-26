import groovy.json.JsonOutput
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardCopyOption

/** Runtime checks and retained diagnostics at the Nextflow execution boundary. */
class AnnotationExecution {
    static void validateExecutionGroup(Object key, List items) {
        if (items.size() == 1) { return }
        def first = items[0][0]
        def fields = ['task_directory', 'tool', 'batch_id', 'action', 'bundle', 'resource', 'container', 'cpus', 'memory_gib']
        def sharedBatch = first.tool == 'eggnog' && first.batch_id && first.batch_id != 'NA'
        def accessions = items.collect { it[0].accession }
        if (!sharedBatch || accessions.unique(false).size() != items.size() ||
            items.any { item -> fields.any { field -> item[0][field] != first[field] } }) {
            throw new IllegalArgumentException("Duplicate or inconsistent functional execution key in annotation plan: ${key}")
        }
    }

    static void requireEngine(Object engine) {
        if (!(engine in ['docker', 'singularity', 'apptainer'])) {
            throw new IllegalArgumentException('Enabled annotation requires an active Docker, Singularity or Apptainer engine; select a container profile')
        }
    }

    static void requireRuntime(Object engine, Object actual, Object expected) {
        requireEngine(engine)
        if (!expected || !actual || actual.toString() != expected.toString()) {
            throw new IllegalArgumentException("Effective annotation container does not match the planned runtime: ${actual}")
        }
    }

    static Path failureDirectory(Object output, Object started, Map meta) {
        def key = meta.task_directory?.toString()
        if (!key || !(key ==~ /[A-Za-z0-9_-]+/)) {
            throw new IllegalArgumentException('Invalid annotation failure receipt identity')
        }
        return Path.of(output.toString()).toAbsolutePath().normalize()
            .resolve('pipeline_info/annotation_failures').resolve(started.toInstant().toEpochMilli().toString()).resolve(key)
    }

    static void preserveFailure(Object task, Map meta, Object output, Object started) {
        try {
            def destination = failureDirectory(output, started, meta)
            Files.createDirectories(destination)
            def work = task.workDir
            def raw = work?.resolve('raw')
            def diagnosticError = null
            if (raw && Files.isDirectory(raw)) {
                try {
                    def target = destination.resolve('raw')
                    def paths = Files.walk(raw).withCloseable { stream -> stream.toList() }
                    if (paths.any { Files.isSymbolicLink(it) || !(Files.isDirectory(it) || Files.isRegularFile(it)) }) {
                        throw new IOException('Unsafe native diagnostic tree')
                    }
                    if (Files.exists(target)) {
                        def copies = Files.walk(target).withCloseable { stream -> stream.toList() }
                        if (paths.size() != copies.size() || paths.any { source ->
                            def copy = target.resolve(raw.relativize(source))
                            Files.isSymbolicLink(copy) ||
                                (Files.isDirectory(source) ? !Files.isDirectory(copy) :
                                    !Files.isRegularFile(copy) || Files.mismatch(source, copy) != -1)
                        }) { throw new IOException('Existing native diagnostics differ') }
                    } else { raw.copyTo(target) }
                } catch (Exception error) { diagnosticError = error.message }
            }
            for (name in ['annotation_storage.json', '.command.out', '.command.err', '.command.log', '.exitcode']) {
                def source = work?.resolve(name)
                if (source && Files.isRegularFile(source)) {
                    def target = destination.resolve(name)
                    if (Files.isSymbolicLink(source) || Files.isSymbolicLink(target)) {
                        diagnosticError = 'Linked diagnostic control file'
                    } else { Files.copy(source, target, StandardCopyOption.REPLACE_EXISTING) }
                }
            }
            def receipt = [task_directory: meta.task_directory, task_hash: task.hash,
                process: task.process, task_exit_status: task.exitStatus,
                work_directory: work?.toString(),
                diagnostics_complete: diagnosticError == null && Files.isRegularFile(destination.resolve('raw/exit_code.txt')),
                diagnostic_error: diagnosticError]
            Files.writeString(destination.resolve('failure.json'), JsonOutput.prettyPrint(JsonOutput.toJson(receipt)) + '\n')
        } catch (Exception error) {
            System.err.println("Unable to preserve controller-side annotation failure evidence: ${error.message}")
        }
    }

    static Path failedNativeRaw(Object output, Object started, Map meta) {
        def raw = failureDirectory(output, started, meta).resolve('raw')
        def status = raw.resolve('exit_code.txt')
        // An abrupt executor failure may have no native exit record. Its task
        // diagnostics remain retained; aggregation reports the missing outcome.
        if (!Files.isRegularFile(status) || Files.readString(status).trim() == '0') {
            return null
        }
        return raw
    }
}
