import groovy.json.JsonOutput
import java.nio.file.Files
import java.nio.file.Path

/** Runtime checks and retained diagnostics at the Nextflow execution boundary. */
class AnnotationExecution {
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
        def destination = failureDirectory(output, started, meta)
        Files.createDirectories(destination)
        def work = task.workDir
        def raw = work?.resolve('raw')
        if (raw && Files.isDirectory(raw)) {
            raw.copyTo(destination.resolve('raw'))
        }
        for (name in ['.command.out', '.command.err', '.command.log', '.exitcode']) {
            def source = work?.resolve(name)
            if (source && Files.isRegularFile(source)) {
                Files.copy(source, destination.resolve(name))
            }
        }
        def receipt = [task_directory: meta.task_directory, task_hash: task.hash,
            process: task.process, task_exit_status: task.exitStatus,
            work_directory: work?.toString()]
        Files.writeString(destination.resolve('failure.json'), JsonOutput.prettyPrint(JsonOutput.toJson(receipt)) + '\n')
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
