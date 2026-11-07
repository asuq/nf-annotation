import groovy.json.JsonSlurper
import java.nio.file.Files
import java.nio.file.Path

/** Explicit recovery keeps the original publication and scientific plan immutable. */
class AnnotationRecovery {
    static Path source(Object value, Object output) {
        def source = Path.of(value.toString()).toRealPath()
        def destination = Path.of(output.toString()).toFile().canonicalFile.toPath()
        if (source.startsWith(destination) || destination.startsWith(source)) {
            throw new IllegalArgumentException('eggNOG recovery requires a separate, non-overlapping --outdir')
        }
        if (Files.isDirectory(destination) && destination.toFile().listFiles()?.any { it.name != 'pipeline_info' }) {
            throw new IllegalArgumentException('eggNOG recovery requires a fresh --outdir; retain each previous attempt')
        }
        source
    }

    static Map targets(Path source, Path currentPlan, Path currentPreflight) {
        def parser = new JsonSlurper()
        def plan = parser.parse(currentPlan.toFile())
        def previousPlan = parser.parse(source.resolve('pipeline_info/annotation_plan.json').toFile())
        def identityFields = ['accession', 'tool', 'input_proteins', 'batch_id', 'search_fingerprint',
            'normalization_fingerprint', 'method_id', 'container', 'resource', 'cpus', 'memory_gib']
        def identities = { record -> record.tasks.collect { row -> identityFields.collect { row[it] } } }
        if (previousPlan.accessions != plan.accessions || identities(previousPlan) != identities(plan) ||
            parser.parse(source.resolve('pipeline_info/annotation_preflight.json').toFile()) != parser.parse(currentPreflight.toFile())) {
            throw new IllegalArgumentException('eggNOG recovery requires the identical cohort, tools, inputs, resources and execution plan')
        }
        def targets = [:]
        plan.tasks.findAll { it.tool == 'eggnog' && it.batch_id }.groupBy { it.batch_id }.each { batchId, members ->
            def packet = source.resolve("annotation_batches/${batchId}")
            def record = parser.parse(packet.resolve('batch_result.json').toFile())
            if (record.exit_code != 0) {
                if (members.any { it.action != 'run' }) {
                    throw new IllegalArgumentException('Failed recovery batch does not match a native run task')
                }
                targets[members[0].task_directory] = [packet, members.collect {
                    source.resolve("samples/${it.accession}/annotation/eggnog")
                }]
            }
        }
        if (!targets) { throw new IllegalArgumentException('No failed eggNOG batches in the recovery source') }
        if (plan.tasks.any { it.action != 'skip' && it.action != 'reuse' &&
            !(it.action == 'run' && targets.containsKey(it.task_directory)) }) {
            throw new IllegalArgumentException('eggNOG recovery cannot rerun or reinterpret other native results')
        }
        targets
    }
}
