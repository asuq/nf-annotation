#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { REFERENCE_ACCEPTANCE } from './modules/local/reference_acceptance'

include { ANNOTATION_RESOURCES; FUNCTIONAL_ANNOTATION } from './subworkflows/local/functional_annotation'
include { IMPORT_ANNOTATION_SOURCE } from './modules/local/import_annotation_source'
include { AGGREGATE_ANNOTATIONS; ANNOTATION_ACCEPTANCE } from './modules/local/functional_annotation'

workflow {
    SamplePolicy.rejectRemovedOptions(params)
    def fingerprintWorkers = params.source_fingerprint_workers.toString()
    if (!fingerprintWorkers.isInteger() || (fingerprintWorkers as int) < 1 || (fingerprintWorkers as int) > 8) {
        error '--source_fingerprint_workers must be an integer from 1 to 8'
    }
    if (SamplePolicy.stage(params) != 'full') { error '--pipeline_stage applies to main.nf, not reannotate.nf' }
    if (params.ani_publication_mode != 'self_contained') { error 'reannotate.nf produces self-contained outputs; --ani_publication_mode applies to ANI-only main.nf runs' }
    if (!params.annotation_from) { error '--annotation_from requires native v0.4 published results' }
    def source = file(params.annotation_from, checkIfExists: true).toFile().canonicalFile.toPath()
    def destination = file(params.outdir).toFile().canonicalFile.toPath()
    def roots = ReferenceInputs.paths([source])
    def referenceSource = new groovy.json.JsonSlurper().parse(source.resolve('annotation_results.json').toFile()).pipeline_contract in ['nf-annotation-ani-reference-v1', 'nf-annotation-ani-reference-v2']
    if (roots.any { it.startsWith(destination) || destination.startsWith(it) }) {
        error '--annotation_from and --outdir must be separate, non-overlapping directories'
    }
    if (destination.exists() && destination.toFile().listFiles()?.any { it.name != 'pipeline_info' } && !workflow.resume) {
        error 'Reannotation requires a fresh output directory; use -resume for an interrupted run'
    }
    ANNOTATION_RESOURCES()
    IMPORT_ANNOTATION_SOURCE(Channel.value(source), ANNOTATION_RESOURCES.out.receipt, Channel.value(roots.findAll { it != source }))
    bundles = IMPORT_ANNOTATION_SOURCE.out.samples
        .splitCsv(header: true, sep: '\t')
        .flatMap { sample ->
            // Import validates and republishes the source. Read its immutable
            // original bundle path so a fresh import does not invalidate resume.
            def bundle = source.resolve("samples/${sample.accession}/annotation/bundle")
            bundle.resolve('bundle.json').exists() ? [tuple([accession: sample.accession], bundle.toRealPath())] : []
        }
    FUNCTIONAL_ANNOTATION(
        IMPORT_ANNOTATION_SOURCE.out.samples, bundles,
        ANNOTATION_RESOURCES.out.receipt, Channel.value(source),
    )
    AGGREGATE_ANNOTATIONS(
        FUNCTIONAL_ANNOTATION.out.plan, IMPORT_ANNOTATION_SOURCE.out.master,
        IMPORT_ANNOTATION_SOURCE.out.sample_status, FUNCTIONAL_ANNOTATION.out.bundle_files,
        FUNCTIONAL_ANNOTATION.out.results.map { item -> item[1] }.toList(),
        FUNCTIONAL_ANNOTATION.out.native_batches.toList(),
        Channel.value([]),
        Channel.value([]),
    )
    if (referenceSource) {
        REFERENCE_ACCEPTANCE(AGGREGATE_ANNOTATIONS.out.acceptance, IMPORT_ANNOTATION_SOURCE.out.samples,
            Channel.value(source), Channel.value(roots.findAll { it != source }), IMPORT_ANNOTATION_SOURCE.out.snapshot)
    } else {
        ANNOTATION_ACCEPTANCE(AGGREGATE_ANNOTATIONS.out.acceptance, IMPORT_ANNOTATION_SOURCE.out.samples)
    }
}
