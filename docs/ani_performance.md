# ANI-only reuse and scheduler overhead

From v0.4.5, `--pipeline_stage ani` uses bounded cohort bookkeeping instead of
per-genome import, bundle-building, BUSCO-summary and functional-planning jobs.
The existing `--domain` and `--update_from` interfaces are unchanged.

## Execution and integrity

- Import uses at most eight allocations by default, with two running at once.
  Source validation uses at most eight groups. These are real batches of work,
  not arrays containing a task for every genome.
- A source-content fingerprint is computed on every run. Successful source
  validation can be reused only with matching content and validator code.
  Nested file edits invalidate that validation, including edits that preserve
  file size and modification time. Sources must remain unchanged during a run.
- The lightweight resume-identity check is separate from expensive cohort
  validation. The first resume no longer changes the latter's inputs merely
  because `cohort_update_run.json` has appeared. Different cohort inputs or
  settings still require a fresh output directory.
- A successful bundle is retained when its producer/contract hashes, accession,
  genetic code and native-input hashes match. Otherwise the existing bundle
  adapter runs inside the batch; gene prediction is not rerun.
- BUSCO TSVs already checked against their native JSON are reused. One disabled
  annotation plan retains the complete sample/tool grid and its explicit skipped
  statuses. Aggregation still validates bundle contents.
- The four 16S views use one pass through their inputs and the existing inclusion
  functions. The original glob order, FASTA identifiers and manifest rows remain
  unchanged. Native fixed-taxdump lineage expansion remains metadata bookkeeping;
  no genomic classifier is introduced.

ANI distances, the default 95% threshold, quality screening, representative
scoring, genetic-code policy, domain routing and scientific image references
are unchanged. Per-genome failures remain identifiable in batch diagnostics and
`pipeline_info/ani_import/batch_*.json`. Source-validation receipts are retained
in `tables/source_validation.json`; these are audit records, not an option to
trust an unchecked publication.

## Publication and storage

Outputs remain self-contained, with ordinary files and the existing native v0.4
publication contract. Each sample is installed atomically in the destination
after its payload is ready. There is no intermediate full payload copy in the
task work directory followed by another `publishDir` copy.

The importer uses independent filesystem clones on supported APFS/Linux
filesystems and falls back to one normal copy elsewhere. It never creates
writable hard links to the source. Published data remain usable after source
publications and old work directories are removed. A resume verifies existing
sample payloads and refuses inconsistent copies rather than overwriting them.

The optimised path requires POSIX shared storage: the destination must be
writable and visible to the task nodes. Publication-root inputs use link
staging; do not override these bookkeeping processes to `stageInMode = 'copy'`.
Node scratch with copy stage-out for task-owned outputs remains usable.

On filesystems without clone support, independent cumulative publications still
need one physical payload copy each. Content verification still reads source
bytes, and report generation still reads canonical bundles. The optimisation
removes unnecessary copies and repeated processing; it does not make integrity
checks or independent publications free. Copy audits distinguish logical bytes
cloned from logical bytes copied; they do not measure physical storage traffic.

## Independent FastANI staging

The republished v0.4.7 uses the same independent clone/copy helper for
FastANI FASTAs and sample publication. FastANI input creation and cleanup no
longer create/remove hard links to QC files, so they do not change QC inode
link counts or ctime. Only the required FASTAs are staged; this does not copy
whole QC payload trees. This change does not migrate or clean up existing v0.4.7 work directories.

Source validation now shares one fingerprint implementation with final checks
and reuses semantic provider validation within each fresh source-group pass.
Each view's own tables and bindings are still checked. Start and end scans
remain separate because they observe different moments; no cross-run trust
cache or validation-disable option is introduced. Ctime remains part of the
mutation guard, including for same-size edits with restored mtime.

## Relocatable reference publications

`--ani_publication_mode self_contained` remains the default. From v0.4.7, creating
an ANI-only reference view requires `--ani_publication_mode reference` and an
explicit `--qc_publication_root`. The latter declares the permanent QC tree;
all external providers must be original self-contained QC publications within
it, and the ANI output must be outside it. A reference view or an ANI/update
publication cannot itself be an external payload provider. Permanence is an
operator declaration, not something inferred from a directory name.

For example, keep `02-qc/batch-a` and `03-ani/view` inside one collection. Supply
the current `02-qc` path as `--qc_publication_root`, its QC publication(s) through
`--update_from`, and a fresh `03-ani/view` as `--outdir`. The v2
`ani_references.json` records the QC root relative to the view and each provider
relative to that QC root. External filesystem links are relative too. They are
created for their final installed location before atomic publication; no
post-publication path repair is performed.

Moving the whole collection to a different parent preserves validation and
content digests when its internal layout is retained. New invocations use their
current input paths. Historical absolute command/input paths in untouched QC
provenance remain text and are not used as fallback locations. Removing the
original parent and old work directories is safe; removing a referenced QC
provider is not. Moving only one side of the relationship is not supported.

`--update_from` can consume a v2 view and produce another v2 view or a default
self-contained materialisation. `reannotate.nf` accepts relocated views and
always produces self-contained output. Successor views link directly to the
permanent QC publications. If an earlier view owns a necessary canonical bundle,
the successor owns a validated local copy/reconstruction of that bundle; it
never retains an external link to the earlier ANI view. This exception does not
copy whole QC payload trees or rerun gene prediction.

Both the relative manifest and exact link targets are checked. Missing/changed
providers, absolute dependency declarations, boundary escapes, retargeted or
undeclared links, duplicate identities and cycles fail. Providers are immutable
throughout a run. Final acceptance freshly checks permanent providers and the
original input snapshots, including transient input views that need not remain
after completion. Reading or editing through a reference link still accesses its
provider; the workflow does not change source permissions.

Original v0.4.4 QC publications remain supported without editing them. Existing
v1 absolute references remain readable at their original locations, including
for materialisation, but new v2 reference creation accepts original QC or v2
inputs only. There is no automatic v1 migration or legacy-path compatibility
link. Older releases reject the v2 contract. Self-contained outputs retain the
native v0.4 contract.

Reference mode still writes reports, small channel files and any necessary local
bundle. Materialisation explicitly creates independent files. Runtime inputs and
integrity reads remain; this is not a promise of zero I/O.

## Parallel fingerprinting and progress

From v0.4.6, the initial source-content scan uses one bounded pool of file readers.
`--source_fingerprint_workers` defaults to 4 and accepts 1 through 8. The process
requests that many CPUs, subject to `max_cpus` and site resource limits; the reader
count cannot exceed the actual task CPUs. There are no nested pools or additional
per-file/per-genome scheduler jobs. Directory inventories and final hash assembly
remain ordered. Only a bounded window of twice the worker count is submitted;
file payloads are streamed, with at most one open payload per reader. The sorted
path inventory is retained for one source at a time.

Original self-contained QC fingerprints and source order remain unchanged. V2
reference digests include relative manifest/link text and remain stable when the
collection moves; resolved absolute roots in runtime receipts are location context. Every run
reads the complete scientific content and checks file identity/timestamps before
and after hashing. Directory membership changes, leaf replacements, undeclared
links and reference cycles fail. Sources must remain immutable throughout the
run; this is not a transactional filesystem snapshot. Shared reference providers
are read once within one invocation. Separate controllers still perform fresh
reads; there is no persistent cache that trusts old file metadata.

While a task runs, its work directory contains atomically updated
`fingerprint_progress.json` and periodic messages in `.command.err`. Progress
includes completed/requested sources, files and bytes hashed, the current source,
its discovered totals, elapsed time and scan/hash durations. Totals for sources
not yet scanned are not guessed. Completed metrics and the final progress record
are published under `pipeline_info/source_fingerprinting/`. Progress and timings
are separate from scientific fingerprints and resume identities.

More readers can overlap storage waits, but can also saturate shared storage.
The defaults use four readers per controller, not all CPUs in a large allocation.
See the [local qualification](development/v0.4.6_qualification.md); its timings
are not a shared-filesystem speed guarantee.

## Parallel cohort validation and final checks

From v0.4.7, `PREPARE_COHORT_UPDATE` uses `--cohort_validation_workers` (default 4,
range 1-8, capped by allocated task CPUs). Global membership, collisions and ID
allocation are resolved first. One process pool validates samples independently;
results, failure reporting and publication are ordered by the requested cohort.
There are at most twice as many pending jobs as workers, no nested reader pools
and no per-genome scheduler jobs. Candidate/source genome mutations are rejected.

`REFERENCE_ACCEPTANCE` uses `--source_fingerprint_workers` for one shared reader
pool, including final reannotation checks. Providers shared by multiple views
are hashed once within that invocation. These are fresh reads, not a cross-run
trust cache. Worker changes are execution-only and do not alter scientific
identity. Both tasks request the corresponding CPUs within global/site limits.

Live work-directory progress is `cohort_validation_progress.json` or
`reference_progress.json`, with periodic `.command.err` messages. Final records
are copied to `pipeline_info/cohort_validation/` and
`pipeline_info/reference_acceptance/`. Counts, elapsed time and effective workers
are recorded; fingerprint phases also report file/byte counts. Progress is not
part of scientific identities.

This parallelises per-sample cohort checks and final dependency content reads.
Global reductions, directory inventories and some other semantic-validation
loops remain serial or use their existing task batching. Independent controllers
still repeat content reads, and FastANI's pairwise cost is unchanged. See the
[v0.4.7 qualification](development/v0.4.7_qualification.md) for measured limits.

## Retaining completed ANI on resume

Keep the original launch directory, Nextflow cache, work directory and publication
when using `-resume`. Fingerprint and cohort-validation worker counts are execution settings and do
not change cohort identity. A code upgrade can revalidate source evidence. Fresh content checks still
run, while matching validation receipts remain reusable.

v0.4.5 collectors can reorder summary rows when import task paths change. For a
resumed ANI publication with an existing matrix, v0.4.6 locates retained ANI
preparation tasks through up to nine rotated Nextflow launch logs. It retains
prior summary-table bytes only after comparing the fresh headers and complete
row multiplicities; BUSCO tables must match one-to-one. This preserves existing
Nextflow cache keys without reading an unchecked result or modifying the cache.
Quoted/multiline tables require unchanged row order. Missing retained tasks or
nonmatching summaries fail explicitly instead of silently repeating completed
ANI. Native tool settings and inputs must still match for normal Nextflow reuse.
These controls were checked through controlled release-upgrade and worker-change
fixtures. Publication relocation does not make Nextflow work/cache relocatable;
use relocated publications as inputs to a fresh invocation rather than restoring
retired work/cache.

## Resource settings

The optional defaults are:

```groovy
params {
    source_fingerprint_workers = 4
    cohort_validation_workers = 4
    ani_import_batches = 8
    ani_import_max_forks = 2
    cohort_time = 8.h
    cohort_memory = 8.GB
}
```

Both batch settings accept integers from 1 to 64. `ani_import_max_forks` bounds
concurrent import and source-validation groups. Keep the defaults unless
measurements justify changing them. Changing the import partition count requires
a fresh output directory; concurrency can be tuned without changing science.

`cohort_time` and `cohort_memory` set requests for `process_cohort` bookkeeping.
`max_time` and `max_memory` cap those requests; raising a cap alone does not
increase an unrelated process's request. Site `resourceLimits` and partition
limits still apply. In particular, a site's two-hour limit is not removed by
setting `cohort_time = 8.h`.

Batching works with the selected scheduler. Where local bookkeeping is permitted,
an execution-only override can route the cohort label locally and set a separate
budget for each controller:

```groovy
process {
    withLabel: process_cohort {
        executor = 'local'
    }
}
executor {
    $local {
        cpus = 12
        memory = '24 GB'
    }
}
```

This does not reroute FastANI or change scientific settings. On the pinned
Nextflow 26.04 runtime these CPU/memory budgets apply to the local executor;
do not assume the grid-budget support added in later Nextflow versions.
See the [executor configuration documentation](https://docs.seqera.io/nextflow/reference/config/executor)
and [executor-specific scopes](https://docs.seqera.io/nextflow/executor#configuring-executors).
Increasing `queueSize` does not reduce task granularity. Polling and queue-status
intervals remain resource-only tuning choices, and no worker-pool dependency is
introduced.

For sequential cumulative cohorts, keep the same launch/cache and work directories
and use `-resume` with each fresh output directory. Independent concurrent
controllers do not automatically share a validation cache. Do not concurrently
open the same local cache session; see [Nextflow caching and resume](https://docs.seqera.io/nextflow/cache-and-resume).

## Evidence and limits

The [qualification record](development/v0.4.5_qualification.md) distinguishes
synthetic local timings, actual task counts, copy audits, output equivalence and
the Linux helper-image smoke test. These measurements are not a Slurm throughput
guarantee or a qualification of real-genome biological accuracy. Native ANI's
pairwise work and clustering memory requirements remain unchanged.
