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

## Optional reference publications

`--ani_publication_mode self_contained` is the default described above. Choose
`--ani_publication_mode reference` explicitly for ANI-only views that retain
links to validated QC payloads instead of independent payload copies. The mode
is rejected for other main-workflow stages. Reannotation always produces
self-contained outputs.

Reference publications contain their own cohort reports and an
`ani_references.json` dependency manifest. The manifest binds each permitted
link to a declared sample in an immutable source snapshot. Missing sources,
changed contents, retargeted/undeclared links and circular dependencies fail.
Reference inputs are checked at the start, and dependencies are checked again
before reference output or reannotation acceptance.

Both `--update_from` and `reannotate.nf` accept reference publications. A later
ANI run can create another reference view or choose the default self-contained
mode to materialise independent data. Links are flattened to their physical
providers, so normal chains of views do not require the intermediate views.
If a source view owns a rebuilt bundle, that view remains a declared dependency
for that bundle. Up to 64 dependency levels are accepted; circular graphs are
always rejected.

Keep every declared provider at its recorded path and unchanged while a reference
view is needed. Treat references and their sources as read-only: editing through
a link edits its source. This option does not change source permissions. Removing
an old work directory is safe; removing a referenced QC publication is not.
Older pipeline versions reject the explicit reference publication contract;
self-contained outputs retain the existing native v0.4 contract.

Reference mode saves persistent payload storage, but still writes reports,
small channel files and any necessary rebuilt bundle. Runtime ANI inputs and
integrity reads also remain. It is not a promise of zero I/O or faster validation.

## Resource settings

The optional defaults are:

```groovy
params {
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
