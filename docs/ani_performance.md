# ANI reuse and performance

## Reporting and reuse

- Runs with functional searches disabled publish genome-level reports. They retain
  canonical sample bundles but do not rebuild cohort protein manifests, gene
  coordinate tables, protein function summaries or empty functional matrices.
- Full functional reporting streams tables in accession order. It retains one
  sample's evidence and small sequential sparse-count files, without SQLite.
- Clustering runs reuse the source publication's `cohort/assembly_stats/assembly_stats.tsv`.
  They do not reread genomes to recalculate these statistics. Source membership
  and the fresh source fingerprint bind the reused table to the requested data.
- One source task fingerprints content and performs semantic validation with the
  same digests. It checks file identity before reusing a digest within that task.
  On resume, semantic receipts are reused only after the fresh content fingerprints
  and validator code match; file metadata alone never authorises reuse.
- Fingerprint readers replenish work as files finish. Final digest ordering stays
  deterministic even when a large file is slow.
- Clustering and representative selection share one matrix load. Parsing is one
  pass, with strict values, dimensions, identity and mutation checks. Complete
  linkage, threshold handling, missing ANI values and scoring rules are unchanged.
- Individual functional reuse stages previous evidence directly into its consumer;
  the planner no longer makes an intermediate payload copy. Final outputs remain
  independent of previous work directories.

`report_scope: genome` in the publication
manifest records the narrower report contract. Canonical bundles remain available
for later reannotation, which can build functional reports when requested.

## Resume and upgrades

Use `-resume` with the same v0.5 source, launch directory, cache, work directory,
inputs and output directory. Genome reports are cacheable from their table and
bundle-manifest inputs. Fingerprinting and final reference checks remain fresh.
Worker counts are execution settings and do not change scientific identities.

v0.5 rejects older publications and caches. Start fresh QC when upgrading.
Both full and clustering reuse use bounded imports of existing canonical bundles;
there is no per-sample bundle rebuilding path. Candidate and source genomes that
refer to the same inode are normalised once during cohort validation.
AAI native evidence uses compact binary arrays; threshold and scoring changes
reuse it. See [AAI clustering](aai_clustering.md).

## Publication and storage

`--publication_mode self_contained` remains the default. It makes independent
clones where supported and ordinary copies otherwise. FastANI staging also uses
independent clones/copies, never writable hard links to QC files.

`--publication_mode reference` requires `--qc_publication_root`. Providers
must be original QC publications inside this declared permanent root. References
record relative locations in `cohort_references.json` and flatten successor
dependencies to those providers.
Moving the entire collection preserves its layout; moving just the view does not.
Materialisation and reannotation produce independent payloads.

Source publications must remain unchanged during a run. Final reference
acceptance freshly hashes providers and input snapshots, with shared roots read
once within that invocation. This still costs a full source read. It detects
changes during computation and is not replaced by trusting size/mtime metadata.
Independent controllers each perform their own reads. No production snapshot or
storage permissions are changed by the pipeline.

## Resource settings

Existing defaults remain:

```groovy
params {
    source_fingerprint_workers = 4
    cohort_validation_workers = 4
    import_batches = 8
    import_max_forks = 2
    cohort_time = 8.h
    cohort_memory = 8.GB
}
```

Worker counts accept 1–8 and are capped by task CPUs. Import settings accept
1–64. More readers can saturate shared storage; raising concurrency is not a
substitute for removing duplicate reads. Site resource limits still apply.

Publication requires POSIX storage visible to task nodes. Do not override import
link staging with `stageInMode = 'copy'`. Native FastANI remains one all-versus-all
allocation, and its scientific pairwise workload is unchanged.

## Measurements

The [v0.4.8 qualification](development/v0.4.8_qualification.md) records local
before/after measurements, regression results and their limits. Reporting speedups
must not be presented as whole-pipeline or shared-filesystem speedups. No production
rerun was used to qualify this release.
