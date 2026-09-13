# Change Log

## v0.4.0 - 2026-09-13

- Published the qualified Linux amd64 annotation images and pinned their OCI
  digests, including the shared Python helper, in the default configuration.
- Replaced `--ani_allow_incomplete_16s` with `--ani_16s_policy`: `complete`
  (default) requires `Yes`, `allow_incomplete` accepts `Yes` and `partial`,
  and `ignore` also admits `No` and `NA`. Other ANI quality filters still apply.
- Required an active container engine and matching effective process runtimes
  for functional annotation. Native wrapper identities invalidate earlier
  searches that did not enforce this execution contract.
- Preserved failed native tasks as failures in the Nextflow cache while retaining
  their diagnostics for reporting. Explicit resume reruns failed searches and
  retains successful cached searches.
- Unified canonical gene-ID encoding between protein bundles and eggNOG batches,
  allowing encoded original accession and protein identifiers to survive batching.

- Set eggNOG's DIAMOND ceiling to `sensitive` with iteration enabled, matching
  the pinned mapper's default and reducing search cost for large-cohort work.
  The changed search identity invalidates earlier ultra-sensitive eggNOG
  results. The later execution-wrapper identity change also invalidates
  earlier unverified searches for all five callers. Five-genome native
  equivalence passes; 10,000-sample performance remains unqualified.
- Added deterministic whole-proteome eggNOG batching with a 4 MiB FASTA target,
  one shared DIAMOND search, separate native annotation per proteome and one
  immutable archive per batch. Seed partitions retain the original native fields
  with explicit derivation receipts. A native confidence difference rejected the
  initial pooled-annotation design; no pooled annotation feeds primary results.
  Batch membership participates in search invalidation, and portable reuse
  includes the shared `annotation_batches/` directory. The corrected workflow
  matches all compared native fields and normalized tables from five independent
  sensitive-mode proteome runs.
- Bounded planning command arguments through JSON path lists and reduced
  aggregation memory through SQLite spooling. Source-table validation projects
  required fields while streaming; ANI matrix reading uses two passes while
  retaining the dense clustering algorithm and its quadratic limits.
- Added shared annotation planning and execution for eggNOG-mapper v3,
  COGclassifier, direct Pfam, KOfamScan and PADLOC, with separate pinned runtimes
  and prepared resource identities.
- Added code-aware protein and coordinate bundles, detailed native evidence,
  authoritative tool status, expanded master fields and six source-specific
  gene-count matrices. Missing or invalid evidence remains unavailable.
- Added `reannotate.nf` and integrated functional reuse into cohort updates.
  Published v0.4 results support reuse and renormalization without an old work
  directory; changed search methods or inputs invalidate the affected analyses.
- Made paired CheckM2 mean gene length the default genetic-code criterion:
  table 4 for a valid ratio above 1.5 and table 11 for a valid ratio at or below
  1.5. Invalid pairs remain unresolved. Retained the two selectable
  completeness-difference rules.
- Removed v0.3 result import assumptions, the historical ANI rescue CLI,
  eggNOG-only accession selection and free-form eggNOG/PADLOC arguments.
  The three-policy 16S ANI gate remains independent of annotation.
- Hardened native container mounts, runtime overrides, task temporary storage,
  PADLOC protein joins and offline BUSCO preparation. Full biological
  qualification remains tracked in the [release specification](development/v0.4.md).
- Recorded native Prokka ambiguity masking in source-contig provenance and
  enabled bundle recovery from retained native files during cohort updates.
  Removed overlapping parent-directory publication in both import paths so
  repeated cached resumes preserve bundles and functional results.
- Preserved large validation-detail cells when reading published annotation
  tables, allowing complete native evidence to survive cohort updates and
  reannotation despite Python's default CSV field-size limit.

## v0.3.1 - 2026-09-10

- Added `--update_from` to reuse validated published sample results, analyse
  added accessions, and rebuild ANI and complete cohort reports without the
  previous work directory or cache. Updates write independent result copies
  with source/genome provenance and preserve retained internal IDs and tool
  outcomes.
- Added explicit empty-ANI outputs for revised cohorts with no eligible
  genomes. Invalid or contradictory ANI inputs still fail.
- Staging now uses a distinct input filename, preventing an input FASTA from
  being overwritten when its basename matches the internal-ID output name.
- Updated the shared `nf-helper` profiles: OIST uses Apptainer by default,
  and GWDG no longer selects the retired `medium` partition.

## v0.3.0 - 2026-08-01

- Added an ANI recovery CLI with detailed logging, reusable published assembly
  statistics, safer manifest handling, and parallelised rescue workers.
- Added sequence-derived GC content, `is_new` reporting, optional inclusion of
  genomes with incomplete 16S calls in ANI clustering, and separate
  low-quality 16S cohort outputs.
- Added reusable `nf-helper` site configuration, GWDG SCC and MPCDF Viper CPU
  profiles, and bounded Viper Slurm submission and status polling behaviour.
- Updated the workflow for the strict Nextflow 26 syntax parser and set
  Nextflow 26 as the minimum supported release line.
- Improved BUSCO dataset preparation, lineage parsing, output compression, and
  summary parsing while avoiding concurrent writes to shared databases.
- Improved assembly-statistics failure reporting, temporary-file handling,
  parallel GC backfilling, and protection against pipeline race conditions.
- Corrected per-sample publication layout, retained Prokka GenBank output, and
  hardened final-output channel handling.
- Added a locked Pixi test environment and expanded configuration, integration,
  HPC, reporting, and recovery regression coverage.

This release preserves the software, container, database, and `nf-helper`
revisions already pinned by the latest commit on `main`.

## 2026-03-07

- Scaffolded the repository skeleton from `docs/design_spec.md`.
- Added placeholder Nextflow entrypoints, config files, docs, and assets.
