# Implemented Pipeline Overview

## Purpose and scope

This document describes the pipeline as it is currently implemented in code.
It is aimed at developers and maintainers who need a concise view of the
workflow structure, data flow, published artefacts, and reporting contracts.

It complements `docs/runbook.md`, which remains the operator-facing source for
execution examples, acceptance workflows, and troubleshooting guidance. This
overview is intentionally centred on implementation, not on the historical v1
design freeze.

## Workflow entrypoints

The repository currently exposes two Nextflow DSL2 entrypoints.

The declared minimum Nextflow version is `>=26.04.0`, matching the strict
syntax parser used by current stable Nextflow releases.

### `main.nf`

Use the main entrypoint for normal analysis runs. It orchestrates:

- input validation and genome staging
- taxonomy expansion from a user-supplied taxdump
- per-sample QC with Barrnap, paired CheckM2 runs, and BUSCO
- per-sample Codetta execution and summary generation for all samples
- cohort 16S aggregation
- gcode-gated annotation with Prokka, CRISPRCasFinder, PADLOC, and eggNOG
- ANI preparation, all-vs-all FastANI, clustering, and representative selection
- final table and provenance reporting

`main.nf` hard-fails unless these parameters are set:

- `sample_csv`
- `metadata`
- `taxdump`
- `checkm2_db`
- `codetta_db`
- `eggnog_db`
- `busco_lineages`

BUSCO lineage datasets are resolved by `BUSCO_DATASET_PREP`. The current
implementation reuses existing lineage directories below `busco_db`. If
`prepare_busco_datasets = true`, missing configured lineages are downloaded
under `busco_db` while existing lineages are reused. The downloader stages the
BUSCO DB parent as a `path` input so container profiles can bind the shared
destination for writes. Otherwise, missing lineages fail during preflight with
a message naming the expected path. Main workflow BUSCO resolution requires
either:

- `busco_db`, or
- `prepare_busco_datasets = true`

### `prepare_databases.nf`

Use the database-preparation entrypoint before analysis when one or more
runtime databases are missing on the target machine. It is separate from the
main workflow and prepares tool-consumable directories in place.

The workflow currently supports:

- curated taxdump preparation
- CheckM2 database preparation
- Codetta profile database preparation
- BUSCO lineage directory preparation
- eggNOG database preparation

It requires a non-empty `busco_lineages` list and at least one destination
parameter among:

- `taxdump`
- `checkm2_db`
- `codetta_db`
- `busco_db`
- `eggnog_db`

## End-to-end execution flow

```text
prepare_databases.nf
  -> RUNTIME_DATABASE_PREP
     -> PREP_TAXDUMP_DATABASE
     -> DOWNLOAD_CHECKM2_DATABASE
     -> FINALISE_RUNTIME_DATABASE
     -> PREP_CODETTA_DATABASE
     -> DOWNLOAD_BUSCO_DATABASES
     -> FINALISE_RUNTIME_DATABASE
     -> DOWNLOAD_EGGNOG_DATABASE
     -> FINALISE_RUNTIME_DATABASE
     -> MERGE_RUNTIME_DATABASE_REPORTS

main.nf
  -> INPUT_VALIDATION_AND_STAGING
     -> VALIDATE_INPUTS
     -> STAGE_INPUTS
  -> BUSCO_DATASET_PREP
  -> COHORT_TAXONOMY
     -> TAXONOMY_EXPAND
  -> PER_SAMPLE_QC
     -> BARRNAP
     -> per-sample `best_16S.fna` and `16S_status.tsv`
     -> CHECKM2 (ttable 4)
     -> CHECKM2 (ttable 11)
     -> ASSIGN_GCODE_AND_QC
     -> BUSCO (each configured lineage)
  -> COHORT_16S
     -> BUILD_COHORT_16S
     -> all_best_16S.fna
     -> all_partial_16S.fna
     -> low_quality_best_16S.fna
     -> low_quality_partial_16S.fna
  -> PER_SAMPLE_ANNOTATION
     -> CODETTA
     -> SUMMARISE_CODETTA
     -> PROKKA
     -> CCFINDER
     -> SUMMARISE_CCFINDER
     -> PADLOC
     -> EGGNOG
  -> COHORT_ANI
     -> SUMMARISE_BUSCO
     -> CALCULATE_ASSEMBLY_STATS
     -> BUILD_FASTANI_INPUTS
     -> FASTANI
     -> CLUSTER_ANI
  -> FINAL_OUTPUTS
     -> SELECT_ANI_REPRESENTATIVES
     -> BUILD_MASTER_TABLE
     -> WRITE_SAMPLE_STATUS
     -> COLLECT_VERSIONS
```

The orchestration layer stays in Nextflow. Parsing, summarisation, join logic,
and final table assembly live in small CLIs under `bin/`.

## Key inputs and runtime parameters

The most important implementation-level parameters are:

| Parameter | Entrypoint | Role |
| --- | --- | --- |
| `sample_csv` | `main.nf` | Input manifest for all requested genomes. |
| `update_from` | `main.nf` | Optional published source cohort; reuse retained samples, analyse additions, and rebuild all cohort reports in a separate output. |
| `metadata` | `main.nf` | Metadata block preserved into the final master table. |
| `taxdump` | `main.nf`, `prepare_databases.nf` | Resolved taxdump directory for taxonomy expansion or preparation target. |
| `checkm2_db` | `main.nf`, `prepare_databases.nf` | CheckM2 database directory. |
| `codetta_db` | `main.nf`, `prepare_databases.nf` | Codetta profile directory rooted at `Pfam-A_enone.hmm`. |
| `busco_db` | `main.nf`, `prepare_databases.nf` | Offline BUSCO lineage root unless datasets are downloaded on demand. |
| `prepare_busco_datasets` | `main.nf` | Switches BUSCO lineage resolution from reuse to download. |
| `busco_lineages` | both | Non-empty lineage list; defaults to `bacillota_odb12` and `mycoplasmatota_odb12`. |
| `eggnog_db` | `main.nf`, `prepare_databases.nf` | eggNOG data directory. |
| `gcode_rule` | `main.nf` | Translation-table assignment rule; defaults to `strict_delta`. |
| `ani_threshold` | `main.nf` | FastANI clustering threshold; defaults to `0.95`. |
| `eggnog_only_accessions` | `main.nf` | Optional accession allow-list for eggNOG execution. |
| `outdir` | both | Published output root; defaults to `results`. |
| `download_missing_databases` | `prepare_databases.nf` | Enables in-place population of missing runtime databases. |
| `force_runtime_database_rebuild` | `prepare_databases.nf` | Forces re-preparation of runtime database destinations. |
| `runtime_db_helper_container` | `prepare_databases.nf` | Dedicated helper image for runtime database preparation and finalisation. |
| `runtime_db_scratch_root` | `prepare_databases.nf` | Optional scratch root for database preparation. |

Version-labelling parameters are reported in `tool_and_db_versions.tsv` but do
not change workflow behaviour:

- `taxdump_label`
- `checkm2_db_label`
- `codetta_db_label`
- `eggnog_db_label`

Supported profiles currently defined in `nextflow.config` and `conf/*.config`
are:

- `local`
- `docker`
- `slurm`
- `singularity`
- `oist`
- `gwdg`
- `test`
- `debug`

Profile-specific behaviour that affects implementation understanding:

- `debug` sets `eggnog_only_accessions = 'GCA_000027325.1'`
- `test` wires the stub fixtures under `assets/testdata/stub/`
- `docker` adds an amd64 override for the CRISPRCasFinder container on ARM hosts
- `slurm` and `oist` raise medium and high resource ceilings
- `singularity` and `oist` enable Singularity and honour cache and run options
- `gwdg` enables Singularity on SCC, keeps task scratch on `LOCAL_TMPDIR`, prefers `SHM_TMPDIR` for first-attempt temp files, and retries with `LOCAL_TMPDIR`

## Published output layout

The default published root is `results/`.

```text
results/
  cohort/
    ani_clusters/
      ani_representatives.tsv
      cluster.tsv
    16s/
      all_best_16S.fna
      all_best_16S_manifest.tsv
      all_partial_16S.fna
      all_partial_16S_manifest.tsv
      low_quality_best_16S.fna
      low_quality_best_16S_manifest.tsv
      low_quality_partial_16S.fna
      low_quality_partial_16S_manifest.tsv
    fastani/
      ani_exclusions.tsv
      ani_metadata.tsv
      fastani.matrix
      fastani.tsv
      fastani_paths.txt
      fastani_inputs/
    taxonomy/
      taxonomy_expanded.tsv
  pipeline_info/
    dag.html
    report.html
    timeline.html
    trace.tsv
  samples/
    <accession>/
      staged/
      barrnap/
      checkm2_gcode4/
      checkm2_gcode11/
      checkm2/
      busco/
        <lineage>/
      codetta/
      prokka/
      ccfinder/
      padloc/
      eggnog/
  tables/
    accession_map.tsv
    master_table.tsv
    sample_status.tsv
    tool_and_db_versions.tsv
    validated_samples.tsv
    validation_warnings.tsv
```

Notes on that layout:

- `samples/<accession>/...` uses the original accession as the published folder
  name, not the internal sanitized ID
- `codetta/` is published for every sample and contains the raw Codetta outputs,
  `codetta.log`, and `codetta_summary.tsv`
- `prokka/`, `ccfinder/`, `padloc/`, and `eggnog/` are only published for
  samples whose assigned gcode is `4` or `11`
- `busco/<lineage>/` is published for every configured lineage per sample
- `checkm2_gcode4/`, `checkm2_gcode11/`, `busco/<lineage>/`, `prokka/`, and
  `ccfinder/` are flat curated publish folders: they keep the published logs
  and summary artefacts at the tool root and do not republish duplicate nested
  raw subdirectories when those would only mirror the same files
- `cohort/16s/` contains the combined eligible-sample and low-quality 16S
  FASTA and manifest outputs
- `pipeline_info/` is enabled through base configuration for both entrypoints

The database-preparation workflow publishes these top-level artefacts under its
own `--outdir`:

- `runtime_database_report.tsv`
- `nextflow_args.txt`

## Final reporting contracts

Final reporting is assembled in `FINAL_OUTPUTS` from stable intermediate
summaries rather than by re-reading raw tool directories.

### `master_table.tsv`

The master table preserves the original metadata block exactly as supplied and
in its original order. Derived columns are appended afterwards in the runtime
order resolved from `params.busco_lineages`. The appended block starts with
`is_new`, sourced from `validated_samples.tsv`. The checked-in
`assets/master_table_append_columns.txt` snapshot records the default lineage
pair only.

```text
is_new
superkingdom
phylum
class
order
family
genus
species
Completeness_gcode4
Completeness_gcode11
Contamination_gcode4
Contamination_gcode11
Coding_Density_gcode4
Coding_Density_gcode11
Average_Gene_Length_gcode4
Average_Gene_Length_gcode11
Total_Coding_Sequences_gcode4
Total_Coding_Sequences_gcode11
GC_Content
Gcode
Codetta_Genetic_Code
Codetta_NCBI_Table_Candidates
Low_quality
16S
BUSCO_bacillota_odb12
BUSCO_mycoplasmatota_odb12
CRISPRS
SPACERS_SUM
CRISPR_FRAC
Cluster_ID
Is_Representative
ANI_to_Representative
Score
```

### `sample_status.tsv`

The status table is created twice:

1. `VALIDATE_INPUTS` seeds one row per accession with identity fields and
   initial validation state.
2. `WRITE_SAMPLE_STATUS` overlays downstream module outcomes and ANI inclusion
   state to produce the final audit table.

The final column contract is resolved at runtime from `params.busco_lineages`.
The checked-in `assets/sample_status_columns.txt` snapshot records the default
lineage pair only.

```text
accession
internal_id
is_new
validation_status
taxonomy_status
barrnap_status
checkm2_gcode4_status
checkm2_gcode11_status
gcode_status
gcode
low_quality
busco_bacillota_odb12_status
busco_mycoplasmatota_odb12_status
codetta_status
prokka_status
ccfinder_status
padloc_status
eggnog_status
ani_included
ani_exclusion_reason
warnings
notes
```

### `tool_and_db_versions.tsv`

The versions table merges per-process `versions.yml` files with workflow-level
runtime context. It reports tool versions, container references, database paths
or labels, pipeline metadata, and the active container engine in one final TSV.

## Implementation notes and current limits

- Original accession and internal ID are handled separately. The original
  accession is kept in tables and published paths; the sanitized internal ID is
  only used for execution-safe filenames and tool prefixes.
- Barrnap runs for every sample and also emits the per-sample `best_16S.fna`
  and `16S_status.tsv` artefacts used by downstream ANI and reporting. The
  cohort 16S branch builds `all_best_16S.fna` from intact hits that are
  non-atypical or atypical only because of `unverified source organism`;
  mixed atypical reasons remain excluded. It builds `all_partial_16S.fna`
  from every non-low-quality partial-only sample, including atypical partial
  samples. Literal `Low_quality = true` samples are moved to
  `low_quality_best_16S.fna` or `low_quality_partial_16S.fna`, including
  atypical samples.
- CheckM2 always runs twice per sample, once with translation table `4` and
  once with `11`. `summarise_checkm2.py` applies `params.gcode_rule` and emits
  the merged per-sample QC summary.
- BUSCO is independent of gcode assignment. It runs offline in genome mode for
  every sample and every configured lineage.
- Codetta is independent of gcode assignment. It runs for every sample and
  adds `Codetta_Genetic_Code`, `Codetta_NCBI_Table_Candidates`, and
  `codetta_status` to the final reporting layer.
- `calculate_assembly_stats.sh` uses `seqtk comp` to emit `n50`, `scaffolds`,
  `genome_size`, and the derived `gc_content` cohort table used by ANI rescue
  and to populate the master-table `GC_Content` column.
- The default Codetta image is `quay.io/asuq1617/codetta:2.0`. The workflow
  reports Codetta as `v2.0` and also records the pinned upstream source commit
  used to build that image in `tool_and_db_versions.tsv`.
- `PER_SAMPLE_ANNOTATION` is only partly gated by the assigned gcode. Codetta
  runs for every sample, while only samples with gcode `4` or `11` run Prokka,
  CRISPRCasFinder, PADLOC, and eggNOG.
- PADLOC and eggNOG outputs are intentionally retained in sample folders but
  intentionally excluded from `master_table.tsv`.
- ANI clustering is driven by `BUILD_FASTANI_INPUTS`, which writes both the
  FastANI path list and the accession-keyed eligibility report
  `ani_exclusions.tsv`. By default only `16S = Yes` samples pass the 16S ANI
  gate; the naked flag `--ani_allow_incomplete_16s` also allows `16S = No`
  and `16S = partial`, while `16S = NA` remains excluded.
- The main workflow does not prepare taxdump, CheckM2, or eggNOG at runtime.
  Those resources must already exist or be prepared ahead of time with
  `prepare_databases.nf`.
- BUSCO lineage datasets are the only runtime database that `main.nf` can
  resolve directly through `prepare_busco_datasets = true`.
- The current debug profile is deliberately not a full biological run because
  it short-circuits eggNOG to a single accession by default.

## Related docs

- `docs/runbook.md`: operator usage, profiles, acceptance runs, and debugging
- `docs/design_spec.md`: frozen v1 behaviour contract and original design scope
- `docs/run_pipeline_test.md`: wrapper-specific usage for `bin/run_pipeline_test.sh`
- `docs/change_log.md`: high-level repository change log
