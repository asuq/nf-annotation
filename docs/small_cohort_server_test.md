# Small-Cohort Server Test Guide

## Summary

Use the tracked acceptance cohort as the first real server test. It is small
enough to run cheaply, but it still exercises the important behaviours:

- gcode 4 and 11
- CRISPR positive and negative
- ANI collisions
- missing metadata
- atypical exclusion
- the atypical exception

The cohort definition lives in `assets/testdata/acceptance/cohort_plan.tsv`
and currently contains 9 samples.

Use the `nextflow` conda environment for every Nextflow command. For small
cohorts, do not start with the 20k storage override. First prove correctness
with the normal SLURM path, then test the 20k profile separately if needed.

## 1. Activate The Environment

On the login node:

```bash
eval "$(mamba shell hook --shell bash)"
mamba activate nextflow
python3 --version
nextflow -version
```

Success criteria:

- `python3` is available and recent enough for the repo tooling
- `nextflow` resolves from the activated environment

## 2. Run The Frozen Stub Smoke Test

Before any real-data server run, confirm the checkout still passes the stub
path:

```bash
nextflow run . -profile test -stub-run
```

Success criteria:

- the pipeline exits `0`
- `results/tables/master_table.tsv` exists
- `results/tables/sample_status.tsv` exists
- `results/tables/tool_and_db_versions.tsv` exists

## 3. Prepare The Tracked Small Real Cohort

Generate the tracked acceptance inputs under a dedicated test root:

```bash
python3 bin/run_acceptance_tests.py prepare \
  --work-root /path/to/test_root
```

This writes the generated real inputs to:

- `/path/to/test_root/generated/sample_sheet.csv`
- `/path/to/test_root/generated/metadata.tsv`

If the source genomes are already cached and the server should avoid network
access, add `--offline`.

## 4. Make Sure Runtime Databases Are Ready

If the runtime databases are not already prepared on the server, build them
first:

```bash
nextflow run prepare_databases.nf -profile oist \
  --taxdump /shared/db/ncbi_taxdump_YYYYMMDD \
  --checkm2_db /shared/db/checkm2/CheckM2_database \
  --codetta_db /shared/db/codetta/Pfam-A_enone \
  --busco_db /shared/db/busco \
  --eggnog_db /shared/db/Eggnog_db/Eggnog_Diamond_db \
  --download_missing_databases true \
  --outdir /path/to/runtime-prep
```

Success criteria:

- the database destinations validate
- `runtime_database_report.tsv` is written under the prep `--outdir`

## 5. Run The Small Real Cohort Through The Acceptance Harness

Use the tracked small cohort on SLURM before any larger run:

```bash
python3 bin/run_acceptance_tests.py slurm \
  --annotation-config annotation.config \
  --work-root /path/to/test_root \
  --taxdump /shared/db/ncbi_taxdump_YYYYMMDD \
  --checkm2-db /shared/db/checkm2/CheckM2_database \
  --codetta-db /shared/db/codetta/Pfam-A_enone \
  --busco-db /shared/db/busco \
  --eggnog-db /shared/db/Eggnog_db/Eggnog_Diamond_db \
  --slurm-queue short \
  --singularity-cache-dir /path/to/singularity-cache
```

The harness uses `slurm,singularity` and analyses the complete declared cohort.
Pass `--annotation-config annotation.config` with the resources and immutable
runtimes described in [functional annotation](functional_annotation.md).

The harness writes the server run under:

- `/path/to/test_root/runs/slurm/results/`

## 6. Inspect The Outputs

Inspect these directories after the SLURM run:

- `/path/to/test_root/runs/slurm/results/tables/`
- `/path/to/test_root/runs/slurm/results/cohort/`
- `/path/to/test_root/runs/slurm/results/samples/`

Check the final tables:

- `master_table.tsv`
  - metadata columns are preserved and ordered correctly
  - Codetta fields are present
  - ANI fields are populated for eligible samples
- `sample_status.tsv`
  - the missing-metadata case is recorded
  - the atypical exclusion and atypical exception are recorded
  - BUSCO, CheckM2, Codetta, Prokka, CCFINDER, PADLOC, and eggNOG statuses are sensible
- `tool_and_db_versions.tsv`
  - database labels and paths are correct
  - Codetta version and source commit are present

Check the published sample folders:

- logs are present
- stable master-table input artefacts are present
- pruned raw directories still contain the retained files

## 7. Rerun With Resume

Rerun the same acceptance command with `--resume`:

```bash
python3 bin/run_acceptance_tests.py slurm \
  --annotation-config annotation.config \
  --work-root /path/to/test_root \
  --taxdump /shared/db/ncbi_taxdump_YYYYMMDD \
  --checkm2-db /shared/db/checkm2/CheckM2_database \
  --codetta-db /shared/db/codetta/Pfam-A_enone \
  --busco-db /shared/db/busco \
  --eggnog-db /shared/db/Eggnog_db/Eggnog_Diamond_db \
  --slurm-queue short \
  --singularity-cache-dir /path/to/singularity-cache \
  --resume
```

Success criteria:

- most or all tasks return as cached
- the final tables are unchanged

## Raw Nextflow Fallback

If you want to bypass the acceptance harness after the cohort is prepared, use
the generated inputs directly:

```bash
nextflow run . -c annotation.config -profile oist \
  --sample_csv /path/to/test_root/generated/sample_sheet.csv \
  --metadata /path/to/test_root/generated/metadata.tsv \
  --taxdump /shared/db/ncbi_taxdump_YYYYMMDD \
  --checkm2_db /shared/db/checkm2/CheckM2_database \
  --codetta_db /shared/db/codetta/Pfam-A_enone \
  --busco_db /shared/db/busco \
  --eggnog_db /shared/db/Eggnog_db/Eggnog_Diamond_db \
  --outdir /path/to/test_root/runs/raw-small/results \
  -work-dir /path/to/test_root/runs/raw-small/work
```

Use `-profile oist -c annotation.config` for the first raw small-cohort run. Add `-resume` only when reusing the exact same
launch directory and `-work-dir`.

Do not start the first small-cohort correctness run with
`conf/oist_20k_storage.config`. Use that opt-in override only after the normal
small-cohort path passes.

## Related Docs

- `docs/runbook.md`
- `docs/run_pipeline_test.md`
- `assets/testdata/acceptance/cohort_plan.tsv`
- `conf/oist_20k_storage.config`
