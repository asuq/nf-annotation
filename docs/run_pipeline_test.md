# Manual Pipeline Test Wrapper

## Purpose

`bin/run_pipeline_test.sh` is a thin operator-facing wrapper around
`bin/run_acceptance_tests.py`.

Use it when you want one stable command for manual pipeline testing without
remembering the lower-level Python entrypoint.

The wrapper does not re-implement pipeline logic. It forwards work to the
acceptance harness, which remains the single source of truth for validation,
Nextflow command construction, and output checks.

The wrapper now preflights the host `python3` interpreter before it calls the
acceptance harness. If `python3` on the login node is older than `3.12`, the
wrapper stops immediately with a direct operator message instead of failing
later with a raw Python `SyntaxError`.

For runtime database preparation on HPC, prefer the wrapper's `dbprep-slurm`
mode. The wrapper's `prepare` mode only prepares the tracked acceptance cohort
and its cached source genomes.

Use raw `nextflow run prepare_databases.nf` only when you need direct operator
control outside the standard SLURM validation gate.

The wrapper now preflights `dbprep-slurm` and `all` against
the root `nextflow.config` helper-image default. If the repo still pins the
stale helper image `quay.io/asuq1617/nf-myco_db:0.2`, the wrapper stops before
submitting work and asks you to move `params.runtime_db_helper_container` to
`quay.io/asuq1617/nf-myco_db:0.3`.

For a scripted SLURM validation of runtime database preparation itself, use the
wrapper's `dbprep-slurm` mode. That mode runs `prepare_databases.nf`, downloads
the curated runtime databases, and validates the prepared database tree.

For OIST or any other full-eggNOG HPC validation, use this sequence:

1. `bin/run_pipeline_test.sh prepare` on the login node
2. `bin/run_pipeline_test.sh dbprep-slurm` on SLURM
3. raw `nextflow run . -c annotation.config -profile oist` for the real pipeline run

Do not start with `all` on a fresh HPC environment. First make `prepare`,
then `dbprep-slurm`, then the tracked real pipeline run pass. Only use `all`
after those stages succeed.

For the broader OIST campaign through the medium real-data case, use the
dedicated wrapper:

```bash
bin/run_oist_hpc_matrix.sh --hpc-root /path/on/hpc/root all
```

That wrapper keeps using the coded OIST resource defaults. It does not add
`--max_cpus`, `--max_memory`, or `--max_time` overrides.
The default genetic-code rule is `mean_gene_length_ratio`: a paired ratio above
1.5 selects code 4; a valid ratio at or below 1.5 selects code 11. Both the
acceptance harness and HPC wrapper accept `--gcode-rule strict_delta` or
`--gcode-rule delta_then_11` to use the retained completeness-based criteria.
Codetta remains helper-prepared in the HPC matrix, so a directory with
`Pfam-A_enone.hmm` and the `.h3*` files but without `.nf_myco_ready.json` is
an expected `db-matrix` failure unless force rebuild is used.
The fixed medium `p2` validator rejects unresolved genetic-code assignments.
It allows an isolated failure in the secondary BUSCO lineage when the row
records `busco_summary_failed`. The primary BUSCO lineage remains strict.

The medium Mycoplasmatota/Bacillota cohort is now prepared in the same way as
the tracked small cohort from repo-tracked `source_catalog.tsv` and
`cohort_plan.tsv` files under `assets/testdata/medium/`.

The wrapper's normal `slurm` path remains centred on the acceptance
cohort and local-baseline comparison.

For the recommended first real server validation path with the tracked small
acceptance cohort, see [`docs/small_cohort_server_test.md`](small_cohort_server_test.md).

## Command

```bash
bin/run_pipeline_test.sh [--dry-run] <prepare|unit|stub|local|slurm|dbprep-slurm|all> [args...]
```

Supported modes:

- `prepare`: download or reuse the tracked acceptance source genomes
- `unit`: run the Python unit-test layer
- `stub`: run the stub smoke test through the acceptance harness
- `local`: run the real-data acceptance cohort with the local profile
- `slurm`: run the real-data acceptance cohort with the SLURM profile and compare against the latest local baseline
- `dbprep-slurm`: run the runtime database prep workflow on SLURM and validate the prepared database tree
- `all`: run `prepare`, `unit`, `stub`, `local`, and `slurm` in sequence

Wrapper options:

- `--dry-run`: print the delegated `python3 bin/run_acceptance_tests.py ...` command and exit
- `-h`, `--help`: show wrapper help before a mode is chosen

Arguments after the mode are forwarded unchanged to
`bin/run_acceptance_tests.py`.

If you need mode-specific help, run:

```bash
bin/run_pipeline_test.sh <mode> --help
```

## Prerequisites

Always required:

- `bash`
- `python3 >= 3.12`

Required by pipeline-running modes:

- `nextflow` for `stub`, `local`, `slurm`, `dbprep-slurm`, and `all`
- `sbatch` for `slurm`, `dbprep-slurm`, and `all`

Required by real-data modes:

- `--taxdump`
- `--checkm2-db`
- `--codetta-db`
- `--busco-db` or `--prepare-busco-datasets`
- `--eggnog-db`

CRISPRCasFinder still uses `params.ccfinder_container` from `nextflow.config`.
The wrapper does not accept a separate container override.

The wrapper does not prepare runtime databases such as taxdump, CheckM2,
Codetta, BUSCO, or eggNOG. Prepare those separately before real-data runs when
they are not already available. PADLOC uses the pinned database bundled in its
explicitly configured immutable image.

If an HPC login node still exposes an older interpreter such as `Python 3.6.8`,
load a newer module or move a Python `3.12` installation to the front of
`PATH` before using the wrapper:

```bash
which python3
python3 --version
module avail python
module load python/3.12
hash -r
python3 --version
```

Acceptance-backed `local`, `slurm`, and `all` runs use the normal execution
profiles and the complete declared cohort. Pass
`--annotation-config annotation.config` with the immutable tool runtimes and
resource paths documented in [functional annotation](functional_annotation.md).

## Examples

Quick stub smoke test:

```bash
bin/run_pipeline_test.sh stub
```

Prepare the cached real-data acceptance cohort:

```bash
bin/run_pipeline_test.sh prepare
```

Preview a real local run without executing it:

```bash
bin/run_pipeline_test.sh --dry-run local \
  --annotation-config annotation.config \
  --taxdump /path/to/pinned-taxdump \
  --checkm2-db /path/to/checkm2-db \
  --codetta-db /path/to/codetta-db \
  --busco-db /path/to/busco \
  --eggnog-db /path/to/eggnog-db
```

Run the local real-data acceptance cohort:

```bash
bin/run_pipeline_test.sh local \
  --annotation-config annotation.config \
  --taxdump /path/to/pinned-taxdump \
  --checkm2-db /path/to/checkm2-db \
  --codetta-db /path/to/codetta-db \
  --busco-db /path/to/busco \
  --eggnog-db /path/to/eggnog-db
```

Run the SLURM real-data acceptance cohort:

```bash
bin/run_pipeline_test.sh slurm \
  --annotation-config annotation.config \
  --taxdump /path/to/pinned-taxdump \
  --checkm2-db /path/to/checkm2-db \
  --codetta-db /path/to/codetta-db \
  --busco-db /path/to/busco \
  --eggnog-db /path/to/eggnog-db \
  --slurm-queue short
```

Run the SLURM runtime database-prep validation:

```bash
bin/run_pipeline_test.sh dbprep-slurm \
  --dbprep-profile oist \
  --work-root /path/to/work-root \
  --taxdump /path/to/db/ncbi_taxdump_20240914 \
  --checkm2-db /path/to/db/checkm2/CheckM2_database \
  --codetta-db /path/to/db/codetta/Pfam-A_enone \
  --busco-db /path/to/db/busco \
  --eggnog-db /path/to/db/Eggnog_db/Eggnog_Diamond_db \
  --slurm-queue short \
  --singularity-cache-dir /path/to/singularity-cache
```

Run the full layered workflow:

```bash
bin/run_pipeline_test.sh all \
  --annotation-config annotation.config \
  --taxdump /path/to/pinned-taxdump \
  --checkm2-db /path/to/checkm2-db \
  --codetta-db /path/to/codetta-db \
  --busco-db /path/to/busco \
  --eggnog-db /path/to/eggnog-db
```

## Output locations

The wrapper uses the acceptance harness defaults unless you override them.

- Stub smoke-test results are written under `assets/testdata/local/acceptance/runs/stub/`
- Generated acceptance inputs, downloads, and run artefacts are written under `assets/testdata/local/acceptance/`
- Final published tables for successful runs are under each run directory's `results/tables/`

For the lower-level interface and the raw Nextflow examples, see
`docs/runbook.md`.
