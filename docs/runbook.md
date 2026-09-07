# Runbook

## Required inputs

The pipeline requires:

- `--sample_csv`: manifest CSV with the locked lower-case headers
- `--metadata`: metadata table with `Accession` or `accession`
- `--taxdump`: pinned NCBI taxdump directory containing `names.dmp` and `nodes.dmp`
- `--checkm2_db`: CheckM2 database directory containing one top-level `.dmnd`
- `--codetta_db`: Codetta profile directory containing `Pfam-A_enone.hmm` and
  its `.h3f`, `.h3i`, `.h3m`, and `.h3p` sidecars
- `--busco_db` or `--prepare_busco_datasets true`
- `--eggnog_db` for real eggNOG runs

## Nextflow version

The pipeline requires Nextflow `>=26.04.0`, matching the strict syntax parser
used by current stable Nextflow releases. Select a v26 runtime explicitly on
systems where multiple versions are available:

```bash
NXF_VER=26.04.0 nextflow run . -profile test -stub-run
```

## BUSCO lineage datasets

`main.nf` reuses existing lineage directories below `--busco_db` by default.
For example, a configured lineage `alcaligenaceae_odb12` is expected at:

```text
${BUSCO_DIR}/alcaligenaceae_odb12
```

If a configured lineage is missing, either prepare the BUSCO root ahead of the
analysis run with the same lineage list:

```bash
nextflow run prepare_databases.nf -profile oist \
  --busco_db "$BUSCO_DIR" \
  --busco_lineages alcaligenaceae_odb12 \
  --download_missing_databases true \
  --singularity_cache_dir "$SINGULARITY_CACHE" \
  --outdir "$DBPREP_OUT"
```

or let the main workflow download missing BUSCO lineages explicitly:

```bash
nextflow run . -profile oist \
  --sample_csv "$SAMPLE_CSV" \
  --metadata "$METADATA_TSV" \
  --taxdump "$TAXDUMP_DIR" \
  --checkm2_db "$CHECKM2_DIR" \
  --codetta_db "$CODETTA_DIR" \
  --busco_db "$BUSCO_DIR" \
  --busco_lineages alcaligenaceae_odb12 \
  --prepare_busco_datasets true \
  --eggnog_db "$EGGNOG_DIR" \
  --singularity_cache_dir "$SINGULARITY_CACHE" \
  --outdir "$RESULT_ROOT/out"
```

When `--prepare_busco_datasets true` is enabled, existing lineage directories
are reused and only missing configured lineages are downloaded under
`--busco_db`. The main downloader stages the BUSCO DB parent as a Nextflow
`path` input so Singularity/Apptainer profiles can bind the shared destination
for writes.

### Preparing metadata.tsv

`accessions.txt` should contain one accession per line. `datasets` and
`dataformat` must be installed and available on `PATH`. This prepares only
`metadata.tsv`; the sample manifest still needs to be prepared separately.
This recipe was checked locally in the `gtdb-genome` environment. The
resulting `metadata.raw.tsv` uses human-readable `dataformat` headers, not the
`--fields` token names.

```bash
datasets summary genome accession \
  --inputfile accessions.txt \
  --assembly-source GenBank \
  --assembly-version latest \
  --as-json-lines | \
  dataformat tsv genome \
    --fields accession,assminfo-name,organism-name,assminfo-atypicalwarnings,assminfo-notes,assmstats-total-sequence-len,type_material-label,organism-tax-id,assminfo-level,assminfo-release-date,assminfo-sequencing-tech,assminfo-assembly-method,assmstats-number-of-contigs,assmstats-number-of-scaffolds,assmstats-scaffold-n50,assmstats-genome-coverage,source_database,assminfo-bioproject,assminfo-biosample-accession,assminfo-submitter \
  > metadata.raw.tsv
```

```bash
python3 - <<'PY'
import csv
from pathlib import Path

input_path = Path("metadata.raw.tsv")
output_path = Path("metadata.tsv")

column_map = [
    ("Assembly Accession", "accession"),
    ("Assembly Name", "Assembly_Name"),
    ("Organism Name", "Organism_Name"),
    ("Assembly Atypical Warnings", "Atypical_Warnings"),
    ("Assembly Notes", "Notes"),
    ("Assembly Stats Total Sequence Length", "Genome_Size"),
    ("Type Material Label", "Type_Material"),
    ("Organism Taxonomic ID", "Tax_ID"),
    ("Assembly Level", "Assembly_Level"),
    ("Assembly Release Date", "Release_Date"),
    ("Assembly Sequencing Tech", "Sequencing_Tech"),
    ("Assembly Assembly Method", "Assembly_Method"),
    ("Assembly Stats Number of Contigs", "Contigs"),
    ("Assembly Stats Number of Scaffolds", "Scaffolds"),
    ("Assembly Stats Scaffold N50", "N50"),
    ("Assembly Stats Genome Coverage", "Genome_Coverage"),
    ("Source Database", "Source_Database"),
    ("Assembly BioProject Accession", "BioProject"),
    ("Assembly BioSample Accession", "BioSample"),
    ("Assembly Submitter", "Submitter"),
]

with input_path.open("r", encoding="utf-8", newline="") as handle:
    reader = csv.reader(handle, delimiter="\t")
    rows = list(reader)

if not rows:
    raise SystemExit("metadata.raw.tsv is empty.")

header = [cell.strip() for cell in rows[0]]
index = {name: position for position, name in enumerate(header)}
missing = [source for source, _target in column_map if source not in index]
if missing:
    raise SystemExit(
        "metadata.raw.tsv is missing expected column(s): " + ", ".join(missing)
    )

with output_path.open("w", encoding="utf-8", newline="") as out_handle:
    writer = csv.writer(out_handle, delimiter="\t")
    writer.writerow([target for _source, target in column_map])
    for row in rows[1:]:
        writer.writerow(
            [
                row[index[source]] if index[source] < len(row) else ""
                for source, _target in column_map
            ]
        )
PY
```

The rename step matters because downstream consumers currently look for
headers such as `Tax_ID`, `Organism_Name`, `Assembly_Level`,
`Atypical_Warnings`, `Notes`, `Genome_Size`, `Scaffolds`, and `N50`.

Optional labels used only in `tool_and_db_versions.tsv`:

- `--taxdump_label`
- `--checkm2_db_label`
- `--codetta_db_label`
- `--eggnog_db_label`

Shared helper image requirement:

- `params.python_container` is the single helper image for Python-based tasks
  such as validation, ANI clustering, ANI representative selection, final
  status assembly, and version collection
- it must provide `biopython`, `numpy`, and `scipy`
- the repo-owned Dockerfile for this image lives under
  `docker/python_helper/Dockerfile`
- local Docker and HPC Singularity runs should reuse that one helper image to
  reduce container pulls

Runtime database prep helper image:

- `prepare_databases.nf` uses the dedicated helper image configured by
  `params.runtime_db_helper_container`
- the default helper image is `quay.io/asuq1617/nf-myco_db:0.3`
- the repo-owned Dockerfile for this image lives under
  `docker/runtime_db_helper/Dockerfile`

## Runtime database preparation

Use `nextflow run prepare_databases.nf` when one or more runtime databases are
missing on the target machine.

- It is a separate operator-facing Nextflow entry point, not part of the normal
  analysis workflow runtime.
- It uses the same database directory params as the main workflow.
- An existing valid database directory is reused in place.
- A missing database directory is populated there when
  `--download_missing_databases true` is set.
- Taxdump uses the shared helper downloader. CheckM2, BUSCO, and eggNOG use
  their own tool-native download or update commands.
- Scratch work uses `--runtime_db_scratch_root` when set and otherwise stays
  beside the destination rather than in system tmp.
- A prepared destination is reusable only when it validates and contains
  `.nf_myco_ready.json`.
- Concurrent preparation is blocked by a per-destination lock sidecar ending in
  `.nf_myco_prepare.lock`.
- The helper logic lives in `bin/prepare_runtime_databases.py`, but that script
  is now the implementation detail behind the prep workflow rather than the
  main operator interface.

Supported runtime databases:

- taxdump
- CheckM2
- Codetta
- BUSCO lineage directories
- eggNOG

Example:

```bash
nextflow run prepare_databases.nf -profile oist \
  --taxdump /shared/db/ncbi_taxdump_20240914 \
  --checkm2_db /shared/db/checkm2/CheckM2_database \
  --codetta_db /shared/db/codetta/Pfam-A_enone \
  --busco_db /shared/db/busco \
  --eggnog_db /shared/db/Eggnog_db/Eggnog_Diamond_db \
  --download_missing_databases true \
  --runtime_db_scratch_root /shared/db/.scratch \
  --outdir /shared/db/runtime-prep
```

The prep workflow writes `runtime_database_report.tsv` and `nextflow_args.txt`
under `--outdir` on success.

## Profiles

- `debug`: composable behaviour profile that defaults eggNOG smoke runs to `GCA_000027325.1`
- `local`: local executor
- `slurm`: SLURM executor with optional `params.slurm_queue`, `params.slurm_qos`, and `params.slurm_cluster_options`
- `singularity`: Singularity execution with optional `params.singularity_cache_dir` and `params.singularity_run_options`
- `oist`: standalone OIST HPC profile with SLURM and Singularity enabled, using the submitting user account
- `gwdg`: standalone GWDG SCC profile with SLURM, Singularity, `LOCAL_TMPDIR` scratch, and SHM-first temp files
- `test`: local fixture profile for `-stub-run`

## Minimal test path

Use the frozen fixtures under `assets/testdata/stub/`:

```bash
nextflow run . -profile test -stub-run
```

This validates the DSL wiring, channel contracts, publish locations, and final
table/version outputs without requiring external databases or tool containers.

## Updating cohort membership

Use the main workflow's `--update_from <published-results-directory>` option
for routine sample additions and removals. The sample CSV is the complete
desired cohort, not a list of changes. An accession is reused only when it is
in the source's `tables/validated_samples.tsv`; the `is_new` CSV flag does not
control reuse. An accession absent from the chosen source is analysed again,
even if it appeared in some other, older results directory.

For example, after a normal run containing A and B, give the update a sample
CSV containing B and C. B's published QC and annotation are reused, C receives
normal analysis, and ANI and cohort reports are rebuilt for B and C. A's files
remain in the previous output and are not copied into the new cohort.

```bash
nextflow run . -profile local,docker \
  --update_from /path/to/results-v1 \
  --sample_csv samples-v2.csv \
  --metadata metadata-v2.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --outdir /path/to/results-v2
```

Use the usual execution profile on HPC. Per-sample database paths and BUSCO
dataset preparation are needed only if the revised cohort contains additions.
For a removal-only update:

```bash
nextflow run . -profile local,docker \
  --update_from /path/to/results-v2 \
  --sample_csv samples-v3.csv \
  --metadata metadata-v3.tsv \
  --taxdump /path/to/pinned-taxdump \
  --outdir /path/to/results-v3
```

Supply the intended `--busco_lineages` in order, plus any non-default ANI
settings, on every update. Requested BUSCO lineages must already be available
for retained samples. Their QC/annotation settings and outcomes remain as
published; newly supplied per-sample tool settings govern additions only.
The ANI threshold, scoring profile, primary BUSCO column and optional
incomplete-16S gate apply to the whole revised cohort. ANI pairs are recomputed,
including retained-to-retained comparisons; this is not a pairwise ANI cache.

Before new analysis starts, preflight requires the source validated manifest,
complete master/status tables, versions report and retained sample artefacts.
Documented failed/skipped tool outcomes remain valid stored outcomes. Missing
files, malformed or inconsistent summaries, nested sample symlinks, and
changed genomes under a retained accession cause an explicit error. Older or
partial results lacking the required published format must be repaired before
they can be used for a routine update; the recovery helper below has a
different purpose.

The updated CSV must resolve each genome, including retained genomes. A
retained input can be moved, rewrapped or gzip-compressed, and its FASTA
description may change. Ordered record IDs and sequence bytes, including
case/soft masking, must match the source's published staged FASTA. When the
original FASTA is unavailable, point the CSV at that staged FASTA instead.
Use a new accession or a full analysis for a changed assembly. At least one
sample must remain in the requested cohort; zero ANI-eligible samples are
supported and produce header-only ANI tables with explicit exclusion reasons.

The new results directory must be separate from the source: neither directory
may contain the other. Retained per-sample artefacts are copied, so a completed
update works as a later source without the original work directory, cache or
older results tree. Existing internal IDs are preserved when membership
changes create or remove sanitisation collisions.

Additional audit outputs under `tables/` are:

- `cohort_update.tsv`: added/reused/removed actions, internal IDs, fingerprints
  of requested genomes, source location and source manifest/versions hashes.
  Removed genomes have fingerprint `NA` because their files are not required.
- `cohort_update_run.json`: the input, source and settings identity used to
  prevent resuming an unrelated update into the same directory.
- `inherited_versions.tsv`: source provenance with explicit reuse labels.
  These rows are also retained in `tool_and_db_versions.tsv` alongside the
  current execution's provenance; reused tools are not reported as rerun.

The existing master/status table columns stay unchanged. Metadata and
cohort-dependent fields are rebuilt from current inputs, so ANI cluster IDs,
representatives and inclusion decisions can change. No stale cohort tables
are copied from the source.

To recover an interrupted update, rerun its exact command with `-resume`,
preserving that update's launch/cache and work directories. Inputs and
scientific settings must match its saved identity. Changing the cohort again
requires another fresh `--outdir` and a completed published source.

## Rescue ANI from published results

Use `bin/rescue_ani_from_results.py` when an older run already published
per-sample Barrnap, CheckM2, and BUSCO outputs, but the ANI branch needs to be
rebuilt without rerunning the heavy per-sample tools.

- `--source-outdir` must point to the failed run's published results tree
- `--outdir` must be a different fresh rescue directory
- by default, rescue rebuilds the validation-time `source_sample_status.tsv`
  seed from `validated_samples.tsv` plus `validation_warnings.tsv`
- the old published `tables/sample_status.tsv` is not used as the default seed,
  because in a completed or partially completed run it may already be the later
  overlaid final audit table
- `--initial-status` remains available as an explicit override when you really
  have the original validation-time seed table
- `--busco-lineage` must be supplied once with the original run order as
  whitespace-separated values because the first lineage remains the primary
  ANI-scoring BUSCO column

Example:

```bash
uv run bin/rescue_ani_from_results.py \
  --source-outdir /path/to/failed-results \
  --metadata /path/to/metadata.tsv \
  --outdir /path/to/rescued-ani \
  --busco-lineage bacillota_odb12 mycoplasmatota_odb12
```

If FastANI is only available through a container wrapper, pass the wrapper as a
quoted command prefix:

```bash
uv run bin/rescue_ani_from_results.py \
  --source-outdir /path/to/failed-results \
  --metadata /path/to/metadata.tsv \
  --outdir /path/to/rescued-ani \
  --busco-lineage bacillota_odb12 mycoplasmatota_odb12 \
  --fastani-binary "singularity exec /path/to/fastani.sif fastANI"
```

By default, rescue reuses the published source table at
`cohort/assembly_stats/assembly_stats.tsv`. If that table is missing or you do
not trust it, point rescue at a different TSV with `--assembly-stats` or opt
into recomputation with `--recalculate-assembly-stats`.

If `seqtk` is only available through a container wrapper, pass that wrapper
when recomputation is requested:

```bash
uv run bin/rescue_ani_from_results.py \
  --source-outdir /path/to/failed-results \
  --metadata /path/to/metadata.tsv \
  --outdir /path/to/rescued-ani \
  --busco-lineage bacillota_odb12 mycoplasmatota_odb12 \
  --recalculate-assembly-stats \
  --seqtk-binary "singularity exec /path/to/seqtk.sif seqtk" \
  --fastani-binary "singularity exec /path/to/fastani.sif fastANI"
```

The rescue command rebuilds:

- per-sample `best_16S.fna` and `16S_status.tsv`
- per-sample `checkm2_summary.tsv`
- per-lineage BUSCO summary TSVs
- `cohort/assembly_stats/assembly_stats.tsv` by copying the published source table, or by recomputing it when `--recalculate-assembly-stats` is used
- `cohort/fastani/` ANI inputs, exclusions, matrix, and logs
- `cohort/ani_clusters/` cluster, ANI summary, and representative tables
- partial `tables/master_table.tsv` and `tables/sample_status.tsv`
- `tables/rescue_provenance.tsv`

## Small-cohort server validation

For the recommended first real server validation path, use the tracked small
acceptance cohort instead of starting with a custom or medium cohort. The full
step-by-step guide is in
[`docs/small_cohort_server_test.md`](small_cohort_server_test.md).

In short:

1. activate the `nextflow` conda environment on the login node
2. run the stub smoke test
3. prepare the tracked small acceptance cohort
4. ensure the runtime databases are ready
5. run the SLURM acceptance cohort
6. inspect the published outputs
7. rerun with `--resume`

Do not start the first small-cohort correctness run with
`conf/oist_20k_storage.config`. Use the normal path first, then test the 20k
override separately.

## Manual wrapper

For a single operator-facing entrypoint, use `bin/run_pipeline_test.sh`.

See `docs/run_pipeline_test.md` for wrapper usage, prerequisites, examples, and
expected output locations. Its `dbprep-slurm` and `all` modes now preflight the
repo config so they fail fast if the root helper-image default still pins the
stale runtime-db helper image tag.

## Acceptance harness

Use `bin/run_acceptance_tests.py` for the layered acceptance workflow:

- `prepare`: download or reuse the tracked acceptance source genomes and build a generated cohort under `assets/testdata/local/acceptance/`
- `unit`: run the Python unit-test layer for fine-grained and minor edge cases
- `stub`: run the full-pipeline `-stub-run` smoke test
- `local`: run the generated positive cohort with `-profile debug,local,docker`
- `slurm`: run the same cohort with `-profile debug,slurm,singularity` and compare its stable outputs against the latest successful local run
- `dbprep-slurm`: run `prepare_databases.nf` on SLURM, download the curated runtime databases, and validate the prepared database tree
- `all`: run `prepare`, `unit`, `stub`, `local`, and `slurm` in sequence

Tracked cohort descriptors live under `assets/testdata/acceptance/`. Large
downloads, generated manifests, and run artefacts stay under the ignored local
tree `assets/testdata/local/acceptance/`.

Real-data `local`, `slurm`, and `all` runs require:

- `--taxdump`
- `--checkm2-db`
- `--codetta-db`
- `--busco-db` or `--prepare-busco-datasets`
- `--eggnog-db`

These acceptance runs also assume the CRISPRCasFinder image is already set via
`params.ccfinder_container` in pipeline config. The harness does not accept a
separate CCFINDER container argument.

PADLOC uses the fixed database bundled in the default PADLOC image, so neither
the pipeline nor the acceptance harness accepts an external PADLOC database
path anymore.

Acceptance runs use the `debug` profile by default, which sets
`params.eggnog_only_accessions = 'GCA_000027325.1'`. Normal raw pipeline runs
still execute eggNOG for every gcode-qualified sample unless you opt into
`debug` or set that parameter explicitly. If you override `--local-profile` or
`--slurm-profile` in the harness, include `debug` yourself if you still want
the smoke-only eggNOG behaviour.

For OIST or any other full-eggNOG HPC validation, use raw `nextflow run .`
instead of the default SLURM acceptance wrapper. The wrapper's SLURM mode is
still centred on the debug acceptance cohort and local-baseline comparison.

Example local acceptance run:

```bash
python3 bin/run_acceptance_tests.py local \
  --taxdump /path/to/pinned-taxdump \
  --checkm2-db /path/to/checkm2-db \
  --codetta-db /path/to/codetta-db \
  --busco-db /path/to/busco \
  --eggnog-db /path/to/eggnog-db
```

Example SLURM acceptance run:

```bash
python3 bin/run_acceptance_tests.py slurm \
  --taxdump /path/to/pinned-taxdump \
  --checkm2-db /path/to/checkm2-db \
  --codetta-db /path/to/codetta-db \
  --busco-db /path/to/busco \
  --eggnog-db /path/to/eggnog-db \
  --slurm-queue short \
  --slurm-qos 2h
```

Example SLURM database-prep run:

```bash
python3 bin/run_acceptance_tests.py dbprep-slurm \
  --dbprep-profile oist \
  --work-root /path/to/work-root \
  --taxdump /path/to/db/ncbi_taxdump_20240914 \
  --checkm2-db /path/to/db/checkm2/CheckM2_database \
  --codetta-db /path/to/db/codetta/Pfam-A_enone \
  --busco-db /path/to/db/busco \
  --eggnog-db /path/to/db/Eggnog_db/Eggnog_Diamond_db \
  --slurm-queue short \
  --slurm-qos 2h \
  --singularity-cache-dir /path/to/singularity-cache
```

## Example real runs

Local:

```bash
nextflow run . \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --outdir results
```

Local debug smoke variant:

```bash
nextflow run . -profile debug,local,docker \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --outdir results
```

Local debug run with an overridden eggNOG smoke accession:

```bash
nextflow run . -profile debug,local,docker \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --eggnog_only_accessions SOME_OTHER_ACCESSION \
  --outdir results
```

SLURM:

```bash
nextflow run . -profile slurm \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --slurm_qos 2h \
  --outdir results
```

Use `--slurm_qos 2h` for QoS. If you need generic SLURM options that begin
with a dash, use equals syntax so Nextflow keeps the value attached, for
example `--slurm_cluster_options='--qos=2h'`.

Singularity:

```bash
nextflow run . -profile singularity \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --outdir results
```

OIST:

```bash
nextflow run . -profile oist \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --singularity_cache_dir /path/to/singularity-cache \
  --outdir results
```

GWDG SCC:

```bash
nextflow run . -profile gwdg \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --singularity_cache_dir /path/to/singularity-cache \
  --outdir results
```

## OIST HPC testing

Use the OIST profile directly for full eggNOG runs:

```bash
nextflow run . -profile oist ...
```

Do not include `debug` and do not set `--eggnog_only_accessions` when the goal
is to validate full eggNOG execution on HPC.

Build everything on HPC-native storage. Do not copy the local database tree to
HPC. Use three stages in order:

1. prepare the tracked acceptance cohort on the HPC login node
2. prepare runtime databases on SLURM with `dbprep-slurm`
3. run the real pipeline with raw `nextflow run . -profile oist`

For the scripted campaign through the medium Mycoplasmatota/Bacillota run, use
the dedicated wrapper:

```bash
bin/run_oist_hpc_matrix.sh --hpc-root /path/on/hpc/root all
```

Optional:

```bash
bin/run_oist_hpc_matrix.sh --hpc-root /path/on/hpc/root --gcode-rule delta_then_11 all
```

That wrapper:

- prepares the tracked 9-sample cohort on HPC-native storage
- prepares the fixed medium Mycoplasmatota/Bacillota cohort from the
  repo-tracked source catalog and cohort plan under `assets/testdata/medium/`
- runs the full runtime DB gate and the disposable DB existence-state cases
- runs the tracked 9-sample real pipeline gate
- runs the medium Mycoplasmatota/Bacillota real-data gate
- writes the medium run `sample_sheet.csv`, `metadata.tsv`, and
  `source_stats.tsv` under `"$HPC_ROOT/medium/generated"`
- keeps Codetta on the helper-prepared path rather than the native
  `download -> finalise` path, so a Codetta directory with valid profile files
  but without `.nf_myco_ready.json` is treated as incomplete in `db-matrix`
  unless force rebuild is requested

It intentionally uses the resource defaults already coded in the OIST profile
and process configuration. It does not pass explicit `--max_cpus`,
`--max_memory`, or `--max_time` overrides.
Nextflow runtime tasks retry up to 3 total attempts. Database download tasks
retry up to 2 total attempts. Sample-level soft-fail modules also retry their
wrapped tools internally up to `soft_fail_attempts` total attempts before
emitting fallback placeholder outputs.

Step-by-step procedure:

1. Start a persistent shell on the login node and move into the repository.

```bash
tmux new -s nf_annotation_hpc
cd /path/to/nf-annotation
```

2. Define one clean HPC root.

```bash
export HPC_ROOT=/scratch/$USER/nf_annotation_hpc_$(date +%Y%m%d)
export ACCEPT_ROOT=$HPC_ROOT/acceptance
export DB_ROOT=$HPC_ROOT/db
export RESULT_ROOT=$HPC_ROOT/results
export SINGULARITY_CACHE=$HPC_ROOT/singularity_cache

mkdir -p "$ACCEPT_ROOT" "$DB_ROOT" "$RESULT_ROOT" "$SINGULARITY_CACHE"
```

3. Run preflight checks.

```bash
which python3
python3 --version
nextflow -version
singularity --version
sbatch --version
nextflow config -profile oist >/dev/null
```

The host-side wrappers and validators require `python3 >= 3.12`. If the login
node still exposes an older interpreter such as `Python 3.6.8`, load a newer
module or adjust `PATH` before you use the wrappers:

```bash
module avail python
module load python/3.12
hash -r
which python3
python3 --version
```

If the login node cannot provide `python3 >= 3.12`, skip the wrappers and run
raw Nextflow directly instead. In that fallback path you lose the wrapper-level
acceptance checks, but the actual pipeline tasks still run inside containers.

4. Prepare the tracked 9-sample cohort on the HPC login node.

```bash
bin/run_pipeline_test.sh prepare \
  --work-root "$ACCEPT_ROOT"
```

Expected result:

- `$ACCEPT_ROOT/generated/sample_sheet.csv`
- `$ACCEPT_ROOT/generated/metadata.tsv`
- `$ACCEPT_ROOT/download_checksums.tsv`
- downloaded FASTA files under `$ACCEPT_ROOT/downloads/`

If this step fails, bring back:

- `$ACCEPT_ROOT/prepare.log` if present
- `$ACCEPT_ROOT/generated/`
- `$ACCEPT_ROOT/download_checksums.tsv`

5. Prepare the runtime databases on SLURM with the wrapper gate.

Set the explicit HPC destination directories first:

```bash
export TAXDUMP_DIR=$DB_ROOT/ncbi_taxdump_20240914
export CHECKM2_DIR=$DB_ROOT/checkm2/CheckM2_database
export CODETTA_DIR=$DB_ROOT/codetta/Pfam-A_enone
export BUSCO_DIR=$DB_ROOT/busco
export EGGNOG_DIR=$DB_ROOT/Eggnog_db/Eggnog_Diamond_db
```

```bash
bin/run_pipeline_test.sh dbprep-slurm \
  --dbprep-profile oist \
  --work-root "$ACCEPT_ROOT" \
  --taxdump "$TAXDUMP_DIR" \
  --checkm2-db "$CHECKM2_DIR" \
  --codetta-db "$CODETTA_DIR" \
  --busco-db "$BUSCO_DIR" \
  --eggnog-db "$EGGNOG_DIR" \
  --singularity-cache-dir "$SINGULARITY_CACHE"
```

Expected result:

- `$ACCEPT_ROOT/runs/dbprep-slurm/results/runtime_database_report.tsv`
- `$ACCEPT_ROOT/runs/dbprep-slurm/results/nextflow_args.txt`
- `$ACCEPT_ROOT/runs/dbprep-slurm/results/pipeline_info/trace.tsv`
- `$ACCEPT_ROOT/runs/dbprep-slurm/results/pipeline_info/report.html`
- `$ACCEPT_ROOT/runs/dbprep-slurm/results/pipeline_info/timeline.html`
- `$ACCEPT_ROOT/runs/dbprep-slurm/results/pipeline_info/dag.html`

Validate the prepared DB directories:

- taxdump:
  - `$TAXDUMP_DIR/names.dmp`
  - `$TAXDUMP_DIR/nodes.dmp`
- CheckM2:
  - exactly one top-level `*.dmnd` in `$CHECKM2_DIR`
- Codetta:
  - `$CODETTA_DIR/Pfam-A_enone.hmm`
  - `$CODETTA_DIR/Pfam-A_enone.hmm.h3f`
  - `$CODETTA_DIR/Pfam-A_enone.hmm.h3i`
  - `$CODETTA_DIR/Pfam-A_enone.hmm.h3m`
  - `$CODETTA_DIR/Pfam-A_enone.hmm.h3p`
- BUSCO:
  - `$BUSCO_DIR/bacillota_odb12/dataset.cfg`
  - `$BUSCO_DIR/mycoplasmatota_odb12/dataset.cfg`
- eggNOG:
  - `$EGGNOG_DIR/eggnog.db`
  - `$EGGNOG_DIR/eggnog_proteins.dmnd`
- ready markers:
  - `.nf_myco_ready.json` in each prepared DB root

For Codetta specifically, the ready marker is required for reuse. A directory
that already contains `Pfam-A_enone.hmm` and the `.h3*` sidecars but lacks
`.nf_myco_ready.json` is intentionally treated as incomplete by the helper and
is an expected failure in the OIST `db-matrix` case until force rebuild is
used.

If this step fails, bring back:

- `$ACCEPT_ROOT/runs/dbprep-slurm/results/`
- failed task dir under `$ACCEPT_ROOT/runs/dbprep-slurm/work/`
- `.nextflow.log`
- `.command.sh`, `.command.out`, `.command.err`, and `.exitcode` from the failed task

Fallback when host Python cannot be fixed:

```bash
nextflow run prepare_databases.nf -profile oist \
  --taxdump "$TAXDUMP_DIR" \
  --checkm2_db "$CHECKM2_DIR" \
  --codetta_db "$CODETTA_DIR" \
  --busco_db "$BUSCO_DIR" \
  --eggnog_db "$EGGNOG_DIR" \
  --download_missing_databases true \
  --runtime_db_scratch_root "$DB_ROOT/.scratch" \
  --singularity_cache_dir "$SINGULARITY_CACHE" \
  --outdir "$ACCEPT_ROOT/dbprep_manual"
```

That path still prepares the databases, but it skips the wrapper-level
acceptance validation.

6. Run one optional structural smoke test.

```bash
nextflow run . -profile test -stub-run --outdir "$RESULT_ROOT/stub"
```

7. Run the tracked 9-sample cohort on OIST with full eggNOG.

Do not jump straight to `all` on a new HPC setup. Make this tracked run pass
first, then move on to the full scripted campaign.

```bash
nextflow run . -profile oist \
  -work-dir "$RESULT_ROOT/p1/work" \
  --sample_csv "$ACCEPT_ROOT/generated/sample_sheet.csv" \
  --metadata "$ACCEPT_ROOT/generated/metadata.tsv" \
  --taxdump "$TAXDUMP_DIR" \
  --checkm2_db "$CHECKM2_DIR" \
  --codetta_db "$CODETTA_DIR" \
  --busco_db "$BUSCO_DIR" \
  --eggnog_db "$EGGNOG_DIR" \
  --singularity_cache_dir "$SINGULARITY_CACHE" \
  --outdir "$RESULT_ROOT/p1/out"
```

Do not include `debug` and do not set `--eggnog_only_accessions`.

For large OIST cohorts where shared-storage pressure matters, add the tracked
opt-in storage override:

```bash
nextflow run . -profile oist \
  -c conf/oist_20k_storage.config \
  -work-dir /flash/path/to/work_nf_annotation_20k \
  --sample_csv "$ACCEPT_ROOT/generated/sample_sheet.csv" \
  --metadata "$ACCEPT_ROOT/generated/metadata.tsv" \
  --taxdump "$TAXDUMP_DIR" \
  --checkm2_db "$CHECKM2_DIR" \
  --codetta_db "$CODETTA_DIR" \
  --busco_db "$BUSCO_DIR" \
  --eggnog_db "$EGGNOG_DIR" \
  --singularity_cache_dir "$SINGULARITY_CACHE" \
  --outdir "$RESULT_ROOT/p1/out"
```

To guarantee `-resume` for that run design:

- keep the same launch directory so `.nextflow/cache` is reused
- keep the same `-work-dir` path
- do not clean the persistent `/flash` work tree between runs
- do not enable `cleanup = true`
- do not change input paths, mtimes, or task-defining config between runs

Resume example:

```bash
nextflow run . -profile oist \
  -c conf/oist_20k_storage.config \
  -resume \
  -work-dir /flash/path/to/work_nf_annotation_20k \
  --sample_csv "$ACCEPT_ROOT/generated/sample_sheet.csv" \
  --metadata "$ACCEPT_ROOT/generated/metadata.tsv" \
  --taxdump "$TAXDUMP_DIR" \
  --checkm2_db "$CHECKM2_DIR" \
  --codetta_db "$CODETTA_DIR" \
  --busco_db "$BUSCO_DIR" \
  --eggnog_db "$EGGNOG_DIR" \
  --singularity_cache_dir "$SINGULARITY_CACHE" \
  --outdir "$RESULT_ROOT/p1/out"
```

Monitor:

```bash
squeue -u "$USER"
tail -f .nextflow.log
```

Expected final outputs:

- `$RESULT_ROOT/p1/out/tables/master_table.tsv`
- `$RESULT_ROOT/p1/out/tables/sample_status.tsv`
- `$RESULT_ROOT/p1/out/tables/tool_and_db_versions.tsv`
- `$RESULT_ROOT/p1/out/pipeline_info/trace.tsv`
- `$RESULT_ROOT/p1/out/pipeline_info/report.html`
- `$RESULT_ROOT/p1/out/pipeline_info/timeline.html`
- `$RESULT_ROOT/p1/out/pipeline_info/dag.html`

Acceptance checks:

- `sample_status.tsv` has no failed samples
- `tool_and_db_versions.tsv` has no `unknown` rows
- `eggnog_mapper` is `2.1.13`
- PADLOC tool version is clean
- DB path rows point at `$DB_ROOT/...`
- for full eggNOG runs, `eggnog_status` is `done` for eligible samples and
  `skipped` only when `gcode = NA`

8. Prepare and run the fixed medium Mycoplasmatota/Bacillota cohort.

The wrapper prepares this cohort automatically from the repo-tracked medium
catalogue and plan before `p2` and `all`, or you can run that step directly:

```bash
bin/run_oist_hpc_matrix.sh --hpc-root "$HPC_ROOT" medium-prepare
```

The generated medium inputs are written under:

- `"$HPC_ROOT/medium/generated/sample_sheet.csv"`
- `"$HPC_ROOT/medium/generated/metadata.tsv"`
- `"$HPC_ROOT/medium/generated/source_stats.tsv"`
- `"$HPC_ROOT/medium/download_checksums.tsv"`

Then run the medium case with:

```bash
bin/run_oist_hpc_matrix.sh --hpc-root "$HPC_ROOT" p2
```

The fixed medium cohort intentionally includes Bacillota genomes that may stay
unresolved at gcode assignment under the default `strict_delta` rule. For `p2`,
the validator accepts those rows only when `gcode_status=failed` is the sole
failed column, `warnings` includes `gcode_na`, and downstream annotation
statuses remain non-failed, typically `skipped`.
The same `p2` validator also accepts an isolated failure in the secondary
BUSCO lineage when `warnings` includes `busco_summary_failed`. Primary BUSCO
lineage failures still fail validation because they affect ANI eligibility and
indicate a materially incomplete QC result.

Only after `prepare`, `dbprep-slurm`, and the tracked `p1` run succeed should
you launch the full scripted campaign:

```bash
bin/run_oist_hpc_matrix.sh --hpc-root "$HPC_ROOT" all
```

9. If a run fails, download the useful artefacts back to local.

For database-prep failures, bring back:

- `$ACCEPT_ROOT/runs/dbprep-slurm/results/`
- failed task dir under `$ACCEPT_ROOT/runs/dbprep-slurm/work/`
- `.nextflow.log`

For pipeline failures, bring back:

- `$RESULT_ROOT/<case>/out/tables/`
- `$RESULT_ROOT/<case>/out/pipeline_info/`
- failed task dir under `$RESULT_ROOT/<case>/work/`
- `.nextflow.log`

Minimum useful artefacts:

- `.command.sh`
- `.command.out`
- `.command.err`
- `.exitcode`
- `versions.yml`
- task-specific logs such as `fastani.log`, `ccfinder.log`, or `result.json`

## Final outputs

Published final tables are written under `results/tables/`:

- `master_table.tsv`
- `sample_status.tsv`
- `tool_and_db_versions.tsv`

The master table preserves the original metadata block, then appends derived
columns in the runtime order resolved from `params.busco_lineages`. The
checked-in `assets/master_table_append_columns.txt` file records the default
lineage pair only. The appended block starts with `is_new`, sourced from
`validated_samples.tsv`. `GC_Content` is derived from the published
`cohort/assembly_stats/assembly_stats.tsv` table produced by `seqtk comp`.
Codetta contributes `Codetta_Genetic_Code` and
`Codetta_NCBI_Table_Candidates`, and `sample_status.tsv` includes
`codetta_status`.

Published sample and cohort folders retain the logs and stable intermediate
files used to assemble the master table, while large raw CheckM2, BUSCO,
Prokka, CRISPRCasFinder, and Codetta artefacts are pruned before task
completion. BUSCO is also run with its native `--tar` option so file-heavy
intermediate BUSCO subdirectories are compressed before the existing curated
BUSCO summary/log view is retained. Under the 20k OIST override,
`results/cohort/fastani/` keeps only
`ani_metadata.tsv`, `ani_exclusions.tsv`, `fastani_paths.txt`, and the FastANI
matrix/log outputs; the staged `fastani_inputs/` directory remains in `work/`.
By default, only `16S = Yes` samples pass the ANI 16S gate. Supply the naked
flag `--ani_allow_incomplete_16s` to also allow `16S = No` and
`16S = partial`; `16S = NA` remains excluded.

## Notes

- PADLOC and eggNOG outputs are retained per sample but are not merged into the final master table.
- Original accessions remain the published sample-folder names. Internal sanitized IDs are execution-only.
- The shared Python helper image now includes `numpy` and `scipy` so ANI clustering and representative selection reuse the same helper container as the other Python tasks.
- Codetta provenance is split deliberately: `tool_and_db_versions.tsv` reports
  Codetta as `v2.0`, uses the default image
  `quay.io/asuq1617/codetta:2.0`, and records the pinned upstream `main`
  source commit used for that image build.
