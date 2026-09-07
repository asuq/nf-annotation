# nf-annotation

[![Nextflow](https://img.shields.io/badge/version-26.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

`nf-annotation` is a Nextflow DSL2 pipeline for cohort-scale genome
annotation, QC, taxonomy expansion, and ANI-based reporting.

The implemented workflow validates the input manifest and metadata, stages
genomes to execution-safe internal IDs, runs Barrnap, paired CheckM2
translation-table predictions, BUSCO, Codetta, Prokka, CRISPRCasFinder,
PADLOC, eggNOG, and FastANI-based clustering, then publishes three final
reporting tables:

- `master_table.tsv`
- `sample_status.tsv`
- `tool_and_db_versions.tsv`

This README is the operator-facing landing page for normal setup and runs. For
full operator detail, acceptance workflows, and troubleshooting, see
[`docs/runbook.md`](docs/runbook.md). For the implemented architecture and
reporting contracts, see
[`docs/implemented_pipeline.md`](docs/implemented_pipeline.md). For the
recommended first real server validation path, see
[`docs/small_cohort_server_test.md`](docs/small_cohort_server_test.md).

## Table of contents

- [Project summary](#project-summary)
- [ANI clustering and representatives](#ani-clustering-and-representatives)
- [Entrypoints](#entrypoints)
- [Prerequisites](#prerequisites)
- [Required inputs and databases](#required-inputs-and-databases)
- [Preparing sample_csv](#preparing-sample_csv)
- [Preparing metadata.tsv](#preparing-metadatatsv)
- [Quick-start workflow](#quick-start-workflow)
- [Main run examples](#main-run-examples)
- [Updating a cohort](#updating-a-cohort)
- [Output summary](#output-summary)
- [Profiles](#profiles)
- [Documentation map](#documentation-map)


## Project summary

The main analysis workflow is designed around these steps:

- validate the sample manifest and metadata table
- stage each genome under a sanitized internal ID for tool-safe execution
- expand taxonomy from a pinned user-supplied taxdump
- run per-sample QC with Barrnap, paired CheckM2 runs, and BUSCO
- retain the logs and stable intermediate files used to assemble
  `master_table.tsv` while pruning bulky raw tool artefacts before task
  completion
- run gcode-gated annotation with Prokka, CRISPRCasFinder, PADLOC, and eggNOG
- build ANI inputs, run all-vs-all FastANI, cluster genomes, and select
  representatives
- publish final cohort tables and a combined versions report

PADLOC and eggNOG outputs are retained in per-sample folders but are not
merged into `master_table.tsv`. Codetta is merged into the final reporting
tables through `GC_Content`, `Codetta_Genetic_Code`,
`Codetta_NCBI_Table_Candidates`, and `codetta_status`. In
`tool_and_db_versions.tsv`, Codetta is reported as
compatible with `v2.0`, the default container image is
`quay.io/asuq1617/codetta:2.0`, and the pinned upstream `main` commit used for
the image build is recorded explicitly.

## ANI clustering and representatives

ANI reporting starts by building the FastANI input set from samples that pass
the cohort-level eligibility gates. A sample is eligible only when it is not
low quality, has an assigned gcode of `4` or `11`, has `16S = Yes` by default,
has the primary BUSCO result needed for scoring, and is not atypical. The one
atypical exception is `unverified source organism`, which is still allowed into
ANI clustering. Supplying `--ani_allow_incomplete_16s` relaxes only the 16S
gate so `16S = No` and `16S = partial` may enter ANI; `16S = NA` remains
excluded.

FastANI then runs all-vs-all on the eligible genomes. The matrix is clustered
with complete linkage at `--ani_threshold`, which defaults to `0.95`. Complete
linkage means every pair of genomes inside a cluster must meet the threshold;
the clustering helper also performs a post-check and fails if a completed
cluster contains a below-threshold or missing ANI pair. Cluster IDs are stable:
clusters are named `C000001`, `C000002`, and so on, ordered first by descending
cluster size and then by the smallest accession in the cluster.

One representative is selected for each ANI cluster. Each cluster member is
scored from six components:

- assembly rank
- CheckM2 quality
- BUSCO completeness-minus-missing
- N50
- scaffold count
- ANI centrality within the cluster

The BUSCO component uses the primary BUSCO lineage, which is the first lineage
in `--busco_lineages`. The CheckM2 component uses the completeness and
contamination values associated with the sample's assigned gcode. Raw BUSCO,
CheckM2, N50, and scaffold values are winsorised for larger clusters,
normalised to the `[0, 1]` range within the cluster, and combined with the
chosen `--ani_score_profile`. ANI centrality is based on the average ANI from
one member to the other members in the same cluster; singleton clusters receive
a centrality value of `1.0`. The scaffold component is inverted before
weighting, so fewer scaffolds contribute a higher component value.

The scoring profiles use these weights in the component order above:

| `--ani_score_profile` | Weights |
| --- | --- |
| `default` | `3.0, 2.0, 1.5, 0.75, 0.5, 0.5` |
| `isolate` | `3.5, 1.5, 1.5, 1.0, 0.25, 0.5` |
| `mag` | `1.0, 3.0, 2.0, 0.5, 1.5, 0.5` |

If the top representative score is tied within the script tolerance, the
pipeline resolves the tie deterministically in this order:

- higher assembly rank
- higher BUSCO complete percentage
- lower BUSCO missing percentage
- lower CheckM2 contamination
- higher CheckM2 completeness
- fewer scaffolds
- higher N50
- lexicographically smallest accession

The published ANI files are under `results/cohort/ani_clusters/`.
`cluster.tsv` records accession-to-cluster membership, `ani_representatives.tsv`
records one representative row per cluster, and `ani_summary.tsv` records the
per-accession ANI fields later merged into `master_table.tsv`: `Cluster_ID`,
`Is_Representative`, `ANI_to_Representative`, and `Score`.

## Entrypoints

The repository exposes two operator-facing Nextflow entrypoints.

### `main.nf`

Use `main.nf` for normal analysis runs. It expects the analysis inputs and
runtime databases to be ready before execution, except that BUSCO lineage
datasets can be downloaded through the workflow if
`--prepare_busco_datasets true` is used.

### `prepare_databases.nf`

Use `prepare_databases.nf` before analysis when one or more runtime databases
are missing on the target machine. It prepares taxdump, CheckM2, Codetta,
BUSCO, and eggNOG resources in place and writes
`runtime_database_report.tsv` plus `nextflow_args.txt` under its `--outdir`.

## Prerequisites

- Nextflow `>= 26.04.0`
- either Docker or Singularity for containerised runs
- access to the configured tool images
- a pinned taxdump directory
- a CheckM2 database directory
- a Codetta profile directory containing `Pfam-A_enone.hmm` and its pressed
  HMMER sidecars
- BUSCO lineage datasets under `--busco_db`, or permission to let the workflow
  download them with `--prepare_busco_datasets true`
- an eggNOG data directory for real eggNOG runs

Typical profile combinations are:

- local container run: `-profile local,docker`
- SLURM container run: `-profile slurm,singularity`
- OIST HPC run: `-profile oist`

## Required inputs and databases

The main workflow expects these inputs:

| Input | Required details |
| --- | --- |
| Sample manifest | CSV with the locked lower-case columns `accession`, `is_new`, `assembly_level`, and `genome_fasta` |
| Metadata table | Delimited table with `Accession` or `accession` as the join key |
| Taxdump | Directory containing at least `names.dmp` and `nodes.dmp` |
| CheckM2 database | Directory containing exactly one top-level `.dmnd` file |
| Codetta database | Directory containing `Pfam-A_enone.hmm` plus `.h3f`, `.h3i`, `.h3m`, and `.h3p` |
| BUSCO source | `--busco_db` with lineage directories, or `--prepare_busco_datasets true` |
| eggNOG database | Directory passed through `--eggnog_db` |

Useful defaults from the implementation:

- `--busco_lineages` defaults to `bacillota_odb12,mycoplasmatota_odb12`
- `--gcode_rule` defaults to `strict_delta`
- `--ani_threshold` defaults to `0.95`
- `--ani_allow_incomplete_16s` defaults to off; when supplied, `16S = No`
  and `16S = partial` samples may enter ANI if all other ANI gates pass
- `--outdir` defaults to `results`

## Preparing sample_csv

Use `sample_csv` as a simple CSV manifest with one genome per row.

```csv
accession,is_new,assembly_level,genome_fasta
GCA_000027325.1,false,NA,/data/genomes/GCA_000027325.1.fasta
MYCO_NEW_001,true,Complete Genome,/data/genomes/MYCO_NEW_001.fasta
```

- headers must stay lower-case and in this exact structure
- `assembly_level` is required when `is_new=true`
- `assembly_level` can be `NA` when `is_new=false`
- `genome_fasta` must point to the input FASTA for that accession

## Preparing metadata.tsv

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

## Quick-start workflow

### 1. Optionally prepare runtime databases

Use the dedicated prep entrypoint when one or more databases do not already
exist in their final tool-consumable locations.

```bash
nextflow run prepare_databases.nf -profile local,docker \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --download_missing_databases true \
  --outdir runtime-prep
```

### 2. Run the main analysis workflow

```bash
nextflow run . -profile local,docker \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --outdir results
```

### 3. Optionally run the stub smoke test

```bash
nextflow run . -profile test -stub-run
```

This validates the DSL wiring, channel contracts, publish locations, and final
report outputs without requiring real external databases or tool containers.

## Main run examples

Local Docker run:

```bash
nextflow run . -profile local,docker \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --outdir results
```

SLURM plus Singularity run:

```bash
nextflow run . -profile slurm,singularity \
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

Use `--slurm_qos 2h` for QoS. For generic SLURM options that start with a
dash, pass the Nextflow parameter with equals syntax, for example
`--slurm_cluster_options='--qos=2h'`.

OIST profile run:

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

MPCDF Viper CPU run, launched from `viper05i`:

```bash
export NXF_APPTAINER_CACHEDIR="/ptmp/$USER/apptainer-cache"

nextflow run . -profile viper-cpu \
  -w "/ptmp/$USER/nf-annotation-work" \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --outdir "/ptmp/$USER/nf-annotation-results"
```

The profile sends compute tasks to Slurm without fixing a partition, account,
or QoS. BUSCO-lineage, taxdump, Codetta, CheckM2, BUSCO-database, and eggNOG
preparation tasks remain intact and run locally because they need internet
access. They use one CPU and at most 16 GB each, with no more than two running
at once. Nextflow keeps at most 250 Slurm tasks outstanding by default, leaving
headroom below Viper's default 300-job per-user submission limit. Lower the
ceiling with `--viper_slurm_queue_size` when other jobs or workflow launches
share that limit. Override the container cache with `--apptainer_cache_dir` and
the default `apptainer/1.4.3` module with
`--viper_apptainer_module <module/name>` when necessary.

Both the `-w` directory and Apptainer cache must be shared `/ptmp` paths.
Viper `/ptmp` is not backed up and inactive files are subject to retention
cleanup, so preserve final results elsewhere when required. `/r` is accessible
only from login nodes; `/tmp` and generic `$TMPDIR` are unsuitable for
Nextflow work or shared container caches. The opt-in `process_local_scratch`
label uses `$JOB_TMPDIR` only for processes explicitly assigned that label.

Large OIST cohort run with the opt-in storage overrides:

```bash
nextflow run . -profile oist \
  -c conf/oist_20k_storage.config \
  -work-dir /flash/path/to/work_nf_annotation_20k \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --busco_db /path/to/busco \
  --eggnog_db /path/to/eggnog-db \
  --singularity_cache_dir /path/to/singularity-cache \
  --outdir /nfs/path/to/results
```

For OIST storage, placing `-work-dir` on Lustre, such as `/flash`, while
placing `--outdir` on NFS can reduce retained Lustre usage because final
published outputs land on NFS. Active Nextflow work still uses the Lustre work
directory during execution.

To resume that run safely, reuse the same launch directory, Lustre-backed
`-work-dir`, and NFS-backed `--outdir`, then append `-resume`. Do not enable
`cleanup = true`, do not delete `.nextflow/cache`, and do not purge the shared
`/flash` work tree between runs.

Main-workflow BUSCO dataset download variant:

```bash
nextflow run . -profile local,docker \
  --sample_csv samples.csv \
  --metadata metadata.tsv \
  --taxdump /path/to/pinned-taxdump \
  --checkm2_db /path/to/checkm2-db \
  --codetta_db /path/to/codetta-db \
  --prepare_busco_datasets true \
  --eggnog_db /path/to/eggnog-db \
  --outdir results
```

## Updating a cohort

Add `--update_from` to reuse published per-sample results when adding or
removing samples. Supply the **complete revised sample CSV**, current metadata,
and a fresh results directory:

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

Samples present in the source manifest reuse their published QC and annotation
results. Only accessions absent from that manifest receive new QC and
annotation; the CSV's `is_new` column keeps its existing metadata meaning.
Omitted samples are excluded from the revised cohort. ANI is recalculated
all-vs-all for the eligible revised cohort, and the workflow rebuilds clusters,
representatives, cohort 16S outputs, taxonomy and final reports.

The original work directory and Nextflow cache are unnecessary. Retained sample
files are copied into the new output, which can itself be used as the source
of a later update. The source results remain unchanged. Allow storage for a
complete copy of the retained sample artefacts.

For removal-only updates, the per-sample database options can be omitted;
`--taxdump` is still required for refreshed taxonomy. Retained FASTAs must be
available through the revised CSV and match the published staged genomes.
Paths and FASTA wrapping may change, but sequence content, case, record IDs
and record order must match. If original inputs are no longer available,
the CSV may point to their published `samples/<accession>/staged/*.fasta` files.

Retained analyses keep their original results and provenance; per-sample tool
settings apply only to additions. ANI options apply to the entire revised
cohort. Cluster IDs and representatives may change with cohort membership.
See [the update runbook](docs/runbook.md#updating-cohort-membership) for
validation requirements, audit files and resume behaviour.

## Output summary

By default the workflow publishes under `results/`:

```text
results/
  cohort/
  pipeline_info/
  samples/
  tables/
```

The most important outputs are:

- `results/tables/master_table.tsv`
- `results/tables/sample_status.tsv`
- `results/tables/tool_and_db_versions.tsv`
- `results/samples/<accession>/...` for the logs and stable intermediate files
  used to build `master_table.tsv`, with bulky raw artefacts pruned before task
  completion
- `results/cohort/16s/` for the intact-only eligible-sample 16S FASTA, the
  partial-only 16S FASTA, and their manifests
- `results/cohort/fastani/` for ANI metadata, exclusions, path lists, and
  matrix outputs; the staged `fastani_inputs/` directory stays in `work/`
- `results/cohort/ani_clusters/` for clustering and representative tables
- `results/cohort/taxonomy/` for expanded taxonomy output
- `results/pipeline_info/` for Nextflow trace, report, timeline, and DAG files

For the full published layout and the exact reporting column contracts, see
[`docs/implemented_pipeline.md`](docs/implemented_pipeline.md).

Run the dependency-locked development test suite with `pixi run test`.

## Profiles

- `local`: local executor
- `docker`: enables Docker execution
- `slurm`: enables the SLURM executor with optional `--slurm_queue`,
  `--slurm_qos`, and `--slurm_cluster_options`
- `singularity`: enables Singularity execution
- `oist`: OIST HPC profile with SLURM and Singularity enabled
- `viper-cpu`: MPCDF Viper CPU profile; launch from `viper05i` with a
  `/ptmp` work directory and Apptainer cache
- `test`: fixture-backed local profile for `-stub-run`
- `debug`: composable profile that sets
  `--eggnog_only_accessions GCA_000027325.1`

Use `debug` only when that eggNOG short-circuit is acceptable. For full eggNOG
validation on HPC, use the normal execution profile instead of `debug`.

## Documentation map

- `README.md`: onboarding, prerequisites, quick-start, and normal run examples
- [`docs/runbook.md`](docs/runbook.md): extended operator detail, runtime
  database preparation notes, acceptance workflows, and troubleshooting
- [`docs/small_cohort_server_test.md`](docs/small_cohort_server_test.md):
  recommended first real server validation path for the tracked small
  acceptance cohort
- [`docs/implemented_pipeline.md`](docs/implemented_pipeline.md): implemented
  workflow structure, published layout, and final reporting contracts
- [`docs/run_pipeline_test.md`](docs/run_pipeline_test.md): wrapper-specific
  usage for `bin/run_pipeline_test.sh`
