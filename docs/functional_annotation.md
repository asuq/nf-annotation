# Functional annotation in v0.4

The shared annotation workflow is implemented in `main.nf` and `reannotate.nf`.
Biological qualification of the complete five-genome cohort remains in progress;
see the [release specification](development/v0.4.md) for the open release gates.

## Configuration

All five tools are enabled by default: `eggnog,cogclassifier,pfam,kofam,padloc`.
Every enabled tool analyses each eligible protein bundle in the declared cohort,
including retained genomes during a cohort update. Low-quality genomes can be
annotated when genetic-code selection and the protein bundle are valid. ANI
eligibility is independent of functional annotation.

Provide a configuration file with immutable runtimes and prepared databases:

```groovy
params {
    annotation_tools = 'eggnog,cogclassifier,pfam,kofam,padloc'
    annotation_cpus = 4
    annotation_memory = '32 GB'
    annotation_max_forks = 2

    python_container = '/path/to/python-scipy.sif'
    eggnog_container = '/path/to/eggnog.sif'
    cogclassifier_container = '/path/to/cogclassifier.sif'
    pfam_container = '/path/to/pfam.sif'
    kofam_container = '/path/to/kofam.sif'
    padloc_container = '/path/to/padloc.sif'

    eggnog_db = '/path/to/prepared/eggnog/7.0'
    cogclassifier_db = '/path/to/prepared/cogclassifier/COG2024'
    pfam_db = '/path/to/prepared/pfam/38.2'
    kofam_db = '/path/to/prepared/kofam/2026-07-02'
}
```

SIF references must be absolute existing paths; preflight hashes their contents.
Docker references must contain `@sha256:<64 hex characters>` or be a local
`sha256:<64 hex characters>` image ID. Mutable annotation image tags are rejected.
The shared Python helper runtime is part of the interpretation identity and
must also be immutable when annotation is enabled. Each tool has its own image.
PADLOC's database is pinned inside its image and checked from native provenance.

Prepare the other four databases with the
[resource preparation workflow](development/annotation_resources.md).
Preflight checks all requested resources before starting expensive sample work.
Files are checked once per run, then individual searches consume the recorded
resource identity. Changed database content needs a new prepared destination.

Select a unique comma-separated subset through `annotation_tools`. To disable
all functional tools for an upstream-only run, set `annotation_tools = ''` in
the configuration file. Do not pass an empty string on the Nextflow command
line: Nextflow interprets that parameter as a boolean flag. The old eggNOG
accession selector and free-form eggNOG/PADLOC argument parameters are removed.

The CPU and memory settings are bounded by `max_cpus` and `max_memory` and
recorded in the method identity. Native task allocations must match the plan.
DIAMOND's block size and index chunks are explicit and derived from that
allocation, with reserved memory for annotation. They do not use host RAM.
The `local` profile caps task memory at 16 GB by default. Full eggNOG runs need
a larger explicitly provisioned budget, for example `--max_memory '32 GB'`
with `annotation_memory = '32 GB'`, when that memory is available to the runtime.

## Entrypoints and reuse

Use the configuration with the normal input and upstream-resource parameters:

```bash
nextflow run main.nf -profile slurm,singularity -c annotation.config \
  --sample_csv samples.csv --metadata metadata.tsv \
  --taxdump /path/to/taxdump --checkm2_db /path/to/checkm2 \
  --codetta_db /path/to/codetta --busco_db /path/to/busco \
  --outdir results -work-dir work
```

Relative `genome_fasta` entries in the sample manifest resolve against the
directory from which Nextflow is launched. Metadata and the existing QC,
genetic-code, CRISPR and ANI fields remain in the master table.

Reannotation requires a published native v0.4 source and a separate destination:

```bash
nextflow run reannotate.nf -profile slurm,singularity -c annotation.config \
  --annotation_from /path/to/results \
  --outdir reannotated -work-dir reannotation-work
```

This entrypoint requires no QC, ANI or taxonomy databases. It validates and
copies the published upstream information and protein bundles. It does not
depend on an old work directory. `main.nf --update_from /path/to/results` uses
the same functional plan while analysing newly added genomes upstream.

| Change | Functional action |
| --- | --- |
| Matching search and interpretation identities; valid evidence checksums | Reuse normalized results. |
| Matching search identity; changed normalizer or interpretation | Normalize the validated native output again. |
| Changed database, runtime, search settings or protein input | Run the affected tool for every affected genome. |
| Previously disabled tool is enabled | Analyse all available bundles, including retained genomes. |
| Previous result failed | Run that tool again. |

The complete accession/tool grid, actions and fingerprints are published in
`pipeline_info/annotation_plan.tsv` and `annotation_plan.json`. A matrix cannot
combine successful results from different method identities. Source tables,
bundle identities, raw evidence and normalized evidence are validated before
their respective reuse. Linked work-directory artefacts are rejected.

## Evidence and counts

Native files, version output, logs and exit codes are retained under
`samples/<accession>/annotation/<tool>/raw/`. Normalized detailed tables and
`result.json` accompany them. `annotation_results.json` records the portable
source contract, table checksums, bundles and result identities.

| Source | Primary acceptance rule |
| --- | --- |
| eggNOG v3 | Native high and medium confidence for each annotation field. GO terms use the exported native confidence for each namespace before the native merge. No additional GO ancestors are added. |
| COGclassifier | Native first reported RPS-BLAST assignment; preserve the full ordered category definition. |
| Pfam | Native gathering thresholds and PfamScan clan-overlap resolution. Retain repeated domains as evidence; count distinct proteins per family. |
| KOfamScan | Native adaptive-threshold pass marker using the frozen prokaryotic profile subset. Preserve multiple accepted KOs. |
| PADLOC | Native systems and members, checked against canonical gene coordinates and native HMMER evidence. Count systems separately from members. |

`tables/protein_manifest.tsv`, `gene_coordinates.tsv` and
`protein_function_summary.tsv` retain every valid input protein, including
unannotated proteins. `tables/annotation_status.tsv` is authoritative for all
five tools per genome; the master table and sample status are derived from it.
The master table includes analysed, mapped and accepted protein counts and
fractions, source-specific assignment totals, COG category counts and fractional
counts, and PADLOC system counts.

The six matrices under `tables/functional_matrices/` are:

- `eggnog_ko_counts.tsv`
- `eggnog_go_counts.tsv`
- `eggnog_ec_counts.tsv`
- `kofam_ko_counts.tsv`
- `cog_counts.tsv`
- `pfam_gene_counts.tsv`

Each cell counts distinct proteins assigned to the feature. Rows follow the
master table; columns are the sorted union of observed accepted features.
`feature_catalogue.tsv` records definitions when available, source, resource
identity and method identity. KO sources remain separate.

Zero means a successful analysis with no accepted assignment to that feature.
`NA` means unavailable, failed, disabled or invalid evidence. An invalid optional
field preserves other valid fields but makes every dependent genome-level count
and matrix cell `NA`, preventing partial undercounts. If no genome has an
accepted feature, a matrix still contains every accession and has no feature
columns. The status table distinguishes the different reasons.

Independent native analyses finish even when one tool fails. The workflow
publishes available results and diagnostics, then the final acceptance check
returns a nonzero exit status if any requested analysis failed or lacked a valid
upstream bundle. A native search terminated before it publishes output receives
`failed` status with reason `missing_planned_result`; the Nextflow log and trace
retain its task error. Failure of the workflow engine itself can still interrupt
publication and never constitutes a completed run.

Resume repeats resource preflight and validates imported source results. Planning
uses content checksums and stable source bundle paths so unchanged native
searches remain cached. Changes to planning code invalidate the plan.

The prepared BUSCO lineage catalogue now pins `2026-05-22` for
`bacillota_odb12` and `mycoplasmatota_odb12`, with the publisher's MD5 values from
[the BUSCO release index](https://busco-data.ezlab.org/v5/data/file_versions.tsv).
The previous `2025-05-14` download URLs returned 404 during qualification.
Existing explicitly supplied lineage directories retain their recorded version;
QC scores from different lineage releases should not be treated as identical.

For cohort execution, prepare BUSCO lineages once in the shared database root
then set `prepare_busco_datasets = false`. Every per-sample BUSCO invocation
uses `--offline` and the staged shared lineage, so sample jobs do not download
lineages or query the release server. Missing lineages fail before those jobs
start. Keep acquisition on the login node, separate from cohort submission.
