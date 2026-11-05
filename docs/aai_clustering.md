# Genus-scale AAI clustering

## Running the analysis

v0.5.0 runs FastAAI alongside ANI in `full` and `clustering` runs. QC-only runs
perform neither comparison. Reannotation retains existing clustering information
and evidence without recomputing it.

```bash
nextflow run main.nf -profile singularity \
  --pipeline_stage clustering --domain bacteria \
  --sample_csv cohort.csv --metadata metadata.tsv \
  --update_from published-qc --outdir clustering-results
```

The shared options are `publication_mode`, `import_batches`, `import_max_forks`
and `score_profile`. Their old `ani_` names and `pipeline_stage ani` are rejected.
`ani_threshold` and `ani_16s_policy` remain specific to ANI. `aai_threshold` is a
fraction from 0.30 to 0.90, defaulting to **0.65**. This is an operational
**genus-scale grouping**, not a taxonomic genus assignment or a universal boundary.
Existing taxonomy names are retained.

## Membership and representatives

AAI accepts valid existing protein bundles with resolved genetic codes, a passing
selected-code CheckM2 quality flag, and no exclusion under the existing atypical
assembly policy. It does not require 16S evidence or complete scoring metadata.
Bacteria and Archaea are compared and clustered separately.

Complete linkage requires every pair in a cluster to meet the cutoff. Missing
comparisons cannot connect clusters. Genomes with no usable markers, and genomes
with no informative comparisons in a multi-genome comparison set, remain explicitly
unclustered. A lone usable genome is a cohort singleton, not evidence of genus novelty.

The representative uses the existing assembly, CheckM2, BUSCO, N50 and scaffold
weights selected by `score_profile`. Normalisation is over fully scoreable
candidates within each AAI cluster; centrality is omitted. The existing tie
cascade is retained. Missing scoring metrics exclude only representative candidacy.
A cluster with no eligible candidate is retained with an unavailable representative.
ANI scoring retains its centrality term and scientific rules.

## Evidence and outputs

- `cohort/fastaai/profiles/`: portable selected proteins, their original HMM hits
  and identities. Profiles bind the accession/protein bytes, upstream source,
  HMM models, environment and extraction code.
- `cohort/fastaai/<domain>/`: condensed `jaccard.npy`, `shared.npy` and
  `possible.npy` arrays with accession ordering in `evidence.json`, input records,
  native logs and comparison provenance. Each unordered pair appears once.
  `python bin/aai_evidence.py <domain-directory> pairs.tsv` exports readable TSV.
- `cohort/aai_clusters/`: cluster membership, representative and per-genome
  summary tables, scoring traces and method metadata.
- The master table adds `AAI_Cluster_ID`, `AAI_Is_Representative`,
  `AAI_to_Representative` and `AAI_Score`. Sample status records exclusions,
  missing evidence and representative eligibility separately.

FastAAI's native Jaccard calculations are unchanged. A writer-only adaptation
retains full precision. Clustering uses `1 - mean Jaccard` and the native
calibration's inverse at the selected AAI cutoff; the monotonic transformation
preserves complete-linkage ordering. Reports retain `<30%`, `>90%` and `NA`
instead of interpreting FastAAI's numeric matrix placeholders as measured AAI.
The representative's self-similarity is 100 by identity.

## Work and reuse

Protein prediction is never repeated. At most eight preprocessing batches run,
with at most two concurrent batches and native workers bounded by task CPUs.
Profiles are installed atomically as independent files in the current publication.
Later cohorts can copy matching profiles from their declared source publications.

One native database is built and queried per domain in task-local scratch. Its
standard protein-plus-HMM input path uses the already selected markers, so database
construction repeats neither protein prediction nor HMM searching. The upstream
crystal merge path is not used because native qualification exposed a SQLite
binary/text merge error in that path.

Native comparisons are separate from clustering and scoring. Cutoff and score
changes can use normal same-version `-resume` without repeating comparisons.
Changes to cohort metadata require a fresh output directory under the existing
cohort identity rules.
The condensed clustering representation avoids a second dense square matrix.
This retains complete linkage's quadratic storage requirements; it does not turn
large-cohort clustering into a linear operation.

The native runtime is pinned to FastAAI commit
`3b4904a223590f8227e9c791b6b7ed6a09c67757` and bundled marker models. Completed
v0.5 QC publications are required; v0.4 publications and caches are unsupported.
Source publications are never modified. This unreleased branch uses the locally
built `nf-annotation-fastaai:0.5.0-dev` image; see `docker/fastaai/Dockerfile`.
See the [upstream software](https://github.com/cruizperez/FastAAI) and
[FastAAI paper](https://doi.org/10.1093/nar/gkaf348).
