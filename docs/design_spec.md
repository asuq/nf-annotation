# Codex procedure manual: Nextflow DSL2 genome-annotation pipeline and master-table builder

This is the historical v1 design. The [implemented pipeline](implemented_pipeline.md)
and [v0.4 specification](development/v0.4.md) describe the current genetic-code,
functional-annotation and published-source contracts where they supersede it.

**Status:** v1 specification frozen on 2026-03-07 from the legacy notebook behaviour plus user clarifications.  
**Blocking ambiguities:** none.  
**A few non-blocking operational defaults were fixed here explicitly so Codex can implement without further questioning.**

---

## 1. Purpose

Refactor the legacy Jupyter-notebook workflow and manual per-sample curation into a single **Nextflow DSL2** pipeline that:

1. runs the required annotation and QC tools on **all input genomes**;
2. uses `is_new` for metadata handling and final reporting;
3. produces per-sample annotation outputs;
4. produces cohort-level outputs (taxonomy expansion, ANI clustering, combined 16S FASTA, status table, versions table);
5. produces one final **master table** with preserved input metadata columns followed by derived columns.

This is a **behaviour-preserving v1 refactor**, not a methodological redesign.

---

## 2. Hard requirements locked by the user

### 2.1 Input policy

- Run all tools on **all** genomes in the sample manifest.
- `is_new` controls how metadata are sourced and merged and must also be emitted in the final master table.
- Required sample-manifest columns are exactly:
  - `accession`
  - `is_new`
  - `assembly_level`
  - `genome_fasta`
- These header names are **lower-case**.
- `assembly_level` is mandatory **only** when `is_new = true`.
- The metadata table must preserve its original column order and content in the master table.
- For `is_new = true` genomes missing from metadata, fill unavailable metadata with `NA`.
- If optional user-supplied metadata for new genomes are present and parseable, use them; if they are invalid, fill `NA` and record a warning in `sample_status.tsv`.

### 2.2 Naming policy

- Keep two identifiers internally:
  - **Original accession**: exact user-provided sample ID.
  - **Internal ID**: sanitised ID used only for tool-safe internal filenames and prefixes.
- Do **not** add the internal ID to the master table.
- Published output folder names must use the **original accession**.
- Internal-ID sanitisation rule:
  1. keep only `A-Z`, `a-z`, `0-9`, `_`;
  2. replace every other character with `_`;
  3. collapse repeated underscores;
  4. trim leading/trailing underscores;
  5. if collisions remain, append a short deterministic suffix.

### 2.3 Taxonomy policy

- The main analysis workflow must **not** fetch taxdump automatically.
- A separate runtime-database prep entry point may download curated taxdump
  releases before the main workflow runs.
- The main analysis workflow still consumes a resolved taxdump path.
- Do **not** infer taxonomy for newly assembled genomes lacking valid `Tax_ID`.

### 2.4 Barrnap / 16S policy

- Run Barrnap on all genomes.
- Only check **rRNA**.
- Output per-sample `best_16S.fna`.
- Produce cohort-level `all_best_16S.fna`.
- Produce cohort-level `all_partial_16S.fna`.
- Move samples with `Low_quality = true` out of the standard cohort 16S
  outputs and into cohort-level `low_quality_best_16S.fna` or
  `low_quality_partial_16S.fna`.
- If a genome is atypical, still run Barrnap and keep the per-sample result, but include its intact 16S in `all_best_16S.fna` only when it is atypical only because of `unverified source organism`.
- Include all non-low-quality partial-only samples in `all_partial_16S.fna`,
  including atypical samples.
- Include all low-quality intact and partial 16S sequences in the matching
  low-quality cohort FASTA, including atypical samples.
- `16S` master-table vocabulary is fixed as:
  - `Yes`
  - `partial`
  - `No`
- If multiple intact 16S hits tie for best Barrnap score:
  1. choose the longest;
  2. if still tied, choose the first hit in GFF order.

### 2.5 CheckM2 / gcode policy

- End user provides the CheckM2 database path.
- Run CheckM2 twice for every sample:
  - forced `--ttable 4`
  - forced `--ttable 11`
- Select code 4 when `Average_Gene_Length_gcode4 / Average_Gene_Length_gcode11 > 1.5`.
- Select code 11 when the valid ratio is at or below 1.5, including exactly 1.5.
- Require valid reports for translation tables 4 and 11, matching genome identities
  and shared statistics, and finite, strictly positive mean lengths.
- Record both means, the ratio, rule, threshold and selection reason. Compare
  decimal values without rounding before selection. Calculate QC from the
  selected code's completeness and contamination.
- `gcode_rule` defaults to `mean_gene_length_ratio`. Retain `strict_delta` and
  `delta_then_11`: select the table whose completeness advantage is strictly
  greater than 10 percentage points. At an absolute difference of 10 or less,
  `strict_delta` leaves the code unresolved and `delta_then_11` selects table 11.
  Invalid report pairs remain unresolved under every rule.
- If gcode is `NA`:
  - record a warning in `sample_status.tsv`;
  - do **not** run Prokka, eggNOG, CRISPRCasFinder, or PADLOC for that sample.
- CheckM2 failure for a sample must **not** kill the whole pipeline; warn in `sample_status.tsv` and continue.

### 2.6 QC policy

- Low-quality rule is fixed as the GTDB-style score:

```text
Low_quality = (chosen_completeness - 5 * chosen_contamination) <= 50
```

- If gcode is `NA`, `Low_quality = NA`.

### 2.7 BUSCO policy

- Use **genome-mode** BUSCO.
- Keep **all** BUSCO outputs.
- Split BUSCO into two separate logical modules:
  1. optional prep workflow to download lineage datasets;
  2. offline BUSCO run using the pre-downloaded datasets.
- Run both lineages for **every** sample:
  - `bacillota_odb12`
  - `mycoplasmatota_odb12`
- The master table must include **both** BUSCO results.
- The **first configured lineage** is the one used later for ANI representative selection.

### 2.8 Prokka policy

- Run Prokka only when gcode is `4` or `11`.
- Use `--compliant`.
- Use `--rfam`.
- Keep the full Prokka output directory.
- Use the **internal ID** as the Prokka `--prefix` / locustag-safe identifier.

### 2.9 CRISPRCasFinder policy

- Separate the CRISPRCasFinder run and CRISPR summary into two processes.
- If CRISPRCasFinder ran for the sample, the old atypical-strain exclusion for CRISPR summarisation is **fully removed**.
- Use the exact staged genome FASTA as input.
- If a suitable BioContainers image does not exist, create a Dockerfile under `docker/` and let the user build/upload it later.

### 2.10 Codetta policy

- Run Codetta on **all** samples, independent of CheckM2 gcode assignment.
- Use a prepared Codetta database directory rooted at `Pfam-A_enone.hmm`.
- Support Codetta database preparation in the separate runtime-database prep workflow.
- Publish all Codetta-generated artefacts for each sample under `samples/<original_accession>/codetta/`.
- Report Codetta publicly as version `v2.0`, while pinning the custom runtime image build to the latest verified upstream `main` commit and recording that source commit in `tool_and_db_versions.tsv`.
- Keep Codetta soft-failing per sample: if Codetta fails, continue the cohort run, write failure status into `sample_status.tsv`, and emit `NA` in the Codetta master-table fields.
- Add two Codetta-derived columns to the master table immediately after `Gcode`:
  - `Codetta_Genetic_Code`
  - `Codetta_NCBI_Table_Candidates`
- Add `codetta_status` to `sample_status.tsv`.
- Codetta table matching must compare the inferred 64-character code against NCBI translation tables while:
  - ignoring `?` positions
  - treating amino-acid versus stop mismatches as incompatible
  - collecting **all** compatible table IDs in ascending order
- If no NCBI table matches, write `unassigned`.

### 2.11 PADLOC policy

- Do not add PADLOC results to the master table.
- Emit a clear warning at pipeline start that PADLOC outputs are produced but not merged into the master table.
- Pre-process Prokka GFF for PADLOC by:
  - removing `gnl|Prokka|` from feature IDs;
  - removing the FASTA tail from the GFF.

### 2.12 eggNOG policy

- Do not add eggNOG results to the master table.
- Keep all final eggNOG outputs except the temporary directory.
- Pre-process FAA headers by replacing whitespace with underscores only on FASTA header lines.

### 2.13 ANI clustering policy

- Use the existing `cluster_ani.py` logic as the starting point.
- Default ANI threshold is **0.95**.
- Exclude from ANI clustering:
  - low-quality genomes;
  - genomes without **complete** 16S (`16S != Yes`);
  - atypical genomes.
- Exception: genomes atypical **only because of** `unverified source organism` must still be included in ANI clustering.
- Optional CLI exception: `--ani_allow_incomplete_16s` relaxes only the 16S
  gate so `16S = No` and `16S = partial` samples may be included; `16S = NA`
  remains excluded.
- The BUSCO lineage prioritised for ANI representative scoring is the **first lineage in the configured list**.

### 2.14 Final reporting policy

- Produce `sample_status.tsv`.
- Add Codetta summary values to the final reporting layer.
- Add `is_new` to the master table as the first appended column after the
  preserved metadata block.
- Add `GC_Content` to the master table immediately after
  `Total_Coding_Sequences_gcode11` and before `Gcode`.
- Derive `GC_Content` from `seqtk comp` as `(G + C) / (A + C + G + T) * 100`,
  ignoring ambiguous bases and formatting the result as a plain decimal without
  a percent sign.
- Produce a final plain-text versions report.
- Prefer BioContainers images.
- If BioContainers do not cover a tool, create a module-specific Dockerfile in `docker/`.

---

## 3. Non-blocking defaults frozen here so implementation can proceed

These were not fully specified, so they are fixed here for v1.

### 3.1 Gcode-NA behaviour in ANI clustering

`cluster_ani.py` requires a concrete gcode to choose the CheckM2 metrics used in representative scoring. Therefore, in v1:

- genomes with `Gcode = NA` are **excluded** from ANI clustering;
- they receive `NA` in the master-table ANI columns:
  - `Cluster_ID`
  - `Is_Representative`
  - `ANI_to_Representative`
  - `Score`
- `sample_status.tsv` must record the exclusion reason as `gcode_na`.
- `Gcode = NA` is reserved for failed, missing, invalid or inconsistent paired
  CheckM2 reports. A valid ratio at or below 1.5 resolves to code 11.

### 3.2 Atypical-genome detection rule

If the metadata table contains `Atypical_Warnings`, then:

- `is_atypical = true` when `Atypical_Warnings` is non-empty and not `NA`;
- `is_atypical_exception = true` only when the parsed atypical-reason list contains one or more entries and every entry is the case-insensitive phrase `unverified source organism`.

If `Atypical_Warnings` is absent entirely:

- treat `is_atypical = false` for clustering;
- record a pipeline-level warning in the versions/reporting area.

### 3.3 Original accession as output folder name

Published output folders must use the original accession. Therefore validation must **hard-fail** if an accession contains path separators or other filesystem-breaking characters that make a directory impossible to create safely on the target platform.

### 3.4 Optional extra columns in the sample manifest

The required sample-manifest columns are exactly the four listed above, but the pipeline may tolerate extra lower-case columns (for example `organism_name`, `tax_id`, `assembly_name`) and use them as supplemental metadata for `is_new = true` genomes.

### 3.5 Status vocabulary

Use the following status vocabulary in `sample_status.tsv`:

- `done`
- `failed`
- `skipped`
- `na`
- `warn`

The main module status columns should use `done`, `failed`, `skipped`, or `na`; a separate semicolon-delimited warning field can carry `warn` messages.

### 3.6 Versions report format

The final versions report will be a tab-separated text file named:

```text
tool_and_db_versions.tsv
```

Recommended columns:

- `component`
- `kind` (`tool`, `database`, `container`, `pipeline`, `runtime`)
- `version`
- `image_or_path`
- `notes`

### 3.7 Execution profiles

Even though the primary runtime target was not locked explicitly, Codex should implement at least:

- `-profile local,docker`
- `-profile slurm,singularity`

This is an implementation convenience, not a scientific rule.

### 3.8 Retry behaviour

- Nextflow task retries remain controlled by `params.task_attempts`.
- Database download task retries remain controlled by
  `params.db_download_attempts`.
- Sample-level modules that intentionally degrade on tool failure instead of
  killing the whole pipeline may also retry the wrapped tool internally before
  emitting placeholder outputs.
- That internal soft-fail retry count is controlled separately by
  `params.soft_fail_attempts`.

---

## 4. Repository layout to implement

```text
.
├── main.nf
├── nextflow.config
├── conf/
│   ├── base.config
│   ├── local.config
│   ├── slurm.config
│   ├── docker.config
│   └── singularity.config
├── modules/local/
│   ├── validate_inputs.nf
│   ├── stage_inputs.nf
│   ├── taxonomy_expand.nf
│   ├── barrnap.nf
│   ├── summarise_16s.nf
│   ├── checkm2.nf
│   ├── assign_gcode_and_qc.nf
│   ├── download_busco_dataset.nf
│   ├── busco.nf
│   ├── summarise_busco.nf
│   ├── codetta.nf
│   ├── summarise_codetta.nf
│   ├── prokka.nf
│   ├── prepare_padloc_inputs.nf
│   ├── ccfinder.nf
│   ├── summarise_ccfinder.nf
│   ├── padloc.nf
│   ├── eggnog.nf
│   ├── build_fastani_inputs.nf
│   ├── fastani.nf
│   ├── cluster_ani.nf
│   ├── build_master_table.nf
│   ├── write_sample_status.nf
│   └── collect_versions.nf
├── subworkflows/local/
│   ├── per_sample_annotation.nf
│   ├── per_sample_qc.nf
│   ├── cohort_taxonomy.nf
│   ├── cohort_16s.nf
│   └── cohort_ani.nf
├── bin/
│   ├── validate_inputs.py
│   ├── sanitise_accessions.py
│   ├── summarise_16s.py
│   ├── summarise_checkm2.py
│   ├── summarise_busco.py
│   ├── summarise_codetta.py
│   ├── summarise_ccfinder.py
│   ├── build_master_table.py
│   ├── build_sample_status.py
│   └── cluster_ani.py
├── docker/
│   └── <tool>/Dockerfile
├── assets/
│   ├── master_table_append_columns.txt
│   ├── sample_status_columns.txt
│   └── testdata/
└── docs/
    ├── design_spec.md
    ├── runbook.md
    └── change_log.md
```

---

## 5. High-level workflow graph

```text
validate_inputs
  -> stage_inputs
  -> taxonomy_expand

per-sample branch (all samples):
  barrnap
  checkm2_ttable4
  checkm2_ttable11
  assign_gcode_and_qc
  busco_bacillota
  busco_mycoplasmatota
  codetta

per-sample gated branch (only if gcode = 4 or 11):
  prokka
  ccfinder
  summarise_ccfinder
  prepare_padloc_inputs
  padloc
  eggnog

cohort branch:
  build_cohort_16s (standard and low-quality cohort 16S FASTA outputs)
  build_fastani_inputs
  fastani_all_vs_all
  cluster_ani

final aggregation:
  build_master_table
  write_sample_status
  collect_versions
```

Key points:

- BUSCO is **not** gated by gcode.
- Prokka, CRISPRCasFinder, PADLOC, and eggNOG **are** gated by gcode.
- ANI clustering is cohort-level and runs after the per-sample summaries required for inclusion filtering.

---

## 6. Input contracts

## 6.1 Sample manifest

Required file: CSV.

Required columns:

| Column | Type | Rule |
|---|---:|---|
| `accession` | string | unique per row |
| `is_new` | boolean-like | accepted forms should be normalised to true/false |
| `assembly_level` | string | mandatory when `is_new=true`; optional otherwise |
| `genome_fasta` | path | must exist |

Validation rules:

- no duplicate `accession`
- every `genome_fasta` must exist
- `is_new=false` requires a metadata row for the accession
- `is_new=true` may omit metadata row
- if `is_new=true`, `assembly_level` must be present and parseable
- if `accession` sanitises to an empty string, hard-fail
- if `accession` contains path separators, hard-fail

## 6.2 Metadata table

Required behaviour:

- preserve the original metadata columns exactly as given and in the same order;
- join by accession;
- append derived columns after the metadata block;
- hard-fail if metadata already contains `is_new`, because that name is reserved
  for the appended reporting column;
- if the metadata key is named `accession`, normalise it internally to `Accession` for joining;
- otherwise require an `Accession` column.

For `is_new=true` genomes:

- if the accession is missing from the metadata table, create a metadata row with `NA` for unknown fields;
- if optional extra fields are provided elsewhere and parse correctly, use them;
- if optional values are malformed, set `NA` and log a warning.

## 6.3 External resources provided by the end user

Mandatory runtime resources:

- pinned NCBI taxdump directory
- CheckM2 database path
- Codetta profile database path rooted at `Pfam-A_enone.hmm`
- BUSCO datasets or ability to run the BUSCO prep workflow
- eggNOG database path
- PADLOC database/resources required by the chosen PADLOC image
- container images or Dockerfiles for tools not covered by BioContainers

---

## 7. Derived identifiers and naming

## 7.1 Original accession

Used for:

- master table
- sample-status table
- published output directory names
- final combined cohort tables

## 7.2 Internal ID

Used for:

- temporary staged FASTA filenames
- tool-safe prefixes
- Prokka `--prefix`
- Prokka locustag-safe identifier
- temporary intermediate filenames

## 7.3 Collision handling

If two accessions sanitise to the same internal ID, append a deterministic suffix, for example a short stable hash of the original accession.

The mapping between original accession and internal ID must be recorded internally and exposed at least in `sample_status.tsv` or an auxiliary mapping file.

---

## 8. Module-by-module implementation specification

## 8.1 `validate_inputs`

Purpose:

- validate the sample manifest and metadata table
- normalise booleans
- build the accession-to-internal-ID mapping
- decide which metadata rows are missing
- generate a validated manifest for downstream modules

Outputs:

- `validated_samples.tsv`
- `accession_map.tsv` with at least:
  - `accession`
  - `internal_id`
  - `is_new`
  - `assembly_level`
  - `genome_fasta`
- `validation_warnings.tsv`

Hard-fail conditions:

- malformed sample manifest
- duplicate accessions
- missing genome FASTA
- invalid `is_new`
- `is_new=true` without `assembly_level`
- `is_new=false` without metadata row
- filesystem-unsafe accession for published output folders

## 8.2 `stage_inputs`

Purpose:

- stage each genome to a tool-safe canonical FASTA using the internal ID
- avoid weird sample naming breaking tools

Outputs per sample:

- staged FASTA named with the internal ID
- optional `.fai` if indexing is useful downstream

Published output naming:

- use the original accession as the visible folder name
- internal filenames may still be used inside the work directory

## 8.3 `taxonomy_expand`

Purpose:

- expand lineage from the pinned taxdump using metadata `Tax_ID`

Rules:

- consume the prepared taxdump path supplied to the main workflow
- do not infer taxonomy for new genomes missing valid `Tax_ID`
- join lineage back by `Tax_ID`

Required derived columns:

- `superkingdom`
- `phylum`
- `class`
- `order`
- `family`
- `genus`
- `species`

If `Tax_ID` is absent or invalid:

- set taxonomy columns to `NA`
- warn in `sample_status.tsv`

## 8.4 `barrnap`

Purpose:

- detect rRNA and enable 16S status calls

Required outputs per sample:

- `rrna.gff`
- `rrna.fa`
- log files as appropriate

Implementation note:

- preserve the user's legacy Barrnap behaviour of producing both GFF and FASTA.
- verify the exact CLI flags against the containerised Barrnap version selected for the pipeline. The required outcome is unambiguous even if packaging-specific flags differ.

## 8.5 `summarise_16s`

Purpose:

- derive the sample-level `16S` status
- select `best_16S.fna`

Rules:

- inspect only Barrnap rRNA results
- classify:
  - `No` = no 16S hit
  - `partial` = 16S exists but all hits are partial
- `Yes` = at least one intact 16S hit
- choose best intact 16S by minimum Barrnap score
- if no intact hit exists, choose best partial 16S by the same score and tie-break rules
- tie-break: longest, then first in GFF order

Required outputs:

- per-sample `best_16S.fna`
- per-sample `16S_status.tsv`

Implementation note:

- cohort-level `all_best_16S.fna`, `all_best_16S_manifest.tsv`,
  `all_partial_16S.fna`, `all_partial_16S_manifest.tsv`,
  `low_quality_best_16S.fna`, `low_quality_best_16S_manifest.tsv`,
  `low_quality_partial_16S.fna`, and
  `low_quality_partial_16S_manifest.tsv` are built in the cohort 16S branch
  from the per-sample Barrnap summaries and CheckM2 QC summaries
- intact cohort inclusion is applied there from metadata `Atypical_Warnings`
  using the locked `unverified source organism` exception
- partial cohort inclusion is applied there from `16S = partial`
- standard cohort outputs exclude only literal `Low_quality = true`; `false`
  and `NA` use the ordinary intact or partial inclusion rules
- low-quality cohort outputs include literal `Low_quality = true` samples with
  matching `16S = Yes` or `16S = partial`, without atypical filtering

## 8.6 `checkm2`

Run as two separate invocations per sample:

- `checkm2_gcode4`
- `checkm2_gcode11`

Required outputs per run:

- full CheckM2 output directory
- `quality_report.tsv`
- tool log

The user-specified legacy command pattern must be preserved semantically:

```bash
${CHECKM2} predict --extension fasta --input "${genome}" --output-directory "checkm2_gcode11" \
    --database_path "${CHECKM2_DB}" --threads $cpus \
    --ttable 11 --force
```

The `ttable 4` run is analogous.

## 8.7 `assign_gcode_and_qc`

Purpose:

- merge the two CheckM2 reports per sample
- verify shared assembly stats are consistent
- assign gcode
- compute `Low_quality`

Required derived columns:

- `Completeness_gcode4`
- `Completeness_gcode11`
- `Contamination_gcode4`
- `Contamination_gcode11`
- `Coding_Density_gcode4`
- `Coding_Density_gcode11`
- `Average_Gene_Length_gcode4`
- `Average_Gene_Length_gcode11`
- `Total_Coding_Sequences_gcode4`
- `Total_Coding_Sequences_gcode11`
- `Gcode`
- `Low_quality`

Rules:

- if one of the two CheckM2 runs is missing or malformed, record `failed` for CheckM2 in `sample_status.tsv` and continue
- if gcode is `NA`, record `warn` and gate downstream annotation modules off for that sample
- `Low_quality = NA` when gcode is `NA`

## 8.8 `download_busco_dataset`

Purpose:

- optional prep workflow for offline BUSCO

Rules:

- not part of the normal mandatory analysis path
- can be invoked when datasets are not yet staged
- datasets must be reused offline later

Outputs:

- downloaded BUSCO lineage directories
- download log

## 8.9 `busco`

Run separately for each configured lineage and sample.

Required behaviour:

- genome mode
- offline mode
- retain full BUSCO output directory

Configured default lineages:

1. `bacillota_odb12`
2. `mycoplasmatota_odb12`

## 8.10 `summarise_busco`

Purpose:

- parse BUSCO machine-readable summary output
- avoid brittle line-number scraping from text summaries

Derived master-table columns:

- `BUSCO_bacillota_odb12`
- `BUSCO_mycoplasmatota_odb12`

Rules:

- one column per configured lineage, named `BUSCO_<lineage>`
- the first configured lineage is the one later used by ANI representative selection
- if a BUSCO run fails, set the relevant column to `NA` and warn in `sample_status.tsv`

## 8.11 `prokka`

Runs only when `Gcode` is `4` or `11`.

Required command semantics:

```bash
prokka "${genome}" -o "${RUN_DIR}/prokka" --prefix "${SAMPLE}" --locustag LOCUS \
    --compliant --gcode "${GCODE}" --cpus 10 --rfam
```

Implementation refinement for v1:

- keep `--compliant`
- keep `--rfam`
- use the internal ID for `--prefix`
- use an internal-ID-safe locustag rather than a globally reused literal `LOCUS`

Required outputs:

- full Prokka output directory
- especially retain `.faa`, `.ffn`, `.gff`, `.gbk`, `.tsv`

## 8.12 `ccfinder`

Runs only when `Gcode` is `4` or `11`.

Required semantics:

- use the staged genome FASTA
- pass the assigned gcode
- preserve the user-specified CRISPRCasFinder invocation pattern
- allow container fallback via custom Dockerfile if BioContainers coverage is missing

Required outputs:

- full `ccfinder` output directory
- logs
- `result.json`

## 8.13 `summarise_ccfinder`

Purpose:

- parse `result.json`
- produce strain-, contig-, and CRISPR-level summary tables

Rules preserved from notebooks:

- ignore arrays with `Evidence_Level == 1`
- compute strain-level summary fields:
  - `CRISPRS`
  - `SPACERS_SUM`
  - `CRISPR_FRAC`

Required outputs:

- `ccfinder_strains.tsv`
- `ccfinder_contigs.tsv`
- `ccfinder_crisprs.tsv`

Master-table contribution:

- only the **strain-level** summary contributes to the master table in v1

## 8.14 `prepare_padloc_inputs`

Purpose:

- convert Prokka outputs to PADLOC-safe inputs

Required transformations:

1. remove `gnl|Prokka|` from GFF content
2. remove the FASTA tail from the GFF

Required outputs:

- cleaned GFF
- original or copied FAA path

## 8.15 `padloc`

Runs only when `Gcode` is `4` or `11`.

Rules:

- use cleaned GFF plus Prokka FAA
- keep all PADLOC outputs
- do not merge PADLOC results into the master table
- pipeline must emit a startup warning stating this explicitly

## 8.16 `eggnog`

Runs only when `Gcode` is `4` or `11`.

Required pre-processing:

```bash
awk '
  /^>/ { gsub(/[[:space:]]+/, "_"); }
  { print }
' "${proteome}" > "${clean_faa}"
```

Rules:

- retain all eggNOG outputs except the temp directory
- create a cleaned tabular annotation file with the comment lines removed, as in the user-provided legacy snippet
- do not merge eggNOG into the master table in v1

## 8.17 `build_fastani_inputs`

Purpose:

- decide which genomes are eligible for ANI clustering
- build the canonical list file(s) for FastANI
- build the metadata subset required by `cluster_ani.py`

Inclusion rule for ANI:

A genome is eligible only if all of the following are true:

- `Low_quality == false`
- `16S == Yes`
- `Gcode` is `4` or `11`
- not atypical, **or** atypical only because of `unverified source organism`

When `params.ani_allow_incomplete_16s` is true, `16S = No` and
`16S = partial` pass the 16S gate. `16S = NA` remains excluded.

For excluded genomes:

- record the exclusion reason(s) in `sample_status.tsv`
- set ANI-derived master-table fields to `NA`

## 8.18 `fastani`

Purpose:

- perform all-vs-all ANI for the eligible genome subset

Required behaviour:

- run FastANI with matrix output compatible with the existing `cluster_ani.py`
- retain the FastANI matrix and raw output files

Default threshold later used for clustering:

- `0.95`

## 8.19 `cluster_ani`

Purpose:

- reuse and adapt the existing `cluster_ani.py`
- cluster genomes
- choose representatives

Required adaptation for v1:

- make the primary BUSCO column configurable based on the first lineage in `params.busco_lineages`
- accept only genomes with valid `Gcode` (`4` or `11`)
- keep representative-scoring logic otherwise unchanged unless forced by input-schema normalisation

Master-table contribution:

- `Cluster_ID`
- `Is_Representative`
- `ANI_to_Representative`
- `Score`

## 8.20 `build_master_table`

Purpose:

- produce the final one-row-per-sample master table

Join strategy:

- start from the validated sample manifest to guarantee every requested sample appears once
- join metadata by accession
- append taxonomy by `Tax_ID`
- append per-sample derived tables by accession
- append ANI results by accession
- reorder columns to the fixed output schema

Hard-fail conditions:

- duplicate final accession rows
- missing required append columns
- inconsistent accession mapping
- any join producing ambiguous duplicates

## 8.21 `write_sample_status`

Purpose:

- create one authoritative audit table of sample-level execution status, warnings, skips, and ANI-inclusion decisions

Recommended columns:

- `accession`
- `internal_id`
- `is_new`
- `validation_status`
- `taxonomy_status`
- `barrnap_status`
- `checkm2_gcode4_status`
- `checkm2_gcode11_status`
- `gcode_status`
- `gcode`
- `low_quality`
- `busco_bacillota_odb12_status`
- `busco_mycoplasmatota_odb12_status`
- `prokka_status`
- `ccfinder_status`
- `padloc_status`
- `eggnog_status`
- `ani_included`
- `ani_exclusion_reason`
- `warnings`
- `notes`

Rules:

- warnings should accumulate rather than overwrite each other
- `ani_exclusion_reason` may be a semicolon-delimited list, e.g. `low_quality;partial_16s;atypical`

## 8.22 `collect_versions`

Purpose:

- gather tool versions, database versions/paths, container references, and runtime information into one final text report

Required output:

- `tool_and_db_versions.tsv`

Recommended rows include:

- Nextflow version
- pipeline version / git commit if available
- each tool version
- each container image tag or URI
- CheckM2 DB path/version
- taxdump path/date label
- BUSCO lineage names and versions
- eggNOG DB path/version label
- PADLOC DB version label if recoverable

---

## 9. Master table specification

## 9.1 General rule

The master table must contain **exactly one row per sample manifest accession**.

The table is built as:

1. preserved metadata columns in original order;
2. appended derived columns in fixed order.

## 9.2 Preserved metadata block

- Preserve all metadata columns exactly as supplied.
- Do not reorder or rename them except for the internal join-key normalisation needed to match the sample manifest.
- Do not inject `is_new` into the metadata block; emit it only in the appended reporting block.
- For `is_new=true` genomes absent from metadata, generate the metadata block with `NA` values unless user-supplied values are available.

## 9.3 Appended derived columns in fixed order

Append after the metadata block:

### Manifest / reporting block

- `is_new`

### Taxonomy block

- `superkingdom`
- `phylum`
- `class`
- `order`
- `family`
- `genus`
- `species`

### CheckM2 block

- `Completeness_gcode4`
- `Completeness_gcode11`
- `Contamination_gcode4`
- `Contamination_gcode11`
- `Coding_Density_gcode4`
- `Coding_Density_gcode11`
- `Average_Gene_Length_gcode4`
- `Average_Gene_Length_gcode11`
- `Total_Coding_Sequences_gcode4`
- `Total_Coding_Sequences_gcode11`
- `GC_Content`

### Gcode / QC block

- `Gcode`
- `Codetta_Genetic_Code`
- `Codetta_NCBI_Table_Candidates`
- `Low_quality`
- `16S`

### BUSCO block

- `BUSCO_bacillota_odb12`
- `BUSCO_mycoplasmatota_odb12`

If more BUSCO lineages are ever configured, append `BUSCO_<lineage>` columns in the same order as the configured lineage list.

### CRISPR block

- `CRISPRS`
- `SPACERS_SUM`
- `CRISPR_FRAC`

### ANI block

- `Cluster_ID`
- `Is_Representative`
- `ANI_to_Representative`
- `Score`

## 9.4 Value conventions

- `Gcode` values: `4`, `11`, `NA`
- `Codetta_Genetic_Code`: 64-character Codetta code string or `NA`
- `Codetta_NCBI_Table_Candidates`: one table ID, multiple IDs joined by `;`, `unassigned`, or `NA`
- `Low_quality` values: `true`, `false`, `NA`
- `16S` values: `Yes`, `partial`, `No`
- `Is_Representative` values: `yes`, `no`, `NA`
- missing values elsewhere: `NA`

---

## 10. ANI representative-selection input policy

The ANI representative-selection metadata table must be built only from ANI-eligible genomes.

Required fields available to `cluster_ani.py`:

- accession
- path
- assembly level
- chosen BUSCO primary lineage column (first lineage in the configured list)
- scaffolds
- genome size
- N50
- chosen CheckM2 completeness/contamination according to gcode
- organism name if available

For `is_new=true` genomes:

- `assembly_level` comes from the sample manifest
- `organism_name` may come from optional supplemental metadata if present, else `NA`
- if `Assembly_Name` is needed but absent, default to the original accession

---

## 11. Gating and failure policy

## 11.1 Pipeline-level hard-fail

Hard-fail the whole pipeline on:

- invalid global inputs
- invalid sample manifest schema
- duplicate accession collisions that cannot be resolved deterministically
- metadata join ambiguity
- final master-table duplicate rows
- FastANI cohort step failure
- ANI clustering failure for the eligible subset
- final table writer failure

## 11.2 Sample-level non-fatal failures

Continue the pipeline but record failure in `sample_status.tsv` for per-sample problems in:

- taxonomy expansion
- Barrnap
- CheckM2
- BUSCO
- Codetta
- Prokka
- CRISPRCasFinder
- PADLOC
- eggNOG

Downstream gating must behave sensibly after failure, for example:

- if CheckM2 fails, gcode becomes `NA`, downstream gcode-dependent modules are skipped
- if Barrnap fails, `16S = NA` or `No` according to the summariser policy and the sample is excluded from ANI
- if BUSCO fails, ANI representative selection may exclude the sample if the primary BUSCO metric required by `cluster_ani.py` is unavailable

For v1, the recommended ANI rule is conservative: missing required ANI-scoring metadata should exclude the genome and record the reason.

---

## 12. Output tree

Recommended published structure:

```text
results/
├── samples/
│   └── <original_accession>/
│       ├── barrnap/
│       ├── checkm2_gcode4/
│       ├── checkm2_gcode11/
│       ├── busco/
│       │   ├── bacillota_odb12/
│       │   └── mycoplasmatota_odb12/
│       ├── codetta/
│       ├── prokka/
│       ├── ccfinder/
│       ├── padloc/
│       └── eggnog/
├── cohort/
│   ├── taxonomy/
│   ├── 16s/
│   │   ├── all_best_16S.fna
│   │   ├── all_best_16S_manifest.tsv
│   │   ├── all_partial_16S.fna
│   │   ├── all_partial_16S_manifest.tsv
│   │   ├── low_quality_best_16S.fna
│   │   ├── low_quality_best_16S_manifest.tsv
│   │   ├── low_quality_partial_16S.fna
│   │   └── low_quality_partial_16S_manifest.tsv
│   ├── fastani/
│   └── ani_clusters/
└── tables/
    ├── master_table.csv
    ├── sample_status.tsv
    └── tool_and_db_versions.tsv
```

Published per-sample tool folders are curated artefact views, not mirrored raw
work directories. For `checkm2_gcode4/`, `checkm2_gcode11/`, `busco/<lineage>/`,
`prokka/`, and `ccfinder/`, keep the published logs and summary artefacts at
the tool root and do not republish nested raw subdirectories when they would
only duplicate those same files.

---

## 13. Nextflow implementation rules for Codex

- Use **DSL2**.
- Put orchestration in Nextflow and business logic in small Python CLIs under `bin/`.
- Avoid writing complex merge logic in Groovy.
- Every module must have explicit inputs and named outputs.
- Prefer BioContainers images pinned by tag or digest.
- If a tool lacks a suitable image, create `docker/<tool>/Dockerfile`.
- Implement `stub` sections for testability.
- Capture module versions in a machine-readable per-module file and merge them at the end.
- Expose database paths and lineage order as params, not hardcoded constants.
- Preserve the order of `params.busco_lineages`; the first lineage controls ANI representative selection.

Recommended params include:

- `params.sample_csv`
- `params.metadata`
- `params.taxdump`
- `params.checkm2_db`
- `params.codetta_db`
- `params.busco_lineages`
- `params.busco_db`
- `params.eggnog_db`
- `params.outdir`
- `params.ani_threshold`

---

## 14. Minimal acceptance tests before declaring v1 finished

Construct a frozen test set that covers at least:

1. one genome assigned gcode 4
2. one genome assigned gcode 11
3. one genome with ambiguous gcode -> `NA`
4. one genome with no 16S
5. one genome with only partial 16S
6. one atypical genome excluded from ANI
7. one atypical genome with `unverified source organism` included in ANI
8. one CRISPR-positive genome
9. one CRISPR-negative genome
10. at least one pair of genomes that cluster together at ANI 0.95
11. one `is_new=true` genome missing from metadata
12. one accessions-sanitisation collision case

Acceptance criteria:

- all requested output files exist
- `master_table.csv` has one row per sample accession
- metadata block is preserved exactly
- append columns are present in the locked order
- samples with `Gcode = NA` skip Prokka/CRISPR/PADLOC/eggNOG
- Codetta runs even when `Gcode = NA`
- ANI inclusion/exclusion reasons match the rules
- `tool_and_db_versions.tsv` is populated
- `sample_status.tsv` explains all warnings and skips

---

## 15. Explicit out-of-scope items for v1

Unless later added deliberately, v1 does **not** include:

- automatic taxdump refresh
- taxonomy inference for new genomes without valid `Tax_ID`
- PADLOC-to-master-table integration
- eggNOG-to-master-table integration
- manual 16S phylogroup curation/integration from notebooks
- any extra manual curation layers not encoded in the current automated rules

---

## 16. Recommended implementation order for Codex

1. Write `docs/design_spec.md` from this manual.
2. Implement and test the standalone Python CLIs first:
   - `validate_inputs.py`
   - `sanitise_accessions.py`
   - `summarise_16s.py`
   - `summarise_checkm2.py`
   - `summarise_busco.py`
   - `summarise_ccfinder.py`
   - `build_master_table.py`
   - `build_sample_status.py`
3. Wrap each tool as an isolated DSL2 module.
4. Build the per-sample subworkflow.
5. Build the cohort-level ANI and taxonomy subworkflows.
6. Integrate final aggregation.
7. Add stub runs and regression tests.
8. Compare against legacy notebook outputs on a frozen reference cohort.

---

## 17. References

These references document the major tools and concepts involved; the operational rules in this manual are still governed first by the locked user specification above.

- CheckM2: Chklovski A, et al. *Nature Methods* (2023). DOI: **10.1038/s41592-023-01940-w**
- FastANI: Jain C, et al. *Nature Communications* (2018). DOI: **10.1038/s41467-018-07641-9**
- BUSCO v5: Manni M, et al. *Molecular Biology and Evolution* (2021). DOI: **10.1093/molbev/msab199**
- Prokka: Seemann T. *Bioinformatics* (2014). DOI: **10.1093/bioinformatics/btu153**
- PADLOC: Payne LJ, et al. *Nucleic Acids Research* (2021). DOI: **10.1093/nar/gkab883**
- eggNOG-mapper v2: Cantalapiedra CP, et al. *Molecular Biology and Evolution* (2021). DOI: **10.1093/molbev/msab293**
- CRISPRCasFinder / related ecosystem: Couvin D, et al. *Nucleic Acids Research* (2018). DOI: **10.1093/nar/gky425**

---

## 18. Final instruction to Codex

Implement exactly this v1 specification unless a later change log explicitly overrides it.

When the specification conflicts with old notebook comments, prefer:

1. the locked user clarifications in this manual;
2. the actual legacy notebook code behaviour;
3. tool defaults.

In other words: **prefer explicit rules over folklore, and prefer reproducible ugliness over elegant nonsense.**
