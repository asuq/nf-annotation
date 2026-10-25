# Explicit domains and separate QC/ANI stages

Every `main.nf` run requires an explicit `--domain` selection, supplied on the
command line or as `params.domain` in a configuration file. There is no default.
The same setting applies to full, QC-only and ANI-only runs.

| Selection | Domain source |
| --- | --- |
| `--domain bacteria` | The whole cohort is bacterial. |
| `--domain archaea` | The whole cohort is archaeal. |
| `--domain metadata` | Read each sample's domain from the exact `Domain` metadata column. |

In metadata mode, every requested sample must have a `Domain` value of `Bacteria`
or `Archaea` (case-insensitive). For example:

```text
accession,Domain
sample_a,Bacteria
sample_b,Archaea
```

Missing, unknown or ambiguous values fail before per-sample QC. In single-domain
mode, a populated `Domain` field must agree with the selected domain. Organism
names, accession prefixes and metabolic groups are never classifiers. Prokka and
Barrnap kingdoms are derived internally; they have no separate CLI options.

The native four-column sample CSV and identifier/mapping rules are unchanged.
Changing a domain can alter biological predictions and requires new QC; native
gene/protein IDs are never rewritten to match genome IDs.

## BUSCO profiles

For metadata mode or an archaeal cohort, supply explicit profiles. A mixed-cohort
configuration is:

```groovy
params {
    domain = 'metadata'
    busco_domain_lineages = [
        Bacteria: ['bacteria_odb12'],
        Archaea: ['archaea_odb12']
    ]
    busco_db = '/path/to/prepared/busco'
}
```

A profile is an ordered, nonempty list of lineage dataset names. Its first
lineage supplies the existing ANI eligibility and representative-scoring BUSCO
value; additional entries are supplementary assessments. Choose biologically
appropriate lineages yourself. The pipeline does not infer subclades, select
whichever lineage yields the largest score, or verify an arbitrary custom
lineage's phylogenetic placement. It rejects shared/duplicate profile lineages
and a bacterial/archaeal root dataset explicitly assigned to the opposite domain.

For a single archaeal cohort, use `--domain archaea` and an `Archaea` BUSCO
profile. For a bacterial cohort, `--domain bacteria` supports either an explicit
`Bacteria` profile or the existing `busco_lineages` and `busco_primary_column`
configuration. A cohort-wide `busco_primary_column` cannot be combined with
domain profiles. Configuring profiles does not supply the required domain choice.

The resolved settings are:

| Domain | Prokka `--kingdom` | Barrnap `--kingdom` | BUSCO |
| --- | --- | --- | --- |
| Bacteria | `Bacteria` | `bac` | Only the declared bacterial profile |
| Archaea | `Archaea` | `arc` | Only the declared archaeal profile |

The existing bacterial BUSCO selection remains effective with `--domain bacteria`
when no domain profiles are configured. QC thresholds, gcode rules and ANI
settings are unchanged.

CheckM2 still evaluates both codes 4 and 11, using the selected existing gcode
rule. Prokka receives that selected code explicitly, including when the declared
domain is Archaea. Missing-code samples still skip gene prediction. Codetta
remains an independent inference/reporting tool and does not rescue missing
CheckM2 code assignments.

ANI distances, clustering thresholds, representative-scoring weights and the
16S policy are unchanged. The change selects the correct declared BUSCO result
for each sample; it does not introduce a domain-specific ANI formula. Missing
primary BUSCO evidence remains unavailable and excludes that sample under the
existing eligibility rules. A lineage not requested for a sample has `NA`
evidence and `na` status rather than a fabricated zero or failed analysis.

Domain-aware master tables add `Annotation_Domain`, `Annotation_Domain_Source`,
`Primary_BUSCO_Lineage` and `Primary_BUSCO_Score`. ANI metadata carries the
selected compact result as `BUSCO_selected`, alongside `annotation_domain` and
`primary_busco_lineage`. The original per-lineage BUSCO columns remain available.
These resolved master fields are reserved and must not overlap input metadata.

`prepare_databases.nf` resolves the union of the configured domain profiles.
It does not require `--domain`, because it prepares resources rather than
analysing samples. Likewise, `reannotate.nf` uses recorded protein evidence and
does not require a new domain declaration.
Keep that configuration consistent with the execution configuration. The same
prepared datasets and native BUSCO genome/offline method are used; no automatic
lineage search or implicit database update is introduced.

## Execution stages

| `pipeline_stage` | New analyses | Publication/use |
| --- | --- | --- |
| `full` (default) | Existing QC, gene prediction, CRISPR, selected functional tools and ANI | Complete existing workflow |
| `qc` | Barrnap, paired CheckM2, BUSCO, Codetta and code-qualified Prokka | Reusable QC/protein-bundle publications; no CRISPR, functional search or ANI |
| `ani` | Cumulative ANI, representative selection and reporting from completed publications | Requires `update_from`; cannot submit new QC or gene prediction |

The QC and ANI stages disable the separate functional callers even when
`annotation_tools` is configured. Prokka retains its existing annotation flags,
including `--rfam` and `--compliant`; QC does not rewrite its native output IDs.
Protein and coordinate mapping tables continue to connect genes to accessions.
CRISPR was not requested for new QC-stage samples, so their status is `na` and
their CRISPR evidence remains `NA`. `QC_Scope` records that scope in the master
table. QC-only publications retain the normal complete publication manifest and
report ANI fields as unavailable, with `ani_included = na`.

For example, with paths and domain profiles in `analysis.config`:

```bash
nextflow run main.nf --domain metadata -profile oist -c analysis.config \
    --pipeline_stage qc --sample_csv batch.csv --metadata metadata.tsv \
    --outdir qc_batch

nextflow run main.nf --domain metadata -profile oist -c analysis.config \
    --pipeline_stage ani --sample_csv cumulative.csv --metadata metadata.tsv \
    --update_from qc_batch,previous_qc_batch --ani_16s_policy ignore \
    --outdir cumulative_ani
```

The cumulative manifest must list every desired genome, using the existing
accession and genome-file contracts. Every requested sample must occur in one
completed source publication, with matching genome content and policy; an absent
sample is an error, not permission to run fresh QC. Source directories must be
distinct from the new output directory. `ignore` bypasses only 16S eligibility;
the quality screen and other ANI exclusions remain active.

The QC stage still needs its normal QC databases and taxonomy input. An ANI-only
update can use the published QC data without supplying CheckM2/Codetta resources
for new analyses, because no new samples are allowed. It still needs the cohort
metadata and taxonomy input for reporting.

## Reuse boundaries and provenance

Domain, declaration source, kingdom choices, ordered BUSCO profile and primary
lineage survive publication and reuse. Changing these settings, or attempting to
assign them to an old publication that did not record them, requires new QC.
Legacy bacterial publications remain reusable with `--domain bacteria` and
their existing BUSCO configuration; this is not an
importer for previously patched or unverified domain assignments.

QC publications support ANI-only reuse and the existing `reannotate.nf` workflow
for functional annotation. Full-stage reuse of a QC-only publication is rejected
because its CRISPR evidence was never requested. Existing full publications can
supply their retained evidence to QC/ANI updates without rerunning those tools.
The original work directory is not needed. Policy-helper code hashes are explicit
validation/reporting task inputs, so changed code is revalidated on resume.
Reused native outputs, identifiers and mapping tables are preserved; no gene-renaming adapter is introduced.

## Codetta staging

The reference directory is staged under the fixed name `codetta_db`, separately
from the `codetta` output directory. A symlink or non-directory at the output
location is rejected before writing. This fixes a reference directory literally
named `codetta` without changing Codetta's command, reference data, thresholds
or published scientific output formats.

## Scientific references

- [Prokka documentation](https://github.com/tseemann/prokka): explicit bacterial
  and archaeal kingdom modes, independent explicit genetic-code selection.
- [Barrnap 0.9 documentation](https://github.com/tseemann/barrnap/tree/0.9): the
  `bac` and `arc` rRNA models used by the pinned runtime.
- [BUSCO user guide](https://busco.ezlab.org/busco_userguide): select an appropriate
  lineage for the organism; matches against an unrelated domain do not establish
  biological completeness. These interfaces retain the pinned BUSCO runtime and
  prepared dataset versions rather than switching to the current guide's newest
  software or OrthoDB release.
