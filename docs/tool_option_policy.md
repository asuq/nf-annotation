# Uniform tool-option policy

All samples being compared or combined use a shared invocation policy for each
tool. The workflow must not select rescue flags, sensitivity, thresholds or
alternate command modes by accession, assembly size, QC result or previous
failure. The same rule applies to new samples and reused evidence.

Sample paths, output identifiers, declared domain, resolved genetic code and
configured domain lineages are inputs to that shared policy. They are not
accession-specific command patches. Run-level extra arguments apply to every
eligible sample. CPU, memory and temporary paths follow common scheduler rules;
resource escalation on retry does not change scientific settings.

Prokka has no pre-QC exclusion gate and no sample-specific `--meta` or
`--metagenome` option. Its released implementation and existing invocation flags
remain unchanged. The workflow may observe execution and report failure, but
must not modify third-party source, replace predictions or rerun a failed sample
under a different scientific method as an automatic rescue.

## Repository audit

The audit covers production command construction in `modules/local/`,
`subworkflows/local/`, `lib/`, `conf/`, and the native command/execution helpers
in `bin/`. Synthetic qualification fixtures are not production commands.

| Tool or path | Shared settings and permitted input/resource differences |
| --- | --- |
| Prokka | One `--compliant --rfam` command template; kingdom and code come from the common resolver, CPU count from the task allocation. No quality/size/accession branch changes flags. |
| CheckM2 | Every sample receives the same two invocations, for tables 4 and 11. Selection between their reports is separate from execution. |
| Barrnap | One template; the declared domain supplies `bac` or `arc`. |
| BUSCO | Genome/offline/tar invocation for configured lineages. The same domain-profile mapping selects each sample's lineage set. |
| Codetta and CRISPRCasFinder | Fixed templates; extra arguments are run-level. CRISPR receives the commonly resolved code and all inputs receive the same reversible header staging. |
| eggNOG, COGclassifier, Pfam, KOfam and PADLOC | Commands are built centrally from tool and resource allocation, not accession. Sensitivity, donor policy and thresholds are shared. eggNOG search and annotation phases are distinct documented operations, not sample-specific methods. |
| FastANI and FastAAI | Shared cohort commands and thresholds. Domain grouping, eligibility and profile reuse change the work required, not the method applied to comparable samples. |
| Retry/site profiles | No accession-specific scientific flags found. CPU/memory/time may increase under common retry rules; input staging and storage location may differ by site. |

The audit found no production accession-specific option overrides. Existing
checks reject mixed annotation methods and mismatched execution resources.
The declared domain/code/lineage inputs and resource rules above are intentional
differences, not a promise that every complete argument vector is byte-identical.

Existing pinned images include previously recorded CRISPR runtime and
eggNOG/Pfam/FastAAI exporter adaptations. They are uniform dependencies, not
per-sample workarounds. This change leaves those image contents frozen and adds
no upstream source patch. External project overlays and arbitrary user-supplied
Nextflow configurations are outside the repository audit and must obey the same
policy before their results are combined.
