# Prokka prediction outcomes

The pipeline still attempts Prokka for genetic codes 4 and 11, including
low-quality QC inputs. Quality-based clustering exclusion is separate from
collecting prediction evidence. CheckM2's own CDS counts are retained separately
and never substituted for Prokka proteins.

## Released software and execution

Every sample uses the same Prokka option set. No quality gate or
`--meta`/`--metagenome` option is added. The released Prokka 1.15.6 code remains
unchanged, including its mode-selection bug: a defined empty training-file
setting can select Prodigal single mode for short input. The workflow does not
override that choice. Existing `--compliant`, `--rfam`, kingdom/code/CPU arguments,
the 200-bp contig filter and native Prodigal `-c -m` arguments remain unchanged.

A workflow-owned observer records the actual Prodigal child outcome while
forwarding the exact argument list and prediction stream to the released binary.
The observer applies to every Prokka invocation. Neither installed software nor
scientific output is modified. A failed child remains a prediction failure even
if Prokka later reports zero CDS and exits zero. The execution receipt retains
parent and child outcomes. Only successful parent and child execution can produce
the `NF_PROKKA_PREDICTION_V1` record used by downstream reporting.

## Zero, unavailable and failed

`prokka_status=done` requires a successful recorded prediction with proteins.
`prokka_status=no_cds` identifies a successful recorded prediction with zero CDS
and an empty protein FASTA. It does not claim biological absence of coding genes.
Nonzero execution, missing prediction provenance or contradictory output remains
`failed`; missing evidence never becomes a zero count.

Bundle construction also checks the retained DNA, GFF and GenBank. A consistent
zero-CDS result records `status=incompatible_input`, `reason=no_cds` and
`protein_count=0`: it is unsuitable for protein-based analysis, not a prediction
failure. Enabled functional tools keep that explicit unavailable-input outcome;
they do not report zero biological hits. Failed predictions use `upstream_failed`
and have no known protein count. Genuine zero counts remain zero in reporting.
Reuse preserves and checks these distinctions.

QC-only `annotation_acceptance.complete` describes completion of the disabled
functional-annotation publication. It does not certify that every sample passed
QC or produced usable proteins. Consumers must retain the per-sample QC,
prediction and bundle outcomes, including failed members of an otherwise useful
batch. Independent samples continue after a prediction failure.

No old output is relabelled or overwritten to appear qualified under the new
observation contract. Existing runs and their original diagnostics remain intact.

Fresh QC publications are required for reuse under this outcome contract.
Older Prokka logs lack the recorded child outcome and cannot establish a verified
zero or successful prediction through this checker. Do not rewrite old logs to
add receipts.
