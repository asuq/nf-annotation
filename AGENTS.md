# Repository policy: uniform tool invocation

- All samples being processed or compared together must use the same pinned tool
  versions, command templates and scientific option settings. Do not add
  accession-specific exceptions, rescue flags or alternate modes based on sample
  identity, assembly size, quality or a failed attempt.
- Input/output paths and identifiers vary by sample. Biological input values
  already supported by the workflow (declared domain, resolved genetic code and
  configured domain lineage) follow the same documented resolution rules for
  every sample. Reviewed genetic-code choices are recorded input data, not
  sample-specific command patches.
- Run-level method settings apply uniformly. Resource arguments follow the same
  scheduler-owned allocation rules. Eligibility filters, requested stages and
  reuse may determine which work is needed; they must not silently change the
  scientific options of work that executes. Retries keep scientific settings;
  resource escalation follows the same declared scheduler policy.
- Prokka keeps the same option set for every sample. Do not add a quality gate
  before Prokka or pass `--meta`/`--metagenome` to rescue particular assemblies.
- Keep the pinned released software unchanged, including known upstream bugs.
  Do not patch third-party source as a recovery mechanism. Workflow-owned
  execution observation may record child status, but must preserve tool arguments,
  input bytes and predictions. Existing pinned runtime contents are the baseline;
  this policy does not authorise changing their recorded adaptations.
- Missing evidence, failed prediction and a verified zero result remain distinct.
  Test the shared invocation and outcome contract; do not repair scientific output
  or fabricate a successful bundle for a failed member.

See [tool-option policy and audit](docs/tool_option_policy.md) and
[Prokka prediction outcomes](docs/prokka_prediction.md).
