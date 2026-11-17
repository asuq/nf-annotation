---
type: decision
title: Assembly length and downstream eligibility
description: Preserve native predictions and ordered master rows while excluding assemblies below 20,000 bp from cohort analyses.
status: stable
generated:
  by: openai/gpt-6.2-astra
  at: 2026-11-17T21:53:00Z
verified:
  by: openai/gpt-6.2-astra
  at: 2026-11-17T21:53:00Z
stale_after: 2027-02-17T00:00:00Z
sources:
  - id: policy
    resource: ../docs/prokka_prediction.md
  - id: boundary
    resource: ../tests/test_short_assembly_integration.py
---

Use the sum of all staged contig lengths. A genome with 19,999 bp is low quality
under this rule and carries the downstream exclusion reason
`assembly_below_20000_bp`. At 20,000 bp its original QC criteria determine
eligibility. The master row retains the CheckM2 classification, selected code
and native evidence. Prokka runs under the released invocation policy.[^policy]

The boundary test exercises both successful and failed predictions below the
cutoff and follows the QC publication into ANI and AAI clustering.[^boundary]

[^policy]: [Prediction and eligibility policy](../docs/prokka_prediction.md)
[^boundary]: [QC and clustering integration test](../tests/test_short_assembly_integration.py)
