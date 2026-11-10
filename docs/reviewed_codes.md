# Reviewed genetic-code choices

Automatic paired-CheckM2 selection remains the default. To apply an explicitly
reviewed choice, pass `--reviewed_codes reviews.json` to `main.nf`:

```json
{
  "schema": "reviewed-genetic-codes-v1",
  "reviews": [
    {
      "accession": "EXAMPLE_1",
      "fasta_identity_md5": "0123456789abcdef0123456789abcdef",
      "genetic_code": 4,
      "review_id": "review-record-id",
      "reason": "Evidence supporting the reviewed choice"
    }
  ]
}
```

The digest above is illustrative. Calculate the actual identity with the pinned
pipeline environment:

```bash
PYTHONPATH=bin pixi run --locked python -c \
  'from pathlib import Path; from reviewed_codes import fasta_identity; print(fasta_identity(Path("genome.fasta")))'
```

Identity includes ordered contig IDs and case-sensitive DNA. Descriptions,
wrapping and gzip compression do not affect it. MD5 is an identity check, not a
security signature. Each accession/identity pair must be unique; only integer
codes 4 and 11 are supported. A supplied accession with no matching FASTA identity
fails explicitly. Accessions absent from the manifest retain automatic selection.

A reviewed choice requires valid, matching CheckM2 reports for both translation
tables. It cannot override a failed report pair. `Gcode_Automatic` and
`Gcode_Automatic_Reason` retain the automatic result; `Gcode` is the effective
choice. `Gcode_Review` records the exact review entry, and
`Gcode_Selection_Reason` becomes `reviewed_code_manifest`. Quality is calculated
from the effective code's report. Downstream native tools receive that code.

Use the same manifest for cohort reuse. A changed review requires regeneration
of the affected sample. Annotation-only reuse needs no new declaration: it
replays the recorded review against the retained genome and paired reports.
No accession-specific overrides are embedded in the workflow. These changes
require new QC publications; old publication/cache compatibility is not promised.
