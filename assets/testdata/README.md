# Test data

- `../fixtures/stub/`: miniature genomes and database placeholders for stub runs.
- `../tables/stub/`: the matching sample sheet and metadata.
- `../tables/headers/`: empty table templates for gated branches.
- `../tables/contracts/`: authoritative report column contracts.
- `../tables/acceptance/` and `../tables/medium/`: frozen cohort descriptors.
- `fastaai/`: native proteome identities and frozen ANI scoring controls.

Downloaded data and generated acceptance outputs remain under ignored `local/`.
Accession validation and internal-ID mapping are implemented in `bin/validate_inputs.py`.
