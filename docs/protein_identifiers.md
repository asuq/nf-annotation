# Protein identifiers

Canonical protein FASTAs use readable IDs such as:

```text
>GCF_000005845_2__p000001
```

The prefix identifies the genome. Accession periods become underscores. Other reserved characters are escaped;
identifiers are never truncated or replaced by hashes. If two accessions would
produce the same readable prefix, the cohort is rejected before analysis rather
than receiving opaque collision suffixes.

Protein numbers start at one within each genome. Ordering follows source-contig
order, then CDS coordinates, strand and original protein ID. Reordering an input
protein FASTA does not renumber unchanged features.

The protein manifest retains the mapping:

- `tool_id`: the readable canonical FASTA/query ID.
- `protein_id`: the exact original protein/locus ID.
- `gene_id`: the accession-qualified source gene key.

Raw source annotations and identifiers remain unchanged. Sequence checksums and
method fingerprints remain hashes because they establish provenance, not protein
names. The v0.5 identifier change requires newly generated canonical bundles;
there is no old-hash ID alias or migration adapter.
