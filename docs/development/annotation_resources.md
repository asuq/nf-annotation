# Annotation resource preparation

This documents the v0.4 preparation implementation. Full biological and
workflow qualification remains governed by [the release specification](v0.4.md).

The canonical source catalogue is
[`assets/runtime/runtime_database_sources.json`](../../assets/runtime/runtime_database_sources.json).
The superseded root-level copy and eggNOG v2 download branch have been removed.

| Component | Selected source | Preparation checks |
| --- | --- | --- |
| eggNOG | eggNOG 7 resources for emapper 3.0 | SQLite schema/integrity, DIAMOND database, GO namespaces, and full comparisons of the native taxid/field-presence caches against the database and ontology |
| COGclassifier | COG2024 definitions from pinned COGclassifier 2.0.0, CDD Cog index and ID mapping | Native RPS-BLAST index reader, COG definitions and complete functional-letter vocabulary |
| Pfam | Pfam 38.2 | Provider MD5s, matched HMM/metadata identities, finite gathering thresholds, HMMER parsing and pressed indexes |
| KOfam | Archive labelled 2026-07-02 | Exact native prokaryotic subset, adaptive thresholds/score types, selected-profile checksums and native parsing of every selected HMM |
| PADLOC | Database v2.0.0 at a fixed commit | Bundled in its dedicated image; source and compiled library checksums recorded during the build |

A source URL or dated directory alone is not a content pin. Acquisition records
store the actual byte count, strong ETag when available, upstream checksum when
provided, content SHA-256 and acquisition time. Pfam's MD5 values come from its
[release checksum file](https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.2/md5_checksums).
The eggNOG server did not provide the manifest described by its downloader when
inspected; its files require recorded content identities and explicit
compatibility checks. The KOfam archive's modification dates differ from its
directory label, so the acquisition record preserves the served metadata.

Malformed optional definition fields are retained with explicit field errors.
The pinned COGclassifier definitions contain a trailing space in the category value
for `COG6144`. Preparation retains the source bytes and records this in
`definition_field_errors.tsv`; it does not silently trim it into a valid category.
A COG assignment can remain usable while any affected category summaries must
be unavailable under the v0.4 field-error policy. This propagation is a required
part of annotation normalization.

Acquisition is separate from readiness. `bin/download_annotation_resources.py`
can download the declared files to an acquisition directory; it never writes a
ready marker. Transfers identify the nf-annotation client. A resume requires a
recorded partial plus a matching strong ETag or pinned upstream checksum. The
helper uses conditional requests and verifies the final size and checksum.
Changed or unrecorded partials fail explicitly. The HTTP response may omit its
own Content-Length; the total must still equal the advertised object size.
There is no automatic transfer retry loop.

`prepare_databases.nf` delegates each annotation resource to
`bin/prepare_runtime_databases.py` in that tool's image. For an existing
acquisition, the same CLI accepts `--pfam-source`, `--cogclassifier-source`,
`--kofam-source` or `--eggnog-source` and the matching `--*-dest` argument.
It also accepts a validated prepared resource directory as a source. It copies
prepared inputs, performs the native checks, and writes
`annotation_resource.json` before the canonical ready marker. Failed builds
remain in a `.preparing` directory for inspection. They are never promoted to
ready data or automatically overwritten.

`annotation_resource.json` inventories the files and defines their content
identity independently of the destination path. A changed source cannot reuse
a destination under the same version label. Whole-resource checksum validation
belongs in preparation or once-per-run preflight, not in every proteome task.
The ontology supports namespace validation and definitions; preparation adds no
GO assignments or ancestor expansion.

The acquisition and immutable-resource controls, four-resource Nextflow stub,
full native resource preparation and positive/zero searches pass. The
[qualification record](v0.4_qualification.md) lists the actual resource and SIF
identities and distinguishes these controls from pending full-cohort acceptance.
