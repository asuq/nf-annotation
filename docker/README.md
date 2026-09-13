# Docker Images

Module-specific Dockerfiles live under `docker/<tool>/Dockerfile` when
BioContainers do not cover a required runtime.

Current repo-owned images:

- `docker/ccfinder/Dockerfile`: CRISPRCasFinder runtime
- `docker/codetta/Dockerfile`: Codetta runtime
- `docker/eggnog/Dockerfile`: eggNOG-mapper v3 beta6, with an export of the
  native GO terms and confidence for each namespace
- `docker/cogclassifier/Dockerfile`: COGclassifier 2.0.0 and RPS-BLAST
- `docker/pfam/Dockerfile`: PfamScan 1.6 and HMMER, retaining search tables and
  results before and after the native clan-overlap resolver
- `docker/kofam/Dockerfile`: KOfamScan 1.3.0 and HMMER
- `docker/padloc/Dockerfile`: PADLOC runtime with a fixed launcher and bundled database
- `docker/python_helper/Dockerfile`: shared Python helper runtime used by
  `params.python_container`, including the ANI scientific stack (`numpy` and
  `scipy`)

## v0.4 annotation runtimes

The five annotation tools have separate images. The eggNOG, COGclassifier,
Pfam and KOfam images use their own committed Linux x86-64 Pixi manifest and
lock file. Their Docker build context is the corresponding tool directory:

```bash
docker build --platform linux/amd64 -t nf-annotation-eggnog:v04-dev docker/eggnog
docker build --platform linux/amd64 -t nf-annotation-cogclassifier:v04-dev docker/cogclassifier
docker build --platform linux/amd64 -t nf-annotation-pfam:v04-dev docker/pfam
docker build --platform linux/amd64 -t nf-annotation-kofam:v04-dev docker/kofam
docker build --platform linux/amd64 -t nf-annotation-padloc:v04-dev docker/padloc
```

The base images are pinned by digest. Source archives are pinned by commit or
checksum. PADLOC bundles database v2.0.0 at commit
`7f99b47b75e232b111c18626badb9ac32e8e0b5a`; its native HMM/CM concatenation is
performed in the C locale and recorded in `/opt/nf-annotation/provenance`.
No database update runs at image build time. The other annotation databases
are prepared separately; building an image does not qualify those resources.

eggNOG preserves the native annotation table and adds
`<annotations>.go_namespaces.tsv`. Its source adaptation only exports the
already computed `gos_mf`, `gos_bp`, `gos_cc` values and their native confidence.
Build-time controls compare the original and adapted eager, lazy and masked
engines, multiprocessing, and native writer bytes. The pipeline uses the
`closest` donor pool; tool-level `--resume` is rejected because it cannot
guarantee a complete sidecar. Nextflow task reuse remains separate.

Pfam's `pfam_scan_evidence.pl` runs the native gathering-threshold search once,
retains HMMER tables and all clan alternatives, then applies the upstream
resolver. Build-time synthetic HMM controls compare its resolved output with
the unmodified upstream program and cover repeated domains, nested domains,
overlapping alternatives and zero hits. Original adapted source files and
checksums remain in `/opt/nf-annotation/provenance`.

Runtime builds, native controls and the five-genome HPC qualification pass.
The development tags above are local build names. The v0.4.0 release publishes
the same qualified image payloads under `quay.io/asuq1617/nf-annotation-<tool>:0.4.0`;
execution uses immutable Linux amd64 manifest digests below, not those tags.
The helper digest matches the OCI manifest used for HPC qualification.

| Runtime | Immutable reference |
| --- | --- |
| python | `quay.io/asuq1617/python-scipy@sha256:30051a8fbc1fddbd6c5e1c10760699f40d2f5711e6640c36d4747a1882ef05c2` |
| eggnog | `quay.io/asuq1617/nf-annotation-eggnog@sha256:4eb430644be8a94e4e752874dc8d646cdc9acb9ddf28c5b2d1bb283987a1027d` |
| cogclassifier | `quay.io/asuq1617/nf-annotation-cogclassifier@sha256:d3900ea264815e527dac4a1c2737a270258bcf6e813cc606d348f7ef674ae10c` |
| pfam | `quay.io/asuq1617/nf-annotation-pfam@sha256:c6ef9dcabe7c494e853a8f3eef3a6ecd987b2a5abafe07ad13df8b406834f2f6` |
| kofam | `quay.io/asuq1617/nf-annotation-kofam@sha256:bfe33293b9d32c6671c20672472d669428458527bd38de958eae2fcb6ed4fcb0` |
| padloc | `quay.io/asuq1617/nf-annotation-padloc@sha256:557718303da941c8f6ac6c55bdb659f0e4c37963657a7da4d4a33e47ae3ccbd8` |

These defaults can be overridden with checksum-verified local SIF paths for
offline HPC use. OCI digests and SIF file checksums are different identities;
retain conversion provenance when preparing SIF files.
