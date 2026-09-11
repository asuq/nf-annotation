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

Runtime builds and these controls pass locally. Compatible full database
preparation, workflow wiring and the five-genome HPC qualification remain
v0.4 release gates; the development tags above are not published releases.
