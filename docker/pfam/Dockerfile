# syntax=docker/dockerfile:1
FROM ghcr.io/prefix-dev/pixi:0.76.2@sha256:8b206ef57005a902cb53f50dbaa47893a4038ca269f0b00038b51f18b1313cd4

WORKDIR /opt/nf-annotation/pfam
COPY pixi.toml pixi.lock ./
RUN --mount=type=cache,target=/root/.cache/rattler pixi install --locked
ADD --checksum=sha256:365c96bc150d5057349c3016d62667c58cb33afcfb6329457ae16ab5aae4f401 \
    https://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/PfamScan.tar.gz /opt/nf-annotation/provenance/PfamScan.tar.gz
RUN tar -xzf /opt/nf-annotation/provenance/PfamScan.tar.gz

ENV PATH="/opt/nf-annotation/pfam/.pixi/envs/default/bin:/opt/nf-annotation/pfam/PfamScan:/opt/nf-annotation/pfam:${PATH}"
ENV PERL5LIB=/opt/nf-annotation/pfam/PfamScan
COPY install_evidence_export.py pfam_scan_evidence.pl verify_evidence_export.py ./
RUN python install_evidence_export.py \
    && chmod +x pfam_scan_evidence.pl \
    && perl -c pfam_scan_evidence.pl \
    && perl -c PfamScan/pfam_scan.pl \
    && hmmscan -h >/dev/null
RUN --network=none python verify_evidence_export.py
WORKDIR /data
