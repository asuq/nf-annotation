# syntax=docker/dockerfile:1
FROM ghcr.io/prefix-dev/pixi:0.76.2@sha256:8b206ef57005a902cb53f50dbaa47893a4038ca269f0b00038b51f18b1313cd4

WORKDIR /opt/nf-annotation/kofam
COPY pixi.toml pixi.lock ./
RUN --mount=type=cache,target=/root/.cache/rattler pixi install --locked
ADD --checksum=sha256:ac4d025d9642d5a979ea10a18cac11df5e468f91eb11e19e1c6ef31e441d3d94 \
    https://github.com/takaram/kofam_scan/archive/62cee3943aa28b96f498150f5e0bb2d7c498e648.tar.gz /opt/nf-annotation/provenance/kofam_scan.tar.gz
RUN mkdir source && tar -xzf /opt/nf-annotation/provenance/kofam_scan.tar.gz --strip-components=1 -C source

ENV PATH="/opt/nf-annotation/kofam/.pixi/envs/default/bin:/opt/nf-annotation/kofam/source:${PATH}"
RUN exec_annotation --version && hmmsearch -h >/dev/null && parallel --version
WORKDIR /data
