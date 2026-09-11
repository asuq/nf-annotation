# syntax=docker/dockerfile:1
FROM ghcr.io/prefix-dev/pixi:0.76.2@sha256:8b206ef57005a902cb53f50dbaa47893a4038ca269f0b00038b51f18b1313cd4

WORKDIR /opt/nf-annotation/eggnog
COPY pixi.toml pixi.lock ./
RUN --mount=type=cache,target=/root/.cache/rattler \
    --mount=type=cache,target=/root/.cache/uv \
    pixi install --locked

ENV PATH="/opt/nf-annotation/eggnog/.pixi/envs/default/bin:${PATH}"
COPY install_go_export.py go_namespace_export.py verify_go_export.py ./
RUN python install_go_export.py \
    && python -c 'import eggnogmapper.annotator._codec, eggnogmapper.annotator._collect_inner' \
    && python verify_go_export.py \
    && emapper.py --version \
    && diamond version \
    && hmmsearch -h >/dev/null

WORKDIR /data
