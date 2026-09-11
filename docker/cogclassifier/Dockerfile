# syntax=docker/dockerfile:1
FROM ghcr.io/prefix-dev/pixi:0.76.2@sha256:8b206ef57005a902cb53f50dbaa47893a4038ca269f0b00038b51f18b1313cd4

WORKDIR /opt/nf-annotation/cogclassifier
COPY pixi.toml pixi.lock ./
RUN --mount=type=cache,target=/root/.cache/rattler \
    --mount=type=cache,target=/root/.cache/uv \
    pixi install --locked

ENV PATH="/opt/nf-annotation/cogclassifier/.pixi/envs/default/bin:${PATH}"
RUN COGclassifier --version && rpsblast -version
WORKDIR /data
