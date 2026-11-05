# syntax=docker/dockerfile:1
FROM ghcr.io/prefix-dev/pixi:0.76.2@sha256:8b206ef57005a902cb53f50dbaa47893a4038ca269f0b00038b51f18b1313cd4
WORKDIR /opt/nf-annotation/fastaai
COPY pixi.toml pixi.lock ./
RUN --mount=type=cache,target=/root/.cache/rattler pixi install --locked
ADD --checksum=sha256:4f4b36102d148e1f47a52a59747e98ebfab69b20e4ac7de45aa7a9358f84975f \
    https://codeload.github.com/cruizperez/FastAAI/tar.gz/3b4904a223590f8227e9c791b6b7ed6a09c67757 /opt/nf-annotation/provenance/fastaai.tar.gz
RUN mkdir source && tar -xzf /opt/nf-annotation/provenance/fastaai.tar.gz --strip-components=1 -C source \
    && pixi run python -m pip install --no-deps --no-build-isolation ./source
ENV PATH="/opt/nf-annotation/fastaai/.pixi/envs/default/bin:${PATH}"
COPY install_export.py native_binary_export.py /opt/nf-annotation/provenance/
RUN python /opt/nf-annotation/provenance/install_export.py \
    && python -c "import fastaai.fastaai as f; assert __import__('pathlib').Path(f.find_hmm()).is_file()"
WORKDIR /data
