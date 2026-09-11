# syntax=docker/dockerfile:1
FROM quay.io/biocontainers/padloc:2.0.0--hdfd78af_1@sha256:57dc14167c536f17c4ef2c9aed9ffa972a68a1087640e9079f128c4fd4493aca

ENV PADLOC_DATA_DIR=/opt/padloc-data
ENV LC_ALL=C

ADD --checksum=sha256:de1c9697ee735c7c8a856295ace5614c9c1959b6010cbc074a29125aa9c04fc9 \
    https://github.com/padlocbio/padloc-db/archive/7f99b47b75e232b111c18626badb9ac32e8e0b5a.tar.gz /opt/nf-annotation/provenance/padloc-db.tar.gz

RUN PADLOC_BIN="$(command -v padloc)" && \
    sed -i "s#mkdir -p \"\${SRC_DIR}/../data\"#mkdir -p \"${PADLOC_DATA_DIR}\"#" "${PADLOC_BIN}" && \
    sed -i "s#DATA=\$(normpath \"\${SRC_DIR}/../data\")#DATA=\$(normpath \"${PADLOC_DATA_DIR}\")#" "${PADLOC_BIN}" && \
    mkdir -p "${PADLOC_DATA_DIR}" && \
    tar -xzf /opt/nf-annotation/provenance/padloc-db.tar.gz --strip-components=1 -C "${PADLOC_DATA_DIR}" && \
    cat "${PADLOC_DATA_DIR}"/hmm/*.hmm > "${PADLOC_DATA_DIR}/padlocdb.hmm" && \
    cat "${PADLOC_DATA_DIR}"/cm/*.cm > "${PADLOC_DATA_DIR}/padlocdb.cm" && \
    mv "${PADLOC_DATA_DIR}/padlocdb.hmm" "${PADLOC_DATA_DIR}/hmm/padlocdb.hmm" && \
    mv "${PADLOC_DATA_DIR}/padlocdb.cm" "${PADLOC_DATA_DIR}/cm/padlocdb.cm" && \
    sha256sum "${PADLOC_DATA_DIR}/hmm/padlocdb.hmm" "${PADLOC_DATA_DIR}/cm/padlocdb.cm" \
      > /opt/nf-annotation/provenance/padloc-compiled.sha256 && \
    printf '%s\n' 'padloc-db v2.0.0' 'commit=7f99b47b75e232b111c18626badb9ac32e8e0b5a' \
      'archive_sha256=de1c9697ee735c7c8a856295ace5614c9c1959b6010cbc074a29125aa9c04fc9' \
      > /opt/nf-annotation/provenance/padloc-db.txt

RUN PADLOC_BIN="$(command -v padloc)" && \
    grep -F "mkdir -p \"${PADLOC_DATA_DIR}\"" "${PADLOC_BIN}" >/dev/null && \
    grep -F "DATA=\$(normpath \"${PADLOC_DATA_DIR}\")" "${PADLOC_BIN}" >/dev/null && \
    ! grep -F 'mkdir -p "${SRC_DIR}/../data"' "${PADLOC_BIN}" >/dev/null && \
    ! grep -F 'DATA=$(normpath "${SRC_DIR}/../data")' "${PADLOC_BIN}" >/dev/null && \
    test -f "${PADLOC_DATA_DIR}/hmm/padlocdb.hmm" && \
    padloc --version >/dev/null
