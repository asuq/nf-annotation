FROM python:3.12-slim

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        aria2 \
        ca-certificates \
        procps \
        tar \
        unzip \
    && rm -rf /var/lib/apt/lists/*

COPY bin/prepare_runtime_databases.py /usr/local/bin/prepare_runtime_databases.py
COPY bin/annotation_resources.py bin/download_annotation_resources.py /usr/local/bin/
COPY bin/finalise_runtime_database.py /usr/local/bin/finalise_runtime_database.py
COPY bin/merge_runtime_database_reports.py /usr/local/bin/merge_runtime_database_reports.py
COPY assets/runtime/runtime_database_sources.json /opt/nf-myco_annotation/runtime_database_sources.json

RUN chmod +x /usr/local/bin/prepare_runtime_databases.py /usr/local/bin/finalise_runtime_database.py /usr/local/bin/merge_runtime_database_reports.py

ENV NF_MYCO_RUNTIME_DB_SOURCE_MANIFEST=/opt/nf-myco_annotation/runtime_database_sources.json
