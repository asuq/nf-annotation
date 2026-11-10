#!/usr/bin/env python3
"""Acquire annotation resources with verifiable restart and content receipts.

This is the acquisition part of prepare_runtime_databases.py. It never writes a
ready marker: extraction, indexes and tool-specific validation are separate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import re
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Any
from urllib.request import Request, urlopen

LOGGER = logging.getLogger(__name__)
CHUNK_SIZE = 8 * 1024 * 1024
USER_AGENT = "nf-annotation/0.4"


class ResourceDownloadError(RuntimeError):
    """An acquisition cannot be completed with verified integrity."""


def checksum_file(path: Path, algorithm: str = "sha256") -> str:
    """Hash a file once using the specified integrity algorithm."""
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, algorithm).hexdigest()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    """Replace a receipt atomically after serializing the full record."""
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n")
    temporary.replace(path)


def expected_checksum(config: dict[str, Any] | None) -> tuple[str, str] | None:
    """Validate a pinned upstream checksum, without guessing its algorithm."""
    if config is None:
        return None
    algorithm, value = config.get("type"), config.get("value")
    lengths = {"md5": 32, "sha256": 64}
    if algorithm not in lengths or not isinstance(value, str):
        raise ResourceDownloadError("Require an explicit md5 or sha256 checksum value")
    if not re.fullmatch(rf"[0-9a-fA-F]{{{lengths[algorithm]}}}", value):
        raise ResourceDownloadError(f"Invalid {algorithm} checksum value")
    return algorithm, value.lower()


def acquire_file(
    url: str, destination: Path, checksum_config: dict[str, Any] | None = None
) -> dict[str, Any]:
    """Download one object, resuming only a recorded and verifiable partial.

    An upstream checksum permits verified range resumption. Without one, a
    strong ETag is required and sent with If-Match on every GET. A changed
    object, unrecorded partial, short response or checksum mismatch fails.
    No failed transfer is promoted to a complete input.
    """
    checksum = expected_checksum(checksum_config)
    destination.parent.mkdir(parents=True, exist_ok=True)
    receipt = destination.with_name(destination.name + ".acquisition.json")
    partial = destination.with_name(destination.name + ".part")
    metadata = partial.with_name(partial.name + ".json")
    request_identity = {"url": url, "checksum": list(checksum) if checksum else None}
    if destination.exists():
        if not receipt.is_file():
            raise ResourceDownloadError(f"Existing file has no receipt: {destination}")
        record = json.loads(receipt.read_text())
        if any(record.get(key) != value for key, value in request_identity.items()):
            raise ResourceDownloadError(
                f"Existing acquisition has a different source: {destination}"
            )
        if record.get("bytes") != destination.stat().st_size or record.get(
            "sha256"
        ) != checksum_file(destination):
            raise ResourceDownloadError(
                f"Existing acquisition has changed: {destination}"
            )
        return record

    with urlopen(
        Request(url, headers={"User-Agent": USER_AGENT}, method="HEAD"), timeout=60
    ) as response:
        raw_size = response.headers.get("Content-Length")
        if raw_size is None or not raw_size.isdecimal() or int(raw_size) <= 0:
            raise ResourceDownloadError(f"Missing positive object size: {url}")
        etag = response.headers.get("ETag")
        strong_etag = etag if etag and re.fullmatch(r'"[^"\r\n]+"', etag) else None
        source = {
            **request_identity,
            "bytes": int(raw_size),
            "etag": strong_etag,
            "last_modified": response.headers.get("Last-Modified"),
            "client": USER_AGENT,
        }
    if partial.exists():
        if not metadata.is_file():
            raise ResourceDownloadError(
                f"Unrecorded partial must be quarantined before retry: {partial}"
            )
        previous = json.loads(metadata.read_text())
        # Last-Modified is provenance only; it is not a strong validator.
        for key in ("url", "checksum", "bytes", "etag"):
            if previous.get(key) != source[key]:
                raise ResourceDownloadError(
                    f"Source changed during acquisition: {url} ({key})"
                )
        if not strong_etag and checksum is None:
            raise ResourceDownloadError(
                f"Cannot verify resumption without a strong ETag or checksum: {url}"
            )
        offset = partial.stat().st_size
        if offset > source["bytes"]:
            raise ResourceDownloadError(f"Partial exceeds source size: {partial}")
    else:
        offset = 0
        write_json(metadata, source)
    if offset < source["bytes"]:
        headers = {"Accept-Encoding": "identity", "User-Agent": USER_AGENT}
        if strong_etag:
            headers["If-Match"] = strong_etag
        if offset:
            headers["Range"] = f"bytes={offset}-"
        LOGGER.info(
            "Downloading %s (%d/%d bytes)", destination.name, offset, source["bytes"]
        )
        with urlopen(Request(url, headers=headers), timeout=60) as response:
            expected_status = 206 if offset else 200
            if response.status != expected_status:
                raise ResourceDownloadError(
                    f"Unexpected HTTP status {response.status}; expected {expected_status}: {url}"
                )
            if strong_etag and response.headers.get("ETag") != strong_etag:
                raise ResourceDownloadError(f"Response changed its strong ETag: {url}")
            if offset:
                expected_range = (
                    f"bytes {offset}-{source['bytes'] - 1}/{source['bytes']}"
                )
                if response.headers.get("Content-Range") != expected_range:
                    raise ResourceDownloadError(f"Unexpected Content-Range: {url}")
            response_length = response.headers.get("Content-Length")
            if response_length is not None and response_length != str(
                source["bytes"] - offset
            ):
                raise ResourceDownloadError(
                    f"Response length does not match the object size: {url}"
                )
            with partial.open("ab" if offset else "wb") as handle:
                while chunk := response.read(CHUNK_SIZE):
                    handle.write(chunk)
        if partial.stat().st_size != source["bytes"]:
            raise ResourceDownloadError(f"Incomplete download: {partial}")
    digests = {"sha256": hashlib.sha256()}
    if checksum:
        digests.setdefault(checksum[0], hashlib.new(checksum[0]))
    with partial.open("rb") as handle:
        while chunk := handle.read(CHUNK_SIZE):
            for digest in digests.values():
                digest.update(chunk)
    if checksum and digests[checksum[0]].hexdigest() != checksum[1]:
        raise ResourceDownloadError(
            f"Upstream {checksum[0]} checksum mismatch: {destination.name}"
        )
    record = {
        **source,
        "file": destination.name,
        "sha256": digests["sha256"].hexdigest(),
        "acquired_at": datetime.now(UTC).isoformat(),
    }
    # A crash between these operations leaves a complete partial that can be
    # revalidated. It never leaves a final input without a receipt.
    write_json(receipt, record)
    partial.replace(destination)
    metadata.unlink()
    return record


def acquire_component(
    component: str, version: str, config: dict[str, Any], destination: Path
) -> dict[str, Any]:
    """Freeze the declared inputs from the canonical runtime source manifest."""
    if config.get("kind") != "annotation_bundle":
        raise ResourceDownloadError(f"Not an annotation resource bundle: {component}")
    files = config.get("files")
    if not isinstance(files, list) or not files:
        raise ResourceDownloadError(f"Empty resource list: {component}")
    names = [entry.get("name") for entry in files]
    if any(
        not isinstance(name, str) or Path(name).name != name or name in (".", "..")
        for name in names
    ):
        raise ResourceDownloadError("Resource names must be plain filenames")
    if len(set(names)) != len(names):
        raise ResourceDownloadError("Duplicate resource names")
    if any(not entry.get("url", "").startswith("https://") for entry in files):
        raise ResourceDownloadError("Resource sources must use HTTPS")
    destination.mkdir(parents=True, exist_ok=True)
    lock = destination / ".acquisition.lock"
    try:
        handle = lock.open("x")
    except FileExistsError as error:
        raise ResourceDownloadError(
            f"Acquisition lock already exists: {lock}"
        ) from error
    try:
        with handle:
            handle.write(f"{component}\n")
        records = [
            acquire_file(
                entry["url"], destination / entry["name"], entry.get("checksum")
            )
            for entry in files
        ]
        payload = {"component": component, "version": version, "files": records}
        write_json(destination / "acquisition.json", payload)
        return payload
    finally:
        lock.unlink()


def main(argv: list[str] | None = None) -> int:
    """Acquire only; the canonical preparation workflow performs validation."""
    from prepare_runtime_databases import (
        DEFAULT_REMOTE_SOURCE_MANIFEST,
        load_remote_source_manifest,
        resolve_remote_version,
    )

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, default=DEFAULT_REMOTE_SOURCE_MANIFEST)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--tools", default="eggnog,cogclassifier,pfam,kofam")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    selected = args.tools.split(",")
    if len(selected) != len(set(selected)) or not set(selected) <= {
        "eggnog",
        "cogclassifier",
        "pfam",
        "kofam",
    }:
        parser.error("Unknown or duplicate resource selection")
    sources = load_remote_source_manifest(args.sources)
    try:
        for tool in selected:
            version, config = resolve_remote_version(
                manifest=sources, component=tool, requested_version=None
            )
            acquire_component(tool, version, config, args.outdir / tool)
            LOGGER.info("Acquired %s; indexing and validation remain required", tool)
    except (ResourceDownloadError, OSError, ValueError) as error:
        LOGGER.error("%s", error)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
