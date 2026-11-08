#!/usr/bin/env python3
"""Prepare MD5 receipts for local annotation SIF images, once after staging."""

import argparse
import hashlib
import json
import re
from pathlib import Path


def receipt_path(image: Path) -> Path:
    return image.with_name(image.name + ".receipt.json")


def metadata(image: Path) -> dict:
    stat = image.stat()
    if not image.is_file() or stat.st_size == 0:
        raise ValueError(f"Missing or empty SIF image: {image}")
    return {"bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def read_receipt(image: Path) -> dict:
    try:
        record = json.loads(receipt_path(image).read_text())
        if (
            not isinstance(record, dict)
            or record.get("schema") != "annotation-sif-md5-v1"
            or not re.fullmatch(r"[0-9a-f]{32}", record.get("md5", ""))
            or record.get("file") != metadata(image)
        ):
            raise ValueError("changed image or malformed receipt")
        return record
    except (OSError, ValueError, TypeError) as error:
        raise ValueError(
            f"Prepare the SIF receipt with prepare_annotation_runtime.py: {image}: {error}"
        ) from error


def prepare(image: Path) -> dict:
    image = image.absolute()
    before = metadata(image)
    with image.open("rb") as stream:
        checksum = hashlib.file_digest(stream, "md5").hexdigest()
    if metadata(image) != before:
        raise ValueError(f"Image changed during preparation: {image}")
    record = {"schema": "annotation-sif-md5-v1", "file": before, "md5": checksum}
    target = receipt_path(image)
    temporary = target.with_name(target.name + ".tmp")
    temporary.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    temporary.replace(target)
    return record


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("images", nargs="+", type=Path)
    args = parser.parse_args()
    for image in args.images:
        prepare(image)
