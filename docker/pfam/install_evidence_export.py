#!/usr/bin/env python3
"""Retain HMMER tables while leaving the pinned PfamScan search unchanged."""

import hashlib
import json
from pathlib import Path

path = Path("/opt/nf-annotation/pfam/PfamScan/Bio/Pfam/Scan/PfamScan.pm")
expected = "f83dfa79d8d3c51dca2ed04d052b86f1e8cb2ef85451529605769860bbeea2ae"
source = path.read_text()
if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
    raise ValueError("Unexpected PfamScan source bytes")
anchor = '    print STDERR "PfamScan::search: hmmscan command: |@params|\\n"'
if source.count(anchor) != 1:
    raise ValueError("Expected one PfamScan command anchor")
addition = """    if (defined $ENV{NF_PFAM_EVIDENCE_DIR}) {
      splice @params, 1, 0,
        '--domtblout', "$ENV{NF_PFAM_EVIDENCE_DIR}/hmmscan.domtblout",
        '--tblout', "$ENV{NF_PFAM_EVIDENCE_DIR}/hmmscan.tblout";
    }

"""
path.write_text(source.replace(anchor, addition + anchor, 1))
provenance = Path("/opt/nf-annotation/provenance")
provenance.mkdir(parents=True, exist_ok=True)
(provenance / "PfamScan.pm").write_text(source)
(provenance / "pfam_export.json").write_text(
    json.dumps(
        {
            "upstream_version": "1.6",
            "archive_sha256": "365c96bc150d5057349c3016d62667c58cb33afcfb6329457ae16ab5aae4f401",
            "original_module_sha256": expected,
            "patched_module_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "installer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
        indent=2,
    )
    + "\n"
)
