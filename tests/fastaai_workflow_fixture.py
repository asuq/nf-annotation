"""Explicit external-tool substitutes for Nextflow stub runs, never native qualification."""

import argparse
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from annotation_common import read_json, write_json


def prepare(request):
    import prepare_aai

    runtime = Path("stub_runtime.json")
    write_json(
        runtime,
        {
            "upstream_commit": "stub",
            "hmm_sha256": "stub",
            "pyhmmer": "stub",
            "numpy": "stub",
            "environment_sha256": "stub",
        },
    )
    prepare_aai.RUNTIME = runtime

    def extract(job):
        _, out = job
        out.mkdir(exist_ok=True)
        (out / "markers.faa").write_text(">stub_protein\nAAAA\n")
        (out / "markers.hmm").write_text(
            "stub_protein - stub PF00000.1 1e-10 100 0 1e-10 100\n"
        )
        return 1, out

    class Pool:
        def __init__(self, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def map(self, unused, jobs):
            return map(extract, jobs)

    for key in ("records", "bundles", "sources"):
        path = Path(key + ".json")
        write_json(path, request[key])
        request[key] = str(path)
    args = SimpleNamespace(
        **{
            k: Path(v) if k not in ("threads", "primary_busco_column", "domain") else v
            for k, v in request.items()
        }
    )
    with patch.object(prepare_aai.concurrent.futures, "ProcessPoolExecutor", Pool):
        prepare_aai.prepare(args)


def compare(records, domain):
    output = Path("comparison")
    output.mkdir()
    ready = sorted(
        (r for r in records if r["status"] == "ready"), key=lambda r: r["accession"]
    )
    from aai_evidence import create

    arrays = create(output, [r["accession"] for r in ready])
    arrays["jaccard"][:] = 0.55
    arrays["shared"][:] = 1
    arrays["possible"][:] = 1
    for array in arrays.values():
        if hasattr(array, "flush"):
            array.flush()
    write_json(output / "records.json", {"domain": domain, "samples": records})
    write_json(
        output / "comparison.json",
        {"runtime": "stub", "accessions": [r["accession"] for r in ready]},
    )
    (output / "native.log").write_text(
        "Synthetic workflow control; FastAAI was not executed.\n"
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("operation", choices=("prepare", "compare"))
    p.add_argument("input", type=Path)
    p.add_argument("--domain")
    args = p.parse_args()
    if args.operation == "prepare":
        prepare(read_json(args.input))
    else:
        compare(read_json(args.input), args.domain)
