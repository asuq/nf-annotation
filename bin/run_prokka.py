#!/usr/bin/env python3
"""Observe released Prokka/Prodigal execution without changing tool arguments."""

from __future__ import annotations

import json
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
from pathlib import Path

PREFIX = "NF_PROKKA_PREDICTION_V1"


def argument(arguments: list[str], flag: str) -> str | None:
    return arguments[arguments.index(flag) + 1] if flag in arguments else None


def observe_child(arguments: list[str]) -> int:
    executable = os.environ["NF_PROKKA_REAL_PRODIGAL"]
    if "-i" not in arguments:
        os.execv(executable, [executable, *arguments])  # Native version probes.
    result = subprocess.run([executable, *arguments], check=False)
    record = {
        "arguments": arguments,
        "mode": argument(arguments, "-p"),
        "genetic_code": argument(arguments, "-g"),
        "returncode": result.returncode,
    }
    with Path(os.environ["NF_PROKKA_CHILD_RECORD"]).open("a") as stream:
        stream.write(json.dumps(record, sort_keys=True) + "\n")
    if result.returncode < 0:
        signum = -result.returncode
        if signum not in (signal.SIGKILL, signal.SIGSTOP):
            signal.signal(signum, signal.SIG_DFL)
        os.kill(os.getpid(), signum)
    return result.returncode


def run(arguments: list[str]) -> int:
    prokka = shutil.which("prokka")
    prodigal = shutil.which("prodigal")
    record = {
        "schema": 1,
        "prokka": prokka,
        "prodigal": prodigal,
        "arguments": arguments,
        "parent_returncode": None,
        "child": None,
        "cds_count": None,
        "status": "failed",
    }
    code = 1
    try:
        if prokka is None or prodigal is None:
            raise RuntimeError("Released Prokka and Prodigal must both be available")
        with tempfile.TemporaryDirectory(
            prefix="prokka-observer-", dir=os.environ.get("TMPDIR")
        ) as directory:
            observer = Path(directory)
            (observer / "prodigal").symlink_to(Path(__file__).resolve())
            child_record = observer / "child.jsonl"
            environment = dict(
                os.environ,
                PATH=str(observer) + os.pathsep + os.environ["PATH"],
                NF_PROKKA_REAL_PRODIGAL=str(Path(prodigal).absolute()),
                NF_PROKKA_CHILD_RECORD=str(child_record),
            )
            counts = []
            with subprocess.Popen(
                [prokka, *arguments],
                env=environment,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            ) as process:
                assert process.stdout is not None
                for line in process.stdout:
                    sys.stdout.buffer.write(line)
                    sys.stdout.buffer.flush()
                    match = re.fullmatch(rb"\[[0-9:]+\] Found ([0-9]+) CDS\s*", line)
                    if match:
                        counts.append(int(match[1]))
                record["parent_returncode"] = process.wait()
            children = (
                [json.loads(line) for line in child_record.read_text().splitlines()]
                if child_record.exists()
                else []
            )
            if len(children) != 1:
                raise RuntimeError("Missing or ambiguous Prodigal prediction execution")
            record["child"] = children[0]
            child = children[0]
            if record["parent_returncode"] != 0 or child["returncode"] != 0:
                code = record["parent_returncode"] or child["returncode"]
                code = 128 - code if code < 0 else code
                raise RuntimeError(
                    f"Prediction failed: Prokka exit={record['parent_returncode']}, Prodigal returncode={child['returncode']}"
                )
            if child["mode"] not in ("single", "meta") or len(counts) != 1:
                raise RuntimeError("Missing or ambiguous native CDS count/mode")
            record.update(status="success", cds_count=counts[0])
            print(f"{PREFIX}\t{child['mode']}\t{counts[0]}", flush=True)
            code = 0
    except (OSError, RuntimeError, ValueError) as error:
        record["reason"] = str(error)
        print(f"PROKKA_PREDICTION_FAILED: {error}", file=sys.stderr)
    Path("prokka_execution.json").write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n"
    )
    return code


if __name__ == "__main__":
    sys.exit(
        observe_child(sys.argv[1:])
        if Path(sys.argv[0]).name == "prodigal"
        else run(sys.argv[1:])
    )
