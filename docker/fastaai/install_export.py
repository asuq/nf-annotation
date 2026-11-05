"""Add binary serialization after native Jaccard calculations; retain original CLI writers."""

import ast
import hashlib
import json
import shutil
from pathlib import Path
import fastaai.fastaai as native
import numpy
import pyhmmer

path = Path(native.__file__)
original = path.read_text()
path.with_suffix(".upstream.py").write_text(original)
lines = original.splitlines(keepends=True)
insertions = []
for function in ast.parse(original).body:
    if isinstance(function, ast.FunctionDef) and function.name in (
        "one_work",
        "on_disk_work_one",
    ):
        blocks = [
            node
            for node in ast.walk(function)
            if isinstance(node, ast.If) and ast.unparse(node.test) == "style == 'tsv'"
        ]
        assert len(blocks) == 1, "Unexpected upstream writer layout"
        line = blocks[0].lineno - 1
        indent = lines[line][: len(lines[line]) - len(lines[line].lstrip())]
        body = (
            indent
            + 'if style == "binary":\n'
            + indent
            + "\tfrom fastaai.binary_export import write_row\n"
            + indent
            + "\twrite_row(out_base, q, jaccard_averages, shared_acc_counts, np.minimum(local_qgak[q][0].shape[0], _tct))\n"
            + indent
            + 'elif style == "tsv":\n'
        )
        insertions.append((line, body))
assert len(insertions) == 2
for line, body in sorted(insertions, reverse=True):
    lines[line] = body
path.write_text("".join(lines))
writer = Path("/opt/nf-annotation/provenance/native_binary_export.py")
shutil.copyfile(writer, path.with_name("binary_export.py"))
models = Path(native.find_hmm())
record = dict(
    upstream_commit="3b4904a223590f8227e9c791b6b7ed6a09c67757",
    native_version="0.1.17",
    environment_sha256=hashlib.sha256(
        Path("/opt/nf-annotation/fastaai/pixi.lock").read_bytes()
    ).hexdigest(),
    numpy=numpy.__version__,
    pyhmmer=pyhmmer.__version__,
    source_sha256=hashlib.sha256(original.encode()).hexdigest(),
    export_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    export_writer_sha256=hashlib.sha256(writer.read_bytes()).hexdigest(),
    hmm_sha256=hashlib.sha256(models.read_bytes()).hexdigest(),
    models=str(models),
    export="native_condensed_v1",
)
Path("/opt/nf-annotation/provenance/fastaai.json").write_text(
    json.dumps(record, sort_keys=True, indent=2) + "\n"
)
