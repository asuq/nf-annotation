"""Synthetic, exact-output comparison of ANI post-processing implementations."""

import argparse
import json
import subprocess
import sys
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--baseline", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
parser.add_argument("--samples", type=int, default=2500)
args = parser.parse_args()
candidate = Path(__file__).resolve().parents[1]
baseline = args.baseline
root = args.output
root.mkdir(parents=True, exist_ok=True)
size = args.samples
path = root / "matrix"
with path.open("w") as handle:
    handle.write(str(size) + "\n")
    for i in range(size):
        handle.write(
            f"SYN{i:06}"
            + "".join("\t" + str(96 + (i + j) % 4) for j in range(i))
            + "\n"
        )
script = """import sys,time,json,resource
from pathlib import Path
sys.path.insert(0,sys.argv[1]+"/bin")
from ani_common import Genome,load_matrix
from cluster_ani import cluster_complete_linkage,assign_cluster_ids,build_cluster_rows
from ani_representatives import select_representative_for_indices
start=time.monotonic()
names,ani,index=load_matrix(Path(sys.argv[2]))
clusters=cluster_complete_linkage(ani,names,.95)
rows=build_cluster_rows(assign_cluster_ids(clusters,names,{n:n for n in names}),names,{n:n for n in names})
if sys.argv[3]=="old": names,ani,index=load_matrix(Path(sys.argv[2]))
meta={n:Genome(n,n,11,95.,1.,100000,30,3000000,"C:98.0%,M:1.0%",98.,1.,"Complete Genome",3,n) for n in names}
result=select_representative_for_indices(list(range(len(names))),names,meta,ani,"default")
import hashlib
identity=hashlib.sha256(json.dumps([rows,result],sort_keys=True).encode()).hexdigest()
print(json.dumps(dict(seconds=time.monotonic()-start,maxrss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,output_sha256=identity)))
"""
runs = []
for source, mode in (
    (baseline, "old"),
    (candidate, "new"),
    (candidate, "new"),
    (baseline, "old"),
):
    result = subprocess.run(
        [sys.executable, "-c", script, str(source), str(path), mode],
        capture_output=True,
        text=True,
        check=True,
    )
    record = json.loads(result.stdout)
    record["mode"] = mode
    runs.append(record)
assert len({r["output_sha256"] for r in runs}) == 1
report = {
    "samples": size,
    "synthetic": True,
    "case": "one complete-linkage cluster with tied quality scores and nonconstant ANI centrality",
    "runs": runs,
    "outputs_identical": True,
}
(root / "benchmark.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
