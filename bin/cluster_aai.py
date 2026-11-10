#!/usr/bin/env python3
"""Complete-linkage genus-scale clusters and quality-only representatives."""

from __future__ import annotations
import argparse
import time
from pathlib import Path
import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage

from aai_common import (
    AAI_COLUMNS,
    REP_COLUMNS,
    SUMMARY_COLUMNS,
    aai_estimate,
    jaccard_cutoff,
)
from ani_common import build_genome_from_row
from build_fastani_inputs import load_busco_index, load_index, load_metadata_index
from build_master_table import resolve_assembly_metric_value
from ani_representatives import select_representative_for_indices
from annotation_common import AnnotationError, read_json, write_json, write_tsv
from annotation_path_lists import read_path_list


def offset(n, i, j):
    if i > j:
        i, j = j, i
    return n * i - i * (i + 1) // 2 + j - i - 1


def cluster_domain(records, pairs, threshold, score_profile, domain):
    records = sorted(records, key=lambda r: r["accession"])
    rows = {
        r["accession"]: dict(
            Accession=r["accession"],
            **{c: "NA" for c in AAI_COLUMNS},
            aai_status=r["status"],
            aai_exclusion_reason=r["reason"],
            aai_representative_status="not_applicable",
        )
        for r in records
    }
    ready = [r for r in records if r["status"] == "ready"]
    names = [r["accession"] for r in ready]
    n = len(names)
    indices = {a: i for i, a in enumerate(names)}
    from aai_evidence import load

    evidence = load(pairs, names)
    jaccard = evidence["jaccard"]
    shared = evidence["shared"]
    valid = shared > 0
    if np.any(
        valid & (~np.isfinite(jaccard) | (jaccard < 0) | (jaccard > 1))
    ) or np.any(shared > evidence["possible"]):
        raise AnnotationError("Invalid AAI pair evidence")
    distances = np.subtract(1.0, jaccard)
    distances[~valid] = 2.0
    comparable = np.zeros(n, dtype=bool)
    for i in range(n - 1):
        start = offset(n, i, i + 1)
        row = valid[start : start + n - i - 1]
        comparable[i] |= row.any()
        comparable[i + 1 :] |= row
    if n == 1:
        comparable[0] = True
    for i, name in enumerate(names):
        if not comparable[i]:
            rows[name].update(
                aai_status="no_pairwise_evidence",
                aai_exclusion_reason="no_shared_markers_with_cohort",
            )
    labels = (
        fcluster(
            linkage(distances, method="complete"),
            t=1 - jaccard_cutoff(threshold),
            criterion="distance",
        )
        if n > 1
        else np.ones(n, dtype=int)
    )
    groups = {}
    for i, label in enumerate(labels):
        if comparable[i]:
            groups.setdefault(int(label), []).append(i)
    groups = sorted(
        groups.values(), key=lambda ids: (-len(ids), min(names[i] for i in ids))
    )
    representatives = []
    members = []
    traces = []
    for serial, ids in enumerate(groups, 1):
        cluster = f"AAI_{domain[0]}{serial:06d}"
        candidates = {}
        missing = {}
        for i in ids:
            name = names[i]
            genome, reasons = build_genome_from_row(
                ready[i]["quality"],
                acc=name,
                busco_column="busco",
                require_existing_path=False,
                require_genome_size=False,
            )
            if genome is not None:
                candidates[name] = genome
            else:
                missing[name] = ";".join(reasons)
        chosen = None
        if candidates:
            candidate_names = sorted(candidates)
            _, chosen, trace = select_representative_for_indices(
                list(range(len(candidate_names))),
                candidate_names,
                candidates,
                None,
                score_profile,
            )
            traces.extend([dict(cluster=cluster, line=line) for line in trace])
        representatives.append(
            dict(
                AAI_Cluster_ID=cluster,
                Domain=domain,
                Cluster_Size=len(ids),
                Representative_Accession=chosen or "NA",
                Representative_Status="selected"
                if chosen
                else "no_eligible_representative",
                AAI_Score=f"{candidates[chosen].Score:.6f}" if chosen else "NA",
            )
        )
        for i in ids:
            name = names[i]
            row = rows[name]
            row.update(
                AAI_Cluster_ID=cluster,
                aai_status="clustered",
                aai_exclusion_reason="NA",
                AAI_Is_Representative=("yes" if name == chosen else "no")
                if chosen
                else "NA",
                AAI_Score=f"{candidates[name].Score:.6f}"
                if name in candidates
                else "NA",
                aai_representative_status=("selected" if name == chosen else "eligible")
                if name in candidates
                else "missing_scoring_metrics:" + missing[name],
            )
            if chosen:
                row["AAI_to_Representative"] = (
                    "100"
                    if name == chosen
                    else aai_estimate(float(jaccard[offset(n, i, indices[chosen])]))
                )
            members.append(dict(Accession=name, AAI_Cluster_ID=cluster, Domain=domain))
    return list(rows.values()), representatives, members, traces


def run(args):
    start = time.monotonic()
    all_rows = []
    all_reps = []
    members = []
    traces = []
    domains = []
    _, metadata = load_metadata_index(args.metadata)
    stats = load_index(args.stats, ("accession",))
    busco = load_busco_index(read_path_list(args.busco))
    for path in read_path_list(args.domains):
        data = read_json(path / "records.json")
        domain = data["domain"]
        if domain not in ("Bacteria", "Archaea") or domain in domains:
            raise AnnotationError("Invalid or duplicate AAI domain")
        domains.append(domain)
        if any(r["domain"] != domain for r in data["samples"]):
            raise AnnotationError("Mixed AAI domains")
        for row in data["samples"]:
            accession = row["accession"]
            for key, column in (
                ("n50", "N50"),
                ("scaffolds", "Scaffolds"),
                ("genome_size", "Genome_Size"),
            ):
                row["quality"][key] = resolve_assembly_metric_value(
                    metadata.get(accession, {}), stats.get(accession), column
                )
            row["quality"]["busco"] = busco.get(accession, {}).get(
                row["primary_busco_column"], "NA"
            )
        rows, reps, group, trace = cluster_domain(
            data["samples"],
            path,
            args.threshold,
            args.score_profile,
            domain,
        )
        all_rows.extend(rows)
        all_reps.extend(reps)
        members.extend(group)
        traces.extend(trace)
    expected = [
        r["accession"] for r in __import__("annotation_common").read_tsv(args.samples)
    ]
    index = {r["Accession"]: r for r in all_rows}
    if len(index) != len(all_rows) or set(index) != set(expected):
        raise AnnotationError("AAI output membership differs from requested samples")
    all_reps.sort(key=lambda r: r["AAI_Cluster_ID"])
    members.sort(key=lambda r: (r["AAI_Cluster_ID"], r["Accession"]))
    traces.sort(key=lambda r: r["cluster"])
    args.output.mkdir()
    write_tsv(
        args.output / "aai_summary.tsv", SUMMARY_COLUMNS, [index[a] for a in expected]
    )
    write_tsv(args.output / "aai_representatives.tsv", REP_COLUMNS, all_reps)
    write_tsv(
        args.output / "cluster.tsv", ("Accession", "AAI_Cluster_ID", "Domain"), members
    )
    write_json(
        args.output / "method.json",
        dict(
            threshold=args.threshold,
            score_profile=args.score_profile,
            linkage="complete",
            distance="one_minus_native_mean_jaccard",
            jaccard_cutoff=jaccard_cutoff(args.threshold),
            centrality=False,
        ),
    )
    write_json(args.output / "scores.json", traces)
    write_json(
        args.output / "metrics.json",
        dict(seconds=time.monotonic() - start, samples=len(expected)),
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("domains", "samples", "metadata", "stats", "busco", "output"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--threshold", type=float, default=0.65)
    p.add_argument(
        "--score-profile", choices=("default", "isolate", "mag"), default="default"
    )
    args = p.parse_args()
    try:
        run(args)
    except (AnnotationError, OSError, ValueError, KeyError) as e:
        p.exit(1, str(e) + "\n")


if __name__ == "__main__":
    main()
