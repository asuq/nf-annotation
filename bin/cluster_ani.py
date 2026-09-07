#!/usr/bin/env python3
"""Cluster ANI-eligible genomes from a FastANI matrix into stable cluster IDs."""

from __future__ import annotations

import argparse
import csv
import importlib
import logging
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, NoReturn, Sequence

from ani_common import (
    AniInputError,
    build_genome_from_row,
    is_empty_ani_cohort,
    load_ani_metadata,
    load_cluster_metadata,
    load_matrix,
    normalize_header,
    subset_matrix,
    try_parse_busco,
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Cluster ANI-eligible genomes from a PHYLIP lower-triangular FastANI matrix "
            "using complete linkage."
        )
    )
    parser.add_argument(
        "-a",
        "--ani-matrix",
        required=True,
        type=Path,
        help="PHYLIP lower-triangular ANI matrix from FastANI --matrix.",
    )
    parser.add_argument(
        "-l",
        "--ani-metadata",
        required=True,
        type=Path,
        help="ANI metadata TSV emitted by build_fastani_inputs.py.",
    )
    parser.add_argument(
        "--matrix-name-column",
        default="matrix_name",
        help="Column in --ani-metadata that matches PHYLIP row names.",
    )
    parser.add_argument(
        "-t",
        "--threshold",
        type=float,
        default=0.95,
        help="ANI threshold as a fraction in (0,1). Default: 0.95.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        type=Path,
        help="Output directory for cluster.tsv.",
    )
    parser.add_argument(
        "-p",
        "--threads",
        type=int,
        default=8,
        help="Upper bound on concurrency for numeric libraries. Default: 8.",
    )
    log_group = parser.add_mutually_exclusive_group()
    log_group.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Show WARNING and ERROR messages only.",
    )
    log_group.add_argument(
        "-d",
        "--debug",
        action="store_true",
        help="Enable verbose debug logging.",
    )
    return parser.parse_args(argv)


def configure_logging(quiet: bool, debug: bool) -> None:
    """Configure logging for the clustering CLI."""
    if quiet:
        level = logging.WARNING
    elif debug:
        level = logging.DEBUG
    else:
        level = logging.INFO

    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
        force=True,
    )


def die(message: str) -> NoReturn:
    """Log one fatal message and terminate."""
    logging.critical(message)
    sys.exit(1)


def _available_cpus() -> int:
    """Return the effective CPU count available to this process."""
    try:
        return len(os.sched_getaffinity(0))  # type: ignore[attr-defined]
    except Exception:
        return os.cpu_count() or 1


def _cgroup_cpu_cap() -> int | None:
    """Return the cgroup CPU cap when detectable."""
    try:
        cpu_max = Path("/sys/fs/cgroup/cpu.max")
        if cpu_max.is_file():
            parts = cpu_max.read_text(encoding="utf-8").strip().split()
            if len(parts) >= 2 and parts[0] != "max":
                quota = int(parts[0])
                period = int(parts[1])
                if quota > 0 and period > 0:
                    return max(1, quota // period)
    except Exception:
        pass

    try:
        quota_file = Path("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")
        period_file = Path("/sys/fs/cgroup/cpu/cpu.cfs_period_us")
        if quota_file.is_file() and period_file.is_file():
            quota = int(quota_file.read_text(encoding="utf-8").strip())
            period = int(period_file.read_text(encoding="utf-8").strip())
            if quota > 0 and period > 0:
                return max(1, quota // period)
    except Exception:
        pass
    return None


def resolve_thread_cap(requested: int) -> int:
    """Cap threads against hardware, scheduler, and cgroup limits."""
    candidates: list[int] = [max(1, requested), _available_cpus()]
    for variable in (
        "NSLOTS",
        "PBS_NUM_PPN",
        "PBS_NP",
        "SLURM_CPUS_ON_NODE",
        "SLURM_CPUS_PER_TASK",
    ):
        value = os.environ.get(variable)
        if value and value.isdigit():
            candidates.append(int(value))
    cgroup_cap = _cgroup_cpu_cap()
    if cgroup_cap is not None:
        candidates.append(cgroup_cap)
    return max(1, min(candidates))


def set_thread_envs(n_threads: int) -> None:
    """Cap math-library thread counts before importing numpy and scipy."""

    def cap(existing: str | None, limit: int) -> str:
        """Return the smaller positive integer thread cap."""
        try:
            if existing is not None:
                parsed = int(existing)
                if parsed > 0:
                    return str(min(parsed, limit))
        except Exception:
            pass
        return str(limit)

    for key in (
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "BLIS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "NUMEXPR_MAX_THREADS",
    ):
        os.environ[key] = cap(os.environ.get(key), n_threads)


def load_scipy_complete_linkage() -> tuple[Any, Any, Any]:
    """Load the SciPy helpers required for ANI complete-linkage clustering."""
    try:
        hierarchy = importlib.import_module("scipy.cluster.hierarchy")
        distance = importlib.import_module("scipy.spatial.distance")
    except ModuleNotFoundError as error:
        raise AniInputError(
            "SciPy is required for ANI clustering. Build params.python_container "
            "with scipy and numpy before running the pipeline."
        ) from error
    return hierarchy.fcluster, hierarchy.linkage, distance.squareform


def cluster_complete_linkage(
    ani: "np.ndarray",
    names: Sequence[str],
    threshold: float,
) -> dict[int, list[int]]:
    """Cluster an ANI matrix by complete linkage and enforce the threshold post-check."""
    if len(names) == 1:
        logging.info("Single ANI sample detected; emitting one singleton cluster.")
        return {1: [0]}

    import numpy as np

    fcluster, linkage, squareform = load_scipy_complete_linkage()

    distance = np.where(np.isnan(ani), 1.0, 1.0 - (ani / 100.0))
    np.fill_diagonal(distance, 0.0)
    condensed = squareform(distance, checks=False)

    cut_distance = 1.0 - threshold
    linkage_matrix = linkage(condensed, method="complete")
    labels = fcluster(linkage_matrix, t=cut_distance, criterion="distance")

    clusters: defaultdict[int, list[int]] = defaultdict(list)
    for index, label in enumerate(labels):
        clusters[int(label)].append(index)

    for label, idxs in clusters.items():
        for left_offset, left_index in enumerate(idxs):
            for right_index in idxs[left_offset + 1 :]:
                value = ani[left_index, right_index]
                if np.isnan(value) or value < (threshold * 100.0):
                    raise AniInputError(
                        "Post-check failed for complete-linkage cluster %s: pair "
                        "(%s, %s) has ANI=%s below %.2f."
                        % (
                            label,
                            names[left_index],
                            names[right_index],
                            value,
                            threshold * 100.0,
                        )
                    )
    return dict(clusters)


def assign_cluster_ids(
    clusters: dict[int, list[int]],
    names: Sequence[str],
    accession_by_matrix_name: dict[str, str],
) -> dict[str, list[int]]:
    """Assign stable cluster IDs ordered by size then smallest accession."""
    ordered_clusters = sorted(
        clusters.values(),
        key=lambda idxs: (
            -len(idxs),
            min(accession_by_matrix_name[names[index]] for index in idxs),
        ),
    )
    return {
        f"C{position:06d}": idxs
        for position, idxs in enumerate(ordered_clusters, start=1)
    }


def build_cluster_rows(
    cluster_ids: dict[str, list[int]],
    names: Sequence[str],
    accession_by_matrix_name: dict[str, str],
) -> list[dict[str, str]]:
    """Build the final cluster membership rows."""
    rows: list[dict[str, str]] = []
    for cluster_id in sorted(cluster_ids):
        idxs = cluster_ids[cluster_id]
        for index in sorted(idxs, key=lambda idx: accession_by_matrix_name[names[idx]]):
            matrix_name = names[index]
            rows.append(
                {
                    "Accession": accession_by_matrix_name[matrix_name],
                    "Cluster_ID": cluster_id,
                    "Matrix_Name": matrix_name,
                }
            )
    return rows


def write_clusters(path: Path, rows: Sequence[dict[str, str]]) -> None:
    """Write one stable ANI cluster membership table."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["Accession", "Cluster_ID", "Matrix_Name"],
            delimiter="\t",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_pipeline(args: argparse.Namespace) -> None:
    """Run ANI complete-linkage clustering and write cluster memberships."""
    if not (0.0 < float(args.threshold) < 1.0):
        raise AniInputError(f"--threshold must be a fraction in (0,1). Got: {args.threshold}")

    if is_empty_ani_cohort(
        args.ani_metadata, args.ani_matrix, matrix_name_column=args.matrix_name_column
    ):
        write_clusters(args.outdir / "cluster.tsv", [])
        logging.info("No ANI-eligible samples; wrote an empty cluster table.")
        return

    names, ani, name_to_idx = load_matrix(args.ani_matrix)
    accession_by_matrix_name, eligible_names = load_cluster_metadata(
        args.ani_metadata,
        names,
        matrix_name_column=args.matrix_name_column,
    )
    names, ani, name_to_idx = subset_matrix(names, ani, eligible_names, name_to_idx)
    del name_to_idx

    clusters = cluster_complete_linkage(ani, names, float(args.threshold))
    stable_clusters = assign_cluster_ids(clusters, names, accession_by_matrix_name)
    rows = build_cluster_rows(stable_clusters, names, accession_by_matrix_name)
    write_clusters(args.outdir / "cluster.tsv", rows)
    logging.info("Wrote ANI cluster memberships for %d sample(s).", len(rows))


def main(argv: Sequence[str] | None = None) -> int:
    """Run the ANI clustering CLI."""
    args = parse_args(argv)
    configure_logging(args.quiet, args.debug)

    requested_threads = max(1, int(args.threads))
    resolved_threads = resolve_thread_cap(requested_threads)
    if resolved_threads < requested_threads:
        logging.info(
            "Capping --threads from %d to %d based on hardware/scheduler limits.",
            requested_threads,
            resolved_threads,
        )
    set_thread_envs(resolved_threads)

    try:
        run_pipeline(args)
    except AniInputError as error:
        logging.error(str(error))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
