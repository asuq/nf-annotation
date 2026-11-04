#!/usr/bin/env python3
"""Representative-selection helpers for ANI clusters.

This module scores all members within one ANI cluster, selects exactly one
representative with a deterministic tie cascade, and formats the per-sample and
per-cluster ANI outputs used downstream.
"""

from __future__ import annotations

import logging
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import TYPE_CHECKING, Sequence

if TYPE_CHECKING:
    import numpy as np

from ani_common import Genome


def select_representative_for_indices(
    idxs: list[int],
    names: Sequence[str],
    meta: dict[str, Genome],
    ani: "np.ndarray",
    score_profile: str,
) -> tuple[int, str, list[str]]:
    """Select one representative for one ANI cluster.

    Each cluster member is scored from six components:

    - assembly rank
    - CheckM2 quality
    - BUSCO completeness-minus-missing
    - N50
    - scaffold count
    - ANI centrality within the cluster

    Raw BUSCO, CheckM2, N50, and scaffold values are winsorized for larger
    clusters, normalized to the [0, 1] range within the cluster, and combined
    with the chosen weight profile. If the top score is tied within a small
    epsilon, the deterministic tie cascade is:

    1. higher assembly rank
    2. higher BUSCO complete percentage
    3. lower BUSCO missing percentage
    4. lower contamination
    5. higher completeness
    6. fewer scaffolds
    7. higher N50
    8. lexicographically smallest accession

    The return value contains the cluster size, the chosen representative
    matrix-name token, and the debug trace lines explaining the decision.
    """
    import numpy as np

    epsilon = 1e-6
    debug_lines: list[str] = []
    matrix_names = [names[index] for index in idxs]
    genomes = [meta[matrix_name] for matrix_name in matrix_names]
    cluster_size = len(genomes)

    def winsorize(values: list[float]) -> "np.ndarray":
        """Winsorize a value vector within one cluster."""
        if cluster_size < 8:
            return np.asarray(values, dtype=float)
        array = np.asarray(values, dtype=float)
        low, high = np.quantile(array, [0.05, 0.95])
        return np.clip(array, low, high)

    def minmax_norm(values: "np.ndarray") -> "np.ndarray":
        """Normalize a vector into [0,1] with a constant fallback."""
        low = float(values.min())
        high = float(values.max())
        if high > low + 1e-12:
            return (values - low) / (high - low)
        return np.full_like(values, 0.5, dtype=float)

    def ani_centrality(cluster_indices: list[int]) -> "np.ndarray":
        """Compute ANI centrality within one cluster."""
        member_count = len(cluster_indices)
        if member_count <= 1:
            return np.ones(member_count, dtype=float)

        means: list[float] = []
        for index in cluster_indices:
            others = ani[index, [other for other in cluster_indices if other != index]]
            means.append(float(np.mean(others)) if len(others) else 100.0)

        values = np.asarray(means, dtype=float)
        spread = float(values.max() - values.min())
        if spread < 0.05:
            return np.full_like(values, 0.5, dtype=float)
        return (values - float(values.min())) / spread

    assembly_values = np.array([genome.Assembly_Rank / 3.0 for genome in genomes], dtype=float)
    busco_raw = np.array(
        [(genome.BUSCO_C - genome.BUSCO_M) / 100.0 for genome in genomes],
        dtype=float,
    )
    busco_values = minmax_norm(winsorize(busco_raw))

    checkm2_raw = np.array(
        [
            (genome.CheckM2_Completeness - 5.0 * genome.CheckM2_Contamination) / 100.0
            for genome in genomes
        ],
        dtype=float,
    )
    checkm2_values = minmax_norm(winsorize(checkm2_raw))

    n50_values = minmax_norm(
        winsorize(np.log10(np.array([genome.N50 for genome in genomes], dtype=float) + 1.0))
    )
    scaffold_values = 1.0 - minmax_norm(
        winsorize(
            np.log10(np.array([genome.Scaffolds for genome in genomes], dtype=float) + 1.0)
        )
    )
    centrality_values = ani_centrality(idxs)

    if score_profile == "isolate":
        weights = (3.5, 1.5, 1.5, 1.0, 0.25, 0.5)
    elif score_profile == "mag":
        weights = (1.0, 3.0, 2.0, 0.5, 1.5, 0.5)
    else:
        weights = (3.0, 2.0, 1.5, 0.75, 0.5, 0.5)
    weight_assembly, weight_checkm2, weight_busco, weight_n50, weight_scaffolds, weight_centrality = weights

    scored: list[tuple[str, str, float, dict[str, float], Genome]] = []
    for offset, (matrix_name, genome) in enumerate(zip(matrix_names, genomes, strict=True)):
        components = {
            "A": float(assembly_values[offset]),
            "Q": float(checkm2_values[offset]),
            "B": float(busco_values[offset]),
            "N": float(n50_values[offset]),
            "S": float(scaffold_values[offset]),
            "C": float(centrality_values[offset]),
            "Q_raw": float(checkm2_raw[offset]),
            "B_raw": float(busco_raw[offset]),
        }
        score = (
            weight_assembly * components["A"]
            + weight_checkm2 * components["Q"]
            + weight_busco * components["B"]
            + weight_n50 * components["N"]
            + weight_scaffolds * components["S"]
            + weight_centrality * components["C"]
        )
        genome.Score = float(score)
        scored.append((matrix_name, genome.Accession, float(score), components, genome))

    debug_lines.append(
        "Cluster size=%d; profile=%s; weights: A=%s, Q=%s, B=%s, N=%s, S=%s, C=%s"
        % (cluster_size, score_profile, *weights)
    )
    for _matrix_name, accession, score, components, genome in sorted(
        scored,
        key=lambda item: (-item[2], item[1]),
    ):
        debug_lines.append(
            "  %s: score=%.6f | A=%.3f Q=%.3f(q_raw=%.3f) B=%.3f(b_raw=%.3f) "
            "N=%.3f S=%.3f C=%.3f | asm=%s(rank=%d) busco(C=%.2f,M=%.2f) "
            "checkm2=%.2f contam=%.2f scaffolds=%d n50=%d"
            % (
                accession,
                score,
                components["A"],
                components["Q"],
                components["Q_raw"],
                components["B"],
                components["B_raw"],
                components["N"],
                components["S"],
                components["C"],
                genome.Assembly_Level,
                genome.Assembly_Rank,
                genome.BUSCO_C,
                genome.BUSCO_M,
                genome.CheckM2_Completeness,
                genome.CheckM2_Contamination,
                genome.Scaffolds,
                genome.N50,
            )
        )

    scored.sort(key=lambda item: (-item[2], item[1]))
    top_score = scored[0][2]
    tied = [
        matrix_name
        for matrix_name, _accession, score, _components, _genome in scored
        if abs(score - top_score) <= epsilon
    ]
    if len(tied) == 1:
        debug_lines.append(f"Chosen by score alone: {meta[tied[0]].Accession} (score={top_score:.6f})")
        return cluster_size, tied[0], debug_lines

    matrix_name_to_genome = {
        matrix_name: genome for matrix_name, _accession, _score, _components, genome in scored
    }

    def reduce_tie(
        matrix_name_list: list[str],
        key_func,
        description: str,
        prefer_lower: bool = False,
    ) -> list[str]:
        """Reduce a tied representative set by one deterministic criterion."""
        if len(matrix_name_list) <= 1:
            return matrix_name_list
        values = [key_func(matrix_name_to_genome[matrix_name]) for matrix_name in matrix_name_list]
        best = min(values) if prefer_lower else max(values)
        winners = [
            matrix_name
            for matrix_name, value in zip(matrix_name_list, values, strict=True)
            if value == best
        ]
        if len(winners) < len(matrix_name_list):
            debug_lines.append(
                f"Tie-break by {description}: {len(matrix_name_list)} -> {len(winners)} (best={best})."
            )
        return winners

    tied = reduce_tie(tied, lambda genome: genome.Assembly_Rank, "Assembly rank (higher)")
    tied = reduce_tie(tied, lambda genome: genome.BUSCO_C, "BUSCO C (higher)")
    tied = reduce_tie(tied, lambda genome: genome.BUSCO_M, "BUSCO M (lower)", True)
    tied = reduce_tie(
        tied,
        lambda genome: genome.CheckM2_Contamination,
        "Contamination (lower)",
        True,
    )
    tied = reduce_tie(tied, lambda genome: genome.CheckM2_Completeness, "CheckM2 (higher)")
    tied = reduce_tie(tied, lambda genome: genome.Scaffolds, "Scaffolds (fewer)", True)
    tied = reduce_tie(tied, lambda genome: genome.N50, "N50 (higher)")

    representative = sorted(
        tied,
        key=lambda matrix_name: matrix_name_to_genome[matrix_name].Accession,
    )[0]
    if len(tied) > 1:
        debug_lines.append(
            "Final tie among %d candidate(s); picking lexicographically smallest: %s"
            % (len(tied), matrix_name_to_genome[representative].Accession)
        )
    else:
        debug_lines.append(f"Chosen by tie cascade: {matrix_name_to_genome[representative].Accession}")
    return cluster_size, representative, debug_lines


def select_representatives_for_clusters(
    cluster_members: dict[str, list[str]],
    names: Sequence[str],
    meta: dict[str, Genome],
    ani: "np.ndarray",
    score_profile: str,
    threads: int = 1,
) -> dict[str, tuple[str, list[str]]]:
    """Select one representative per cluster and return debug traces."""
    name_to_idx = {name: index for index, name in enumerate(names)}
    results: dict[str, tuple[str, list[str]]] = {}

    def run_one(cluster_id: str, matrix_names: list[str]) -> tuple[str, tuple[str, list[str]]]:
        """Select one representative for a single cluster."""
        idxs = [name_to_idx[matrix_name] for matrix_name in matrix_names]
        _cluster_size, representative, debug_lines = select_representative_for_indices(
            idxs,
            names,
            meta,
            ani,
            score_profile,
        )
        return cluster_id, (representative, debug_lines)

    worker_count = max(1, threads)
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        future_to_cluster = {
            executor.submit(run_one, cluster_id, matrix_names): cluster_id
            for cluster_id, matrix_names in cluster_members.items()
        }
        for future in as_completed(future_to_cluster):
            cluster_id = future_to_cluster[future]
            try:
                resolved_cluster_id, payload = future.result()
            except Exception as error:
                raise RuntimeError(
                    f"Representative selection failed for cluster {cluster_id}"
                ) from error
            results[resolved_cluster_id] = payload
    return results


def warn_gcode_mixture(
    cluster_id: str,
    matrix_names: Sequence[str],
    meta: dict[str, Genome],
) -> None:
    """Log a warning when one ANI cluster mixes gcode 4 and 11 genomes."""
    gcode4 = 0
    gcode11 = 0
    for matrix_name in matrix_names:
        gcode = meta[matrix_name].Gcode
        if gcode == 4:
            gcode4 += 1
        elif gcode == 11:
            gcode11 += 1
    if gcode4 > 0 and gcode11 > 0:
        logging.warning(
            "Gcode mixture in %s: gcode4=%d, gcode11=%d",
            cluster_id,
            gcode4,
            gcode11,
        )


def warn_soft_screen_offenders(
    cluster_id: str,
    matrix_names: Sequence[str],
    meta: dict[str, Genome],
) -> None:
    """Log soft-screen ANI warnings for low-quality cluster members."""
    offenders: list[str] = []
    for matrix_name in matrix_names:
        genome = meta[matrix_name]
        reasons: list[str] = []
        if genome.CheckM2_Completeness < 90:
            reasons.append(f"Completeness={genome.CheckM2_Completeness:.2f}")
        if genome.CheckM2_Contamination > 5:
            reasons.append(f"Contam={genome.CheckM2_Contamination:.2f}%")
        if reasons:
            offenders.append(f"{genome.Accession}[{';'.join(reasons)}]")
    if offenders:
        logging.warning("Soft-screen offenders in %s: %s", cluster_id, ", ".join(offenders))


def build_ani_summary_from_clusters(
    cluster_members: dict[str, list[str]],
    names: Sequence[str],
    ani: "np.ndarray",
    name_to_idx: dict[str, int],
    meta: dict[str, Genome],
    score_profile: str = "default",
    threads: int = 1,
) -> tuple[dict[str, dict[str, str]], list[dict[str, str]]]:
    """Build ANI summary rows and representative rows from cluster memberships.

    The first return value is an accession-keyed mapping for the master-table
    ANI columns: `Cluster_ID`, `Is_Representative`, `ANI_to_Representative`,
    and `Score`. The second return value is the ordered representative table
    with one row per cluster for the published `ani_representatives.tsv`.
    """
    import numpy as np

    representative_rows: list[dict[str, str]] = []
    summary_index: dict[str, dict[str, str]] = {}
    representative_map = select_representatives_for_clusters(
        cluster_members=cluster_members,
        names=names,
        meta=meta,
        ani=ani,
        score_profile=score_profile,
        threads=threads,
    )

    for cluster_id in sorted(cluster_members):
        matrix_names = cluster_members[cluster_id]
        representative_matrix_name, debug_lines = representative_map[cluster_id]
        for line in debug_lines:
            logging.debug("%s (%s): %s", cluster_id, representative_matrix_name, line)

        warn_gcode_mixture(cluster_id, matrix_names, meta)
        warn_soft_screen_offenders(cluster_id, matrix_names, meta)

        representative = meta[representative_matrix_name]
        representative_rows.append(
            {
                "Cluster_ID": cluster_id,
                "Representative_Accession": representative.Accession,
                "Organism_Name": re.sub(r"\s+", "_", representative.Organism_Name.strip()),
                "CheckM2_Completeness": f"{representative.CheckM2_Completeness:.2f}",
                "CheckM2_Contamination": f"{representative.CheckM2_Contamination:.2f}",
                "BUSCO": representative.BUSCO_str,
                "Assembly_Level": representative.Assembly_Level,
                "N50": str(representative.N50),
                "Cluster_Size": str(len(matrix_names)),
            }
        )

        representative_index = name_to_idx[representative_matrix_name]
        for matrix_name in sorted(matrix_names, key=lambda item: meta[item].Accession):
            genome = meta[matrix_name]
            if genome.Score is None:
                raise RuntimeError(
                    f"Representative score was not set for accession {genome.Accession!r}."
                )
            if matrix_name == representative_matrix_name:
                ani_to_representative = "100.0000"
                is_representative = "yes"
            else:
                value = ani[name_to_idx[matrix_name], representative_index]
                if np.isnan(value):
                    raise RuntimeError(
                        "Internal error: NA ANI inside a completed ANI cluster for "
                        f"{genome.Accession!r} versus {representative.Accession!r}."
                    )
                ani_to_representative = f"{value:.4f}"
                is_representative = "no"

            summary_index[genome.Accession] = {
                "Cluster_ID": cluster_id,
                "Is_Representative": is_representative,
                "ANI_to_Representative": ani_to_representative,
                "Score": f"{genome.Score:.6f}",
            }

    return summary_index, representative_rows
