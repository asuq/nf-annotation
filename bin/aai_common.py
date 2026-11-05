"""FastAAI calibration, evidence and output contracts."""

from __future__ import annotations
import math

UPSTREAM_COMMIT = "3b4904a223590f8227e9c791b6b7ed6a09c67757"
AAI_COLUMNS = (
    "AAI_Cluster_ID",
    "AAI_Is_Representative",
    "AAI_to_Representative",
    "AAI_Score",
)
AAI_STATUS_COLUMNS = ("aai_status", "aai_exclusion_reason", "aai_representative_status")
SUMMARY_COLUMNS = ("Accession", *AAI_COLUMNS, *AAI_STATUS_COLUMNS)
PAIR_COLUMNS = (
    "query",
    "target",
    "avg_jacc_sim",
    "num_shared_SCPs",
    "poss_shared_SCPs",
    "AAI_estimate",
)
REP_COLUMNS = (
    "AAI_Cluster_ID",
    "Domain",
    "Cluster_Size",
    "Representative_Accession",
    "Representative_Status",
    "AAI_Score",
)


def jaccard_cutoff(threshold: float) -> float:
    if not math.isfinite(threshold) or not 0.30 <= threshold <= 0.90:
        raise ValueError("AAI threshold must be a finite fraction from 0.30 to 0.90")
    return math.exp(
        -((-math.log((threshold + 0.3087057) / 1.810741)) ** 3.435) / 0.2607023
    )


def aai_estimate(jaccard: float) -> str:
    if not math.isfinite(jaccard) or not 0 <= jaccard <= 1:
        raise ValueError("Invalid native Jaccard similarity")
    if jaccard == 0:
        return "<30%"
    value = (
        -0.3087057
        + 1.810741 * math.exp(-((-0.2607023 * math.log(jaccard)) ** (1 / 3.435)))
    ) * 100
    return "<30%" if value < 30 else ">90%" if value > 90 else format(value, ".15g")
