#!/usr/bin/env python3
"""Cluster and select representatives from one parsed ANI matrix."""

import sys

from ani_common import AniInputError
from cluster_ani import (
    configure_logging,
    parse_args,
    run_pipeline,
    set_thread_envs,
    resolve_thread_cap,
)
from select_ani_representatives import (
    ANI_REPRESENTATIVE_COLUMNS,
    ANI_SUMMARY_COLUMNS,
    RepresentativeSelectionError,
    build_ani_outputs,
    write_tsv,
)


def main():
    args = parse_args()
    configure_logging(args.quiet, args.debug)
    threads = resolve_thread_cap(max(1, args.threads))
    set_thread_envs(threads)
    try:
        args.outdir.mkdir(parents=True, exist_ok=True)
        matrix = run_pipeline(args)
        summaries, representatives = build_ani_outputs(
            ani_clusters=args.outdir / "cluster.tsv",
            ani_metadata=args.ani_metadata,
            ani_matrix=args.ani_matrix,
            score_profile=args.ani_score_profile,
            matrix_data=matrix,
            threads=threads,
        )
        write_tsv(args.outdir / "ani_summary.tsv", ANI_SUMMARY_COLUMNS, summaries)
        write_tsv(
            args.outdir / "ani_representatives.tsv",
            ANI_REPRESENTATIVE_COLUMNS,
            representatives,
        )
    except (AniInputError, RepresentativeSelectionError, OSError) as error:
        sys.exit(str(error))


if __name__ == "__main__":
    main()
