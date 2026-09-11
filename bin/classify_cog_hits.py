#!/usr/bin/env python3
"""Run pinned COGclassifier's native classification on an offline RPS-BLAST result."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def main() -> int:
    """Export the native dataframe; database acquisition belongs to preparation."""
    from cogclassifier.blast import BlastAlignmentRecord
    from cogclassifier.cog import (
        CogCddIdTable,
        CogClassifyStats,
        CogDefinitionRecord,
        CogFuncCategoryRecord,
    )

    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("input", "hits", "database", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    stats = CogClassifyStats(
        args.input,
        BlastAlignmentRecord(args.hits),
        CogFuncCategoryRecord(args.database / "cog_func_category.tsv"),
        CogDefinitionRecord(args.database / "cog_definition.tsv"),
        CogCddIdTable(args.database / "cddid.tbl"),
    )
    stats.query_classify_df.to_csv(args.output, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
