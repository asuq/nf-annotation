"""Compact, ordered, non-pickle pair evidence for v0.5 FastAAI consumers."""

from pathlib import Path
import numpy as np
from annotation_common import AnnotationError, read_json, write_json

FORMAT = "fastaai-condensed-v1"
ARRAYS = {
    "jaccard": np.dtype("<f8"),
    "shared": np.dtype("u1"),
    "possible": np.dtype("u1"),
}


def offset(n, i, j):
    if i > j:
        i, j = j, i
    return n * i - i * (i + 1) // 2 + j - i - 1


def create(root, names):
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    count = len(names) * (len(names) - 1) // 2
    arrays = {}
    for name, dtype in ARRAYS.items():
        path = root / (name + ".npy")
        if count:
            arrays[name] = np.lib.format.open_memmap(
                path, mode="w+", dtype=dtype, shape=(count,)
            )
            arrays[name][:] = np.nan if name == "jaccard" else 0
        else:
            arrays[name] = np.empty(0, dtype=dtype)
            np.save(path, arrays[name], allow_pickle=False)
    write_json(
        root / "evidence.json", dict(format=FORMAT, accessions=list(names), pairs=count)
    )
    return arrays


def load(root, names=None):
    root = Path(root)
    record = read_json(root / "evidence.json")
    accessions = record["accessions"]
    count = len(accessions) * (len(accessions) - 1) // 2
    if (
        record.get("format") != FORMAT
        or len(set(accessions)) != len(accessions)
        or record.get("pairs") != count
        or (names is not None and list(names) != accessions)
    ):
        raise AnnotationError("Invalid AAI evidence identity or format")
    arrays = {
        name: np.load(
            root / (name + ".npy"), mmap_mode="r" if count else None, allow_pickle=False
        )
        for name in ARRAYS
    }
    if any(
        array.shape != (count,) or array.dtype != ARRAYS[name]
        for name, array in arrays.items()
    ):
        raise AnnotationError("Invalid AAI evidence array")
    return arrays


def from_pairs(root, names, pairs):
    arrays = create(root, names)
    index = {a: i for i, a in enumerate(names)}
    for (left, right), value in pairs.items():
        k = offset(len(names), index[left], index[right])
        if value is not None:
            arrays["jaccard"][k] = value
            arrays["shared"][k] = 60
            arrays["possible"][k] = 80
    for array in arrays.values():
        if hasattr(array, "flush"):
            array.flush()


def iter_pairs(root):
    from aai_common import aai_estimate

    root = Path(root)
    names = read_json(root / "evidence.json")["accessions"]
    arrays = load(root, names)
    for i, left in enumerate(names):
        for j in range(i + 1, len(names)):
            k = offset(len(names), i, j)
            shared = int(arrays["shared"][k])
            jaccard = float(arrays["jaccard"][k])
            yield dict(
                query=left,
                target=names[j],
                avg_jacc_sim=format(jaccard, ".17g") if shared else "NA",
                num_shared_SCPs=shared,
                poss_shared_SCPs=int(arrays["possible"][k]),
                AAI_estimate=aai_estimate(jaccard) if shared else "NA",
            )


if __name__ == "__main__":
    import argparse
    from annotation_common import write_tsv
    from aai_common import PAIR_COLUMNS

    parser = argparse.ArgumentParser(
        description="Stream a v0.5 AAI evidence directory to TSV."
    )
    parser.add_argument("evidence", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    write_tsv(args.output, PAIR_COLUMNS, iter_pairs(args.evidence))
