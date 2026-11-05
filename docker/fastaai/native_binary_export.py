"""Serialize native all-versus-all results without text or duplicate pair records."""

from pathlib import Path
import numpy as np


def write_row(output, index, jaccard, shared, possible):
    index = int(index)
    if index + 1 >= len(jaccard):
        return
    selection = slice(index + 1, None)
    np.savez(
        Path(output) / "results" / f"{index:08d}.npz",
        jaccard=np.asarray(jaccard[selection], dtype=np.float64),
        shared=np.asarray(shared[selection], dtype=np.uint8),
        possible=np.asarray(possible[selection], dtype=np.uint8),
    )
