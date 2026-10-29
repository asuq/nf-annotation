"""Bounded sample validation with ordered results and coordinator-owned progress."""

import multiprocessing
import sys
import time
from collections import deque
from concurrent.futures import ProcessPoolExecutor, TimeoutError
from pathlib import Path

from annotation_common import write_json


def ordered_validation(
    function, jobs, *, workers: int, total: int, progress: Path | None = None
):
    if type(workers) is not int or not 1 <= workers <= 8:
        raise ValueError("Cohort validation workers must be an integer from 1 to 8")
    start = time.monotonic()
    last = 0.0
    metrics = {
        "phase": "validating",
        "workers": workers,
        "total_samples": total,
        "completed_samples": 0,
        "peak_pending": 0,
    }

    def report(force=False):
        nonlocal last
        now = time.monotonic()
        metrics["seconds"] = now - start
        if not force and now - last < 10:
            return
        last = now
        if progress is not None:
            temporary = progress.with_name(progress.name + ".tmp")
            write_json(temporary, metrics)
            temporary.replace(progress)
            print(f"Cohort validation: {metrics}", file=sys.stderr, flush=True)

    report(True)
    pool = None
    try:
        if workers == 1:
            for job in jobs:
                result = function(job)
                metrics["completed_samples"] += 1
                metrics["peak_pending"] = 1
                report()
                yield result
        else:
            pool = ProcessPoolExecutor(
                max_workers=workers, mp_context=multiprocessing.get_context("spawn")
            )
            pending = deque()
            iterator = iter(jobs)
            exhausted = False
            while pending or not exhausted:
                while not exhausted and len(pending) < 2 * workers:
                    try:
                        job = next(iterator)
                    except StopIteration:
                        exhausted = True
                        break
                    pending.append(pool.submit(function, job))
                    metrics["peak_pending"] = max(metrics["peak_pending"], len(pending))
                if not pending:
                    break
                future = pending.popleft()
                while True:
                    try:
                        result = future.result(timeout=1)
                        break
                    except TimeoutError:
                        if future.done():
                            raise
                        report()
                metrics["completed_samples"] += 1
                report()
                yield result
        metrics["phase"] = "complete"
    except BaseException:
        metrics["phase"] = "failed"
        raise
    finally:
        if pool is not None:
            pool.shutdown(wait=True, cancel_futures=True)
        report(True)
