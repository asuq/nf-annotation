"""Real child processes obey allocation, queue and deterministic-failure bounds."""

import json
import os
import sys
import tempfile
import time
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from cohort_validation_pool import ordered_validation


def controlled_job(job):
    if job.get("barrier"):
        root = Path(job["barrier"])
        (root / str(os.getpid())).touch()
        deadline = time.monotonic() + 15
        while len(list(root.iterdir())) < job["workers"]:
            if time.monotonic() > deadline:
                raise ValueError("workers did not overlap")
            time.sleep(0.01)
    time.sleep(job.get("delay", 0))
    if job.get("fail"):
        raise ValueError(str(job["value"]))
    return job["value"], os.getpid()


class CohortValidationPoolTests(unittest.TestCase):
    def test_bounded_real_processes_and_ordered_results(self):
        with tempfile.TemporaryDirectory(dir="/tmp") as name:
            root = Path(name)
            (root / "barrier").mkdir()
            jobs = [
                {
                    "value": i,
                    "barrier": str(root / "barrier"),
                    "workers": 4,
                    "delay": (15 - i) * 0.001,
                }
                for i in range(16)
            ]
            results = list(
                ordered_validation(
                    controlled_job,
                    jobs,
                    workers=4,
                    total=16,
                    progress=root / "progress.json",
                )
            )
            self.assertEqual([item[0] for item in results], list(range(16)))
            self.assertEqual(len({item[1] for item in results}), 4)
            metrics = json.loads((root / "progress.json").read_text())
            self.assertEqual(metrics["peak_pending"], 8)
            self.assertEqual(metrics["completed_samples"], 16)
            self.assertEqual(metrics["phase"], "complete")

    def test_serial_equivalence_and_deterministic_first_failure(self):
        jobs = [{"value": i} for i in range(10)]
        for workers in (1, 2, 8):
            results = list(
                ordered_validation(controlled_job, jobs, workers=workers, total=10)
            )
            self.assertEqual([value for value, _ in results], list(range(10)))
            if workers == 1:
                self.assertEqual({pid for _, pid in results}, {os.getpid()})
            with self.assertRaisesRegex(ValueError, "^0$"):
                list(
                    ordered_validation(
                        controlled_job,
                        [
                            {"value": 0, "fail": True, "delay": 0.1},
                            {"value": 1, "fail": True},
                        ],
                        workers=workers,
                        total=2,
                    )
                )
        for workers in (0, 9, True):
            with self.assertRaises(ValueError):
                list(
                    ordered_validation(controlled_job, jobs, workers=workers, total=10)
                )


class CohortValidationEquivalenceTests(unittest.TestCase):
    def test_real_sample_validation_is_byte_identical_for_all_worker_counts(self):
        import copy

        import test_prepare_cohort_update as fixture

        case = fixture.PrepareCohortUpdateTestCase()
        case.setUp()
        self.addCleanup(case.doCleanups)
        case.published_cohort(["ID-X", "B"], failed=True, skipped=True)
        base = case.arguments(["B", "ID-X", "ID.X"])
        expected = None
        for workers in (1, 2, 4, 8):
            args = copy.deepcopy(base)
            args.workers = workers
            args.progress = case.root / f"progress-{workers}.json"
            args.outdir = case.root / f"prepared-{workers}"
            fixture.update.run_prepare(args)
            contents = {
                p.name: p.read_bytes() for p in args.outdir.iterdir() if p.is_file()
            }
            if expected is None:
                expected = contents
            self.assertEqual(contents, expected)
            metrics = json.loads(args.progress.read_text())
            self.assertEqual(metrics["completed_samples"], 3)
            self.assertLessEqual(metrics["peak_pending"], 2 * workers)
        victim = base.genome_inputs / "genome000001"
        victim.write_text(">contig1\nTTTT\n")
        for workers in (1, 4):
            args = copy.deepcopy(base)
            args.workers = workers
            args.outdir = case.root / f"failed-{workers}"
            with self.assertRaisesRegex(
                fixture.update.CohortUpdateError, "Genome content changed"
            ):
                fixture.update.run_prepare(args)
            self.assertFalse(args.outdir.exists())
