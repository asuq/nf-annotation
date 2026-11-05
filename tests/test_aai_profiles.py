"""Eligibility and portable marker caches do not repeat native preprocessing."""

from __future__ import annotations
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import prepare_aai
from annotation_common import AnnotationError, digest, read_json, write_json, write_tsv


class AaiProfileTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(dir="/tmp")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        runtime = self.root / "runtime.json"
        write_json(
            runtime,
            dict(
                upstream_commit="fixture",
                hmm_sha256="fixture",
                pyhmmer="fixture",
                numpy="fixture",
                environment_sha256="fixture",
            ),
        )
        self.patch = patch.object(prepare_aai, "RUNTIME", runtime)
        self.patch.start()
        self.addCleanup(self.patch.stop)
        self.bundle = self.root / "bundle"
        self.bundle.mkdir()
        (self.bundle / "proteins.faa").write_text(">g\nMWA\n")
        write_json(
            self.bundle / "bundle.json",
            dict(
                accession="A",
                status="success",
                files={"proteins.faa": digest(self.bundle / "proteins.faa")},
            ),
        )
        self.args = SimpleNamespace(
            **{
                key: self.root / (key + ".json")
                for key in (
                    "records",
                    "metadata",
                    "checkm2",
                    "busco",
                    "stats",
                    "bundles",
                    "sources",
                    "output",
                )
            },
            primary_busco_column="BUSCO_fixture",
            domain="bacteria",
            threads=2,
            destination=self.root / "result",
        )
        self.args.destination.mkdir()
        write_json(
            self.args.records,
            [dict(accession="A", is_new="false", assembly_level="Scaffold")],
        )
        write_tsv(
            self.args.metadata,
            ("accession", "Assembly_Level", "Atypical_Warnings", "16S"),
            [
                dict(
                    accession="A",
                    Assembly_Level="Scaffold",
                    Atypical_Warnings="NA",
                    **{"16S": "No"},
                )
            ],
        )
        write_tsv(
            self.args.checkm2,
            (
                "accession",
                "Gcode",
                "Low_quality",
                "Completeness_gcode4",
                "Contamination_gcode4",
            ),
            [
                dict(
                    accession="A",
                    Gcode="4",
                    Low_quality="false",
                    Completeness_gcode4="95",
                    Contamination_gcode4="1",
                )
            ],
        )
        write_tsv(
            self.args.stats,
            ("accession", "n50", "scaffolds", "genome_size", "gc_content"),
            [
                dict(
                    accession="A",
                    n50="NA",
                    scaffolds="1",
                    genome_size="1200",
                    gc_content="50",
                )
            ],
        )
        busco = self.root / "busco.tsv"
        write_tsv(
            busco,
            ("accession", "BUSCO_fixture"),
            [dict(accession="A", BUSCO_fixture="NA")],
        )
        write_json(self.args.busco, [str(busco)])
        write_json(self.args.bundles, [str(self.bundle)])
        write_json(self.args.sources, [])
        self.calls = 0
        owner = self

        class Pool:
            def __init__(self, **kwargs):
                pass

            def __enter__(self):
                return self

            def __exit__(self, *args):
                pass

            def map(self, function, jobs):
                for _, directory in jobs:
                    owner.calls += 1
                    (directory / "markers.faa").write_text(">g\nMWA\n")
                    (directory / "markers.hmm").write_text(
                        "native fixture HMM evidence\n"
                    )
                    yield 1, directory

        self.pool = patch.object(
            prepare_aai.concurrent.futures, "ProcessPoolExecutor", Pool
        )
        self.pool.start()
        self.addCleanup(self.pool.stop)

    def test_no_16s_or_scoring_gate_and_same_inputs_do_not_repeat_search(self):
        prepare_aai.prepare(self.args)
        first = read_json(self.args.output)["samples"][0]
        self.assertEqual(first["status"], "ready")
        self.assertEqual(first["quality"]["n50"], "NA")
        self.assertEqual(self.calls, 1)
        prepare_aai.prepare(self.args)
        self.assertEqual(self.calls, 1)
        self.assertEqual(first, read_json(self.args.output)["samples"][0])
        self.assertEqual(
            read_json(self.args.output.parent / "preprocessing_metrics.json")["reused"],
            1,
        )

    def test_successor_owns_profiles_after_source_removal(self):
        import shutil

        prepare_aai.prepare(self.args)
        source = self.args.destination
        self.args.destination = self.root / "successor"
        self.args.destination.mkdir()
        write_json(self.args.sources, [str(source)])
        prepare_aai.prepare(self.args)
        self.assertEqual(self.calls, 1)
        path = Path(read_json(self.args.output)["samples"][0]["profile"])
        shutil.rmtree(source)
        self.assertIsNotNone(prepare_aai.read_profile(path.parent, path.name))
        self.assertFalse(any(p.is_symlink() for p in path.iterdir()))

    def test_interrupted_profile_copy_can_resume_without_marker_extraction(self):
        prepare_aai.prepare(self.args)
        source = self.args.destination
        self.args.destination = self.root / "successor"
        self.args.destination.mkdir()
        write_json(self.args.sources, [str(source)])
        with patch.object(prepare_aai, "copy_independent", side_effect=OSError("interrupted copy")):
            with self.assertRaisesRegex(OSError, "interrupted copy"):
                prepare_aai.prepare(self.args)
        profile_root = self.args.destination / "cohort/fastaai/profiles"
        self.assertEqual(list(profile_root.iterdir()), [])
        prepare_aai.prepare(self.args)
        self.assertEqual(self.calls, 1)
        self.assertEqual(read_json(self.args.output.parent / "preprocessing_metrics.json")["reused"], 1)

    def test_changed_protein_or_method_recomputes_and_changed_cache_fails(self):
        prepare_aai.prepare(self.args)
        p = Path(read_json(self.args.output)["samples"][0]["profile"])
        (p / "markers.faa").write_text("changed")
        with self.assertRaises(AnnotationError):
            prepare_aai.prepare(self.args)
        (p / "markers.faa").write_text(">g\nMWA\n")
        record = read_json(self.bundle / "bundle.json")
        (self.bundle / "proteins.faa").write_text(">g\nMWG\n")
        record["files"]["proteins.faa"] = digest(self.bundle / "proteins.faa")
        write_json(self.bundle / "bundle.json", record)
        prepare_aai.prepare(self.args)
        self.assertEqual(self.calls, 2)
        runtime = read_json(prepare_aai.RUNTIME)
        runtime["hmm_sha256"] = "different"
        write_json(prepare_aai.RUNTIME, runtime)
        prepare_aai.prepare(self.args)
        self.assertEqual(self.calls, 3)


if __name__ == "__main__":
    unittest.main()
