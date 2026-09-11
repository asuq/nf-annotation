#!/usr/bin/env python3
"""Compare retained Pfam evidence with an unmodified native PfamScan run."""

import os
import random
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


def run(command, **kwargs):
    result = subprocess.run(command, text=True, capture_output=True, **kwargs)
    if result.returncode:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {command}\n{result.stdout}\n{result.stderr}"
        )
    return result


def data_rows(path):
    return sorted(
        line.split()
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    )


class EvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name)
        cls.database = cls.root / "database"
        cls.database.mkdir()
        cls.native_modules = cls.root / "native_modules"
        native_module = cls.native_modules / "Bio/Pfam/Scan/PfamScan.pm"
        native_module.parent.mkdir(parents=True)
        shutil.copyfile("/opt/nf-annotation/provenance/PfamScan.pm", native_module)
        rng = random.Random(20260911)
        sequence = "".join(rng.choices("ACDEFGHIKLMNPQRSTVWY", k=110))
        models = []
        metadata = []
        for index, name in enumerate(("Model_A", "Model_B", "Model_C"), 1):
            residues = list(sequence)
            # Lower-scoring overlapping alternatives; C is declared nested in A.
            for position in range((index - 1) * 4):
                residues[position * 7] = "A" if residues[position * 7] != "A" else "V"
            alignment = cls.root / f"{name}.sto"
            alignment.write_text(f"# STOCKHOLM 1.0\n{name} {''.join(residues)}\n//\n")
            model = cls.root / f"{name}.hmm"
            run(
                [
                    "hmmbuild",
                    "--amino",
                    "--cpu",
                    "1",
                    "--seed",
                    "42",
                    str(model),
                    str(alignment),
                ]
            )
            text = model.read_text()
            text = text.replace("LENG ", f"ACC   PF{index:05d}.1\nLENG ", 1)
            text = text.replace("STATS ", "GA    5.00 5.00\nSTATS ", 1)
            models.append(text)
            metadata.append(
                f"#=GF ID {name}\n#=GF AC PF{index:05d}.1\n#=GF DE Synthetic control {name}\n#=GF GA 5.00; 5.00;\n#=GF TP Domain\n#=GF ML 110\n#=GF CL CL0001\n"
                + ("#=GF NE Model_C\n" if index == 1 else "")
                + "//\n"
            )
        (cls.database / "Pfam-A.hmm").write_text("".join(models))
        (cls.database / "Pfam-A.hmm.dat").write_text("".join(metadata))
        run(["hmmpress", str(cls.database / "Pfam-A.hmm")])
        cls.fasta = cls.root / "proteins.faa"
        cls.fasta.write_text(
            f">single\n{sequence}\n>repeated\n{sequence}{'P' * 30}{sequence}\n>zero\n{'P' * 110}\n"
        )
        native = cls.root / "native.tsv"
        env = {
            **os.environ,
            "PERL5LIB": f"{cls.native_modules}:{os.environ['PERL5LIB']}",
        }
        env.pop("NF_PFAM_EVIDENCE_DIR", None)
        run(
            [
                "pfam_scan.pl",
                "-fasta",
                str(cls.fasta),
                "-dir",
                str(cls.database),
                "-cpu",
                "1",
                "-outfile",
                str(native),
            ],
            env=env,
        )
        cls.native = native
        cls.evidence = cls.root / "evidence"
        run(
            [
                "pfam_scan_evidence.pl",
                "--fasta",
                str(cls.fasta),
                "--database",
                str(cls.database),
                "--cpus",
                "1",
                "--outdir",
                str(cls.evidence),
            ]
        )

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_native_equivalence(self):
        self.assertEqual(
            data_rows(self.native), data_rows(self.evidence / "pfam.resolved.tsv")
        )

    def test_raw_alternatives_and_hmmer_tables_retained(self):
        self.assertGreater(
            len(data_rows(self.evidence / "pfam.raw.tsv")), len(data_rows(self.native))
        )
        for name in ("hmmscan.domtblout", "hmmscan.tblout"):
            text = (self.evidence / name).read_text()
            self.assertIn("# [ok]", text)
            self.assertIn("Model_B", text)

    def test_repeated_and_nested_domains_preserved(self):
        repeated = [row for row in data_rows(self.native) if row[0] == "repeated"]
        # Two A domains, with the declared C domain nested at each occurrence.
        self.assertEqual([row[6] for row in repeated].count("Model_A"), 2)
        self.assertEqual([row[6] for row in repeated].count("Model_C"), 2)
        self.assertNotIn("Model_B", [row[6] for row in repeated])

    def test_zero_hits_remain_empty(self):
        self.assertNotIn("zero", [row[0] for row in data_rows(self.native)])
        zero = self.root / "zero.faa"
        zero.write_text(f">zero\n{'P' * 110}\n")
        evidence = self.root / "zero-evidence"
        run(
            [
                "pfam_scan_evidence.pl",
                "--fasta",
                str(zero),
                "--database",
                str(self.database),
                "--cpus",
                "1",
                "--outdir",
                str(evidence),
            ]
        )
        self.assertEqual(data_rows(evidence / "pfam.raw.tsv"), [])
        self.assertEqual(data_rows(evidence / "pfam.resolved.tsv"), [])
        self.assertEqual((evidence / "pfam.resolved.json").read_text(), "[]\n")


if __name__ == "__main__":
    unittest.main(verbosity=2)
