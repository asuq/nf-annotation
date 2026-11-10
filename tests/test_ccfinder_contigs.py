"""Uniform, byte-preserving FASTA aliases and strict native result restoration."""

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from typing import Any

SCRIPT = Path(__file__).resolve().parents[1] / "bin/ccfinder_contigs.pl"
sys.path.insert(0, str(SCRIPT.parent))
from ccfinder_mapping import validate_publication


class ContigInputTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / "source.fa"
        self.source.write_bytes(
            b">same.1 description\twith tabs\r\nAcgTN\r\n>same.2 another\nACGTN\n"
        )

    def call(self, mode, *args, ok=True):
        result = subprocess.run(
            ["perl", str(SCRIPT), mode, *map(str, args)],
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode == 0, ok, result.stderr.decode())
        return result

    def stage(self):
        self.call(
            "stage",
            self.source,
            self.root / "tool.fa",
            self.root / "map.json",
            self.root,
            self.root / "files.txt",
        )
        return json.loads((self.root / "map.json").read_text())

    def restore(self, result, ok=True):
        path = self.root / "native.json"
        path.write_text(json.dumps(result))
        self.call(
            "restore",
            self.source,
            self.root / "tool.fa",
            self.root / "map.json",
            path,
            self.root / "restored.json",
            ok=ok,
        )

    def test_tabs_dots_sequence_bytes_and_true_zero(self):
        original = self.source.read_bytes()
        mapping = self.stage()
        self.assertEqual(
            (self.root / "tool.fa").read_bytes(),
            b">ctg00000001\nAcgTN\r\n>ctg00000002\nACGTN\n",
        )
        self.assertEqual(mapping["source_fasta_md5"], hashlib.md5(original).hexdigest())
        self.restore(
            {
                "Sequences": [
                    {
                        "Id": r["tool_id"],
                        "Length": r["length"],
                        "Crisprs": [],
                        "Cas": [],
                    }
                    for r in mapping["contigs"]
                ]
            }
        )
        result = json.loads((self.root / "restored.json").read_text())
        self.assertEqual([r["Id"] for r in result["Sequences"]], ["same.1", "same.2"])
        self.assertEqual(
            result["Sequences"][0]["Description"], "description\twith tabs"
        )
        self.assertIsInstance(result["Sequences"][0]["Length"], int)
        self.assertEqual(self.source.read_bytes(), original)

    def test_omission_duplicates_wrong_lengths_and_invalid_json_fail(self):
        self.stage()
        good = {"Id": "ctg00000001", "Length": 5, "Crisprs": [], "Cas": []}
        for rows in (
            [],
            [good],
            [good, good],
            [good, {"Id": "ctg00000002", "Length": 4, "Crisprs": [], "Cas": []}],
        ):
            self.restore({"Sequences": rows}, ok=False)
        path = self.root / "native.json"
        path.write_text('{"Sequences": [')
        self.call(
            "restore",
            self.source,
            self.root / "tool.fa",
            self.root / "map.json",
            path,
            self.root / "restored.json",
            ok=False,
        )

    def test_changed_mapping_or_sequence_fails(self):
        mapping = self.stage()
        rows: list[dict[str, Any]] = [
            {"Id": r["tool_id"], "Length": r["length"], "Crisprs": [], "Cas": []}
            for r in mapping["contigs"]
        ]
        (self.root / "tool.fa").write_bytes(
            (self.root / "tool.fa").read_bytes().replace(b"AcgTN", b"AcgTA")
        )
        self.restore({"Sequences": rows}, ok=False)

    def test_native_numeric_lengths_keep_their_json_type(self):
        self.source.write_text(">large\n" + "A" * 32150 + "\n")
        self.stage()
        self.restore(
            {
                "Sequences": [
                    {"Id": "ctg00000001", "Length": 32150, "Crisprs": [], "Cas": []}
                ]
            }
        )
        value = json.loads((self.root / "restored.json").read_text())["Sequences"][0][
            "Length"
        ]
        self.assertIs(type(value), int)

    def test_duplicate_or_invalid_source_records_fail(self):
        for data in (b">a\nACGT\n>a\nACGT\n", b">a\nAC GT\n", b">a\n"):
            self.source.write_bytes(data)
            self.call(
                "stage",
                self.source,
                self.root / "tool.fa",
                self.root / "map.json",
                self.root,
                self.root / "files.txt",
                ok=False,
            )

    def test_nested_identifiers_order_and_publication_replay(self):
        self.stage()
        rows: list[dict[str, Any]] = [
            {
                "Id": f"ctg{i:08d}",
                "Length": 5,
                "Cas": [],
                "Crisprs": [
                    {"Name": f"ctg{i:08d}_1", "Start": 1, "End": 5, "Sequence": "ACGTN"}
                ],
                "Summary_CRISPR-Cas": f"(ctg{i:08d}_1[1;5])",
            }
            for i in (2, 1)
        ]
        self.restore({"Sequences": rows})
        result = json.loads((self.root / "restored.json").read_text())
        self.assertEqual(
            [row["Id"] for row in result["Sequences"]], ["same.1", "same.2"]
        )
        self.assertEqual(result["Sequences"][0]["Crisprs"][0]["Name"], "same.1_1")
        self.assertEqual(
            result["Sequences"][0]["Summary_CRISPR-Cas"], "(same.1_1[1;5])"
        )
        for original, published in (
            ("restored.json", "result.json"),
            ("native.json", "native_result.json"),
            ("map.json", "contig_mapping.json"),
        ):
            (self.root / published).write_bytes((self.root / original).read_bytes())
        validate_publication(self.root)
        result["Sequences"][0]["Crisprs"][0]["Start"] = 2
        (self.root / "result.json").write_text(json.dumps(result))
        with self.assertRaisesRegex(ValueError, "differs"):
            validate_publication(self.root)
        rows[0]["Crisprs"][0]["Name"] = "ctg00000001_1"
        self.restore({"Sequences": rows}, ok=False)
