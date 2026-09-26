"""Safety and accounting controls for annotation-owned temporary storage."""

from __future__ import annotations

import json
import sys
import subprocess
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
def finalize(args):
    command = ['bash', str(Path(__file__).resolve().parents[1] / 'bin/annotation_storage.sh'), 'finalize']
    for name in ('execution_root', 'durable_work', 'task_directory', 'tool', 'batch_id', 'native_exit', 'policy'):
        command.extend(['--' + name.replace('_', '-'), str(getattr(args, name))])
    result = subprocess.run(command, text=True, capture_output=True)
    if result.returncode:
        raise RuntimeError(result.stderr)
    return result.returncode


class Args:
    command = "finalize"
    task_directory = "task00000001"
    tool = "kofam"
    batch_id = "NA"
    policy = "success"
    native_exit = 0


class AnnotationStorageTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.execution = self.root / "execution"
        self.durable = self.root / "durable"
        self.execution.mkdir()
        self.durable.mkdir()
        self.raw = self.execution / "raw"
        self.raw.mkdir()
        (self.raw / "exit_code.txt").write_text("0\n")
        (self.raw / "evidence.tsv").write_text("header\nvalue\n")
        for name in ('tool.log', 'versions.txt', 'kofam.tsv'):
            (self.raw / name).write_text('')
        self.scratch = self.execution / "scratch"
        self.scratch.mkdir()
        (self.scratch / "temporary.bin").write_bytes(b"x" * 4096)

    def args(self, **changes):
        args = Args()
        args.execution_root = self.execution
        args.durable_work = self.durable
        for key, value in changes.items():
            setattr(args, key, value)
        return args

    def test_success_removes_only_owned_scratch_and_records_bytes(self):
        self.assertEqual(finalize(self.args()), 0)
        self.assertFalse(self.scratch.exists())
        self.assertEqual((self.raw / "evidence.tsv").read_text(), "header\nvalue\n")
        report = json.loads((self.execution / "annotation_storage.json").read_text())
        self.assertEqual(report["cleanup_status"], "removed")
        self.assertGreaterEqual(report["scratch_before"]["allocated_bytes"], 4096)
        self.assertEqual(report["scratch_after"]["allocated_bytes"], 0)

    def test_cleanup_off_preserves_scratch(self):
        self.assertEqual(finalize(self.args(policy="off")), 0)
        self.assertTrue((self.scratch / "temporary.bin").is_file())
        report = json.loads((self.execution / "annotation_storage.json").read_text())
        self.assertEqual(report["cleanup_status"], "disabled")
        self.assertEqual(report["scratch_after"], report["scratch_before"])

    def test_native_failure_exports_raw_and_preserves_scratch(self):
        (self.raw / "exit_code.txt").write_text("7\n")
        self.assertEqual(finalize(self.args(native_exit=7)), 0)
        self.assertTrue((self.scratch / "temporary.bin").is_file())
        self.assertEqual((self.durable / "raw/evidence.tsv").read_text(), "header\nvalue\n")
        report = json.loads((self.durable / "annotation_storage.json").read_text())
        self.assertEqual(report["failure_export"], "exported")
        self.assertEqual(report["native_exit_status"], 7)

    def test_symlink_under_scratch_is_refused_without_touching_target(self):
        sentinel = self.root / "external.txt"
        sentinel.write_text("keep")
        (self.scratch / "escape").symlink_to(sentinel)
        self.assertEqual(finalize(self.args()), 0)
        report = json.loads((self.execution / "annotation_storage.json").read_text())
        self.assertEqual(report["cleanup_status"], "failed")
        self.assertIn("Unsafe", report["cleanup_error"])
        self.assertEqual(sentinel.read_text(), "keep")
        self.assertTrue(self.scratch.exists())

    def test_exit_disagreement_fails_before_cleanup(self):
        (self.raw / "exit_code.txt").write_text("9\n")
        with self.assertRaisesRegex(RuntimeError, "disagrees"):
            finalize(self.args())
        self.assertTrue(self.scratch.exists())

    def test_failure_export_is_idempotent(self):
        (self.raw / 'exit_code.txt').write_text('7\n')
        (self.raw / 'empty').mkdir()
        finalize(self.args(native_exit=7))
        finalize(self.args(native_exit=7))
        report = json.loads((self.durable / 'annotation_storage.json').read_text())
        self.assertEqual(report['failure_export'], 'exported')
        self.assertTrue((self.durable / 'raw/empty').is_dir())

    def test_failure_export_is_not_blocked_by_unsafe_scratch(self):
        (self.raw / 'exit_code.txt').write_text('7\n')
        (self.scratch / 'external').symlink_to(self.durable, target_is_directory=True)
        finalize(self.args(native_exit=7))
        self.assertTrue((self.durable / 'raw/evidence.tsv').is_file())
        self.assertTrue((self.scratch / 'external').is_symlink())

    def test_missing_exit_record_exports_partial_diagnostics(self):
        (self.raw / 'exit_code.txt').unlink()
        finalize(self.args(native_exit=143))
        report = json.loads((self.durable / 'annotation_storage.json').read_text())
        self.assertFalse(report['diagnostics_complete'])
        self.assertEqual(report['native_exit_status'], 143)
        self.assertTrue((self.durable / 'raw/evidence.tsv').is_file())

    def test_linked_scratch_root_is_never_removed(self):
        import shutil
        shutil.rmtree(self.scratch)
        self.scratch.symlink_to(self.durable, target_is_directory=True)
        (self.durable / 'sentinel').write_text('keep')
        finalize(self.args())
        self.assertEqual((self.durable / 'sentinel').read_text(), 'keep')
        self.assertTrue(self.scratch.is_symlink())

    def test_missing_artifact_keeps_scratch_without_masking_parser_failure(self):
        (self.raw / 'kofam.tsv').unlink()
        finalize(self.args())
        self.assertTrue(self.scratch.is_dir())
        report = json.loads((self.execution / 'annotation_storage.json').read_text())
        self.assertEqual(report['cleanup_status'], 'failed')
        self.assertEqual(report['native_exit_status'], 0)

    def test_sparse_file_accounting_and_many_files(self):
        with (self.scratch / 'sparse').open('wb') as handle:
            handle.truncate(64 * 1024 * 1024)
        for index in range(1000):
            (self.scratch / str(index)).write_text('x')
        finalize(self.args())
        report = json.loads((self.execution / 'annotation_storage.json').read_text())
        self.assertEqual(report['scratch_before']['files'], 1002)
        self.assertGreater(report['scratch_before']['apparent_bytes'], report['scratch_before']['allocated_bytes'])


if __name__ == "__main__":
    unittest.main()
