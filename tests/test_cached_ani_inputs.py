"""Only proven-equivalent summary rows may retain a completed ANI cache entry."""

import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")


@unittest.skipUnless(NEXTFLOW, "Nextflow is required")
class CachedAniInputsTests(unittest.TestCase):
    def test_content_order_guards_and_bounded_cache_discovery(self):
        with tempfile.TemporaryDirectory(dir="/tmp") as directory:
            root = Path(directory)
            (root / "lib").mkdir()
            shutil.copy2(
                ROOT / "lib/CachedAniInputs.groovy", root / "lib/CachedAniInputs.groovy"
            )
            script = root / "main.nf"
            probe = """
class CachedAniProbe {
static void run(String fixture) {
    def root = java.nio.file.Paths.get(fixture).toRealPath()
    def work = root.resolve('work')
    def previous = work.resolve('ab/1234567890')
    previous.toFile().mkdirs()
    previous.resolve('.exitcode').toFile().text = '0'
    def old = previous.resolve('summary.tsv')
    def fresh = root.resolve('fresh.tsv')
    old.toFile().text = 'accession\\tvalue\\nA\\t1\\nB\\t2\\n'
    fresh.toFile().text = 'accession\\tvalue\\nB\\t2\\nA\\t1\\n'
    assert CachedAniInputs.table(fresh, [previous], 'summary.tsv') == old
    for (def changed : ['accession\\tvalue\\nA\\t9\\nB\\t2\\n', 'accession\\tvalue\\nA\\t1\\nA\\t1\\nB\\t2\\n', 'id\\tvalue\\nA\\t1\\nB\\t2\\n']) {
        fresh.toFile().text = changed
        assert !CachedAniInputs.sameRows(fresh, old)
        def failed = false
        try { CachedAniInputs.table(fresh, [previous], 'summary.tsv') } catch (IllegalArgumentException ignored) { failed = true }
        assert failed
    }
    old.toFile().text = 'accession\\tvalue\\nA\\t"x\\ny"\\nB\\tz\\n'
    fresh.toFile().text = old.toFile().text
    assert CachedAniInputs.sameRows(fresh, old)
    fresh.toFile().text = 'accession\\tvalue\\nB\\tz\\nA\\t"x\\ny"\\n'
    assert !CachedAniInputs.sameRows(fresh, old)
    root.resolve('.nextflow.log.1').toFile().text = "Task completed > TaskHandler[name: COHORT_ANI:BUILD_FASTANI_INPUTS; workDir: ${previous}]\\n"
    assert CachedAniInputs.find(root, work, true, true) == [previous]
    assert CachedAniInputs.find(root, work, false, true) == []
    assert CachedAniInputs.find(root, work, true, false) == []
    previous.resolve('.exitcode').toFile().text = '1'
    def failed = false
    try { CachedAniInputs.find(root, work, true, true) } catch (IllegalArgumentException ignored) { failed = true }
    assert failed
    def a = root.resolve('a.tsv'); def b = root.resolve('b.tsv')
    a.toFile().text = 'accession\\tlineage_a\\nA\\t1\\n'
    b.toFile().text = 'accession\\tlineage_b\\nA\\t2\\n'
    def busco = previous.resolve('busco_tables'); busco.toFile().mkdirs()
    def first = busco.resolve('busco_table01.tsv'); def second = busco.resolve('busco_table02.tsv')
    first.toFile().text = b.toFile().text; second.toFile().text = a.toFile().text
    assert CachedAniInputs.busco([a,b], [previous]) == [first,second]
}
}
"""
            (root / "lib/CachedAniProbe.groovy").write_text(probe)
            script.write_text("workflow { CachedAniProbe.run(params.fixture.toString()) }\n")
            result = subprocess.run(
                [NEXTFLOW, "run", str(script), "--fixture", str(root)],
                cwd=root,
                capture_output=True,
                text=True,
                check=False,
                timeout=120,
                env=dict(
                    os.environ, NXF_DISABLE_CHECK_LATEST="true", NXF_ANSI_LOG="false"
                ),
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
