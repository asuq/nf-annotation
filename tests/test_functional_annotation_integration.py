"""Exercise real Nextflow annotation orchestration with synthetic native output.

Database-scale native-tool qualification is separate. These controls validate
staging, identity-dependent reuse, failure publication and the final gate.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
import test_annotation_bundle as bundle_fixture
from annotation_common import (
    digest,
    identity,
    read_json,
    read_tsv,
    write_json,
    write_tsv,
)
from annotation_resources import file_records
from annotation_source import validate_source
from annotation_tasks import PADLOC_RESOURCE, plan_task, preflight
from annotation_workflow_fixture import publish_disabled

NEXTFLOW = shutil.which("nextflow")


@unittest.skipUnless(NEXTFLOW, "Nextflow is required for annotation integration")
class FunctionalAnnotationIntegrationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(
            prefix="nf-annotation-functional-", dir="/tmp"
        )
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.project = self.root / "pipeline"
        self.project.mkdir()
        for name in ("bin", "modules", "subworkflows", "conf", "assets"):
            shutil.copytree(ROOT / name, self.project / name)
        for name in ("reannotate.nf", "nextflow.config"):
            shutil.copyfile(ROOT / name, self.project / name)
        shutil.copytree(
            ROOT / "external/nf-helper/conf", self.project / "external/nf-helper/conf"
        )
        self.source = self.root / "source"
        (self.source / "tables").mkdir(parents=True)
        fixture = bundle_fixture.AnnotationBundleTests()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        bundles = []
        for accession in ("B", "A"):
            bundle = self.source / "samples" / accession / "annotation/bundle"
            shutil.copytree(fixture.build(accession, name=accession), bundle)
            bundles.append(bundle)
        write_tsv(
            self.source / "tables/master_table.tsv",
            ("Accession", "Gcode", "upstream_fixture"),
            [dict(Accession=acc, Gcode=4, upstream_fixture=acc) for acc in ("B", "A")],
        )
        write_tsv(
            self.source / "tables/sample_status.tsv",
            ("accession", "gcode_status"),
            [dict(accession=acc, gcode_status="done") for acc in ("B", "A")],
        )
        write_tsv(
            self.source / "tables/validated_samples.tsv",
            ("accession",),
            [dict(accession=acc) for acc in ("B", "A")],
        )
        publish_disabled(self.source, ["B", "A"], bundles)
        self.database = self.root / "kofam"
        self.database.mkdir()
        write_tsv(
            self.database / "ko_list",
            ("knum", "threshold", "score_type", "definition"),
            [
                dict(
                    knum="K00001",
                    threshold="100",
                    score_type="full",
                    definition="Synthetic control",
                )
            ],
        )
        write_tsv(
            self.database / "prokaryote_profiles.tsv",
            ("ko", "profile", "threshold", "score_type"),
            [
                dict(
                    ko="K00001",
                    profile="K00001.hmm",
                    threshold="100",
                    score_type="full",
                )
            ],
        )
        (self.database / "profiles").mkdir()
        (self.database / "profiles/prokaryote.hal").write_text("K00001.hmm\n")
        self.resource("fixture1")
        wrapper = self.project / "bin/exec_annotation"
        wrapper.write_text(
            f"#!{sys.executable}\n"
            + """import json
import sys
from pathlib import Path
if sys.argv[1:] == ['--version']:
    print('Synthetic KOfam workflow control')
    sys.exit(0)
task = json.loads(Path('task/task.json').read_text())
if Path(__file__).with_name('fail_b').exists() and task['accession'] == 'B':
    sys.exit(4)
query = Path(sys.argv[-1]).read_text().splitlines()[0][1:].split()[0]
text = '#\\tgene name\\tKO\\tthrshld\\tscore\\tE-value\\tKO definition\\n'
text += '#\\t---------\\t------\\t-------\\t------\\t---------\\t-------------\\n'
if task['accession'] == 'A':
    text += f'*\\t{query}\\tK00001\\t100.00\\t101.0\\t1e-20\\tSynthetic control\\n'
Path(sys.argv[sys.argv.index('-o') + 1]).write_text(text)
"""
        )
        wrapper.chmod(0o755)
        hmmsearch = self.project / "bin/hmmsearch"
        hmmsearch.write_text("#!/bin/sh\nprintf 'Synthetic HMMER version control\\n'\n")
        hmmsearch.chmod(0o755)

    def resource(self, version):
        contract = dict(
            component="kofam",
            version=version,
            settings={"fixture": True},
            files=file_records(self.database),
        )
        write_json(
            self.database / "annotation_resource.json",
            dict(
                schema_version=1,
                component="kofam",
                contract=contract,
                resource_id=identity(contract),
            ),
        )

    def run_pipeline(self, name, source, *, enabled=True, succeeds=True, resume=False):
        output = self.root / name
        configuration = self.root / f"{name}.config"
        configuration.write_text(
            "params {\n"
            + f" annotation_tools = '{'kofam' if enabled else ''}'\n"
            + f" kofam_db = '{self.database}'\n"
            + " kofam_container = 'sha256:"
            + "1" * 64
            + "'\n"
            + " python_container = 'sha256:"
            + "2" * 64
            + "'\n"
            + " annotation_cpus = 1\n annotation_memory = '1 GB'\n}\n"
        )
        result = subprocess.run(
            [
                NEXTFLOW,
                "run",
                str(self.project / "reannotate.nf"),
                "-profile",
                "test",
                "-c",
                str(configuration),
                "--annotation_from",
                str(source),
                "--outdir",
                str(output),
                "-work-dir",
                str(self.root / f"work-{name}"),
            ]
            + (["-resume"] if resume else []),
            cwd=self.root,
            text=True,
            capture_output=True,
            timeout=120,
            env=dict(os.environ, NXF_ANSI_LOG="false", NXF_DISABLE_CHECK_LATEST="true"),
        )
        (self.root / f"{name}.log").write_text(result.stdout + result.stderr)
        self.assertEqual(
            result.returncode == 0, succeeds, result.stdout + result.stderr
        )
        manifest = validate_source(output)
        self.assertEqual(manifest["complete"], succeeds)
        return output

    def actions(self, source):
        return [
            row["action"]
            for row in read_tsv(source / "pipeline_info/annotation_plan.tsv")
            if row["tool"] == "kofam"
        ]

    def test_padloc_stages_the_bundle_without_an_external_resource(self):
        entry = preflight(
            dict(
                annotation_tools="padloc",
                helper_container="sha256:" + "2" * 64,
                tools=dict(
                    padloc=dict(container="sha256:" + "1" * 64, cpus=1, memory_gib=1)
                ),
            )
        )["tools"]["padloc"]
        bundle = self.source / "samples/A/annotation/bundle"
        directory = self.root / "padloc-task"
        plan_task(bundle, "padloc", entry, directory)
        raw = directory / "native_control"
        raw.mkdir()
        (raw / "input.domtblout").write_text("# Program: hmmsearch\n# [ok]\n")
        (raw / "tool.log").write_text("Nothing found for input\n")
        (raw / "exit_code.txt").write_text("0\n")
        (raw / "versions.txt").write_text(
            f"commit={PADLOC_RESOURCE['commit']}\n"
            f"{PADLOC_RESOURCE['archive_sha256']}  /opt/nf-annotation/provenance/padloc-db.tar.gz\n"
        )
        (directory / "run.sh").write_text(
            "#!/bin/bash\nset -euo pipefail\ncp -R task/native_control raw\n"
        )
        script = self.project / "padloc_control.nf"
        script.write_text("""nextflow.enable.dsl=2
include { ANNOTATION_SEARCH; NORMALIZE_ANNOTATION } from './modules/local/functional_annotation'
workflow {
    meta = [accession:'A', tool:'padloc', cpus:1, memory_gib:1, container:'fixture']
    ANNOTATION_SEARCH(Channel.of(tuple(meta, file(params.fixture_task), file(params.fixture_bundle), [])))
    NORMALIZE_ANNOTATION(ANNOTATION_SEARCH.out.raw_results)
}
""")
        output = self.root / "padloc-output"
        result = subprocess.run(
            [
                NEXTFLOW,
                "run",
                str(script),
                "-profile",
                "test",
                "--fixture_task",
                str(directory),
                "--fixture_bundle",
                str(bundle),
                "--outdir",
                str(output),
                "-work-dir",
                str(self.root / "padloc-work"),
            ],
            cwd=self.root,
            text=True,
            capture_output=True,
            timeout=120,
            env=dict(os.environ, NXF_ANSI_LOG="false", NXF_DISABLE_CHECK_LATEST="true"),
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        record = read_json(output / "samples/A/annotation/padloc/result.json")
        self.assertEqual(record["status"], "success", record["reason"])
        self.assertEqual(record["accepted_proteins"], 0)

    def test_run_reuse_renormalize_resource_change_and_failed_gate(self):
        original = {
            str(path.relative_to(self.source)): digest(path)
            for path in self.source.rglob("*")
            if path.is_file()
        }
        run = self.run_pipeline("initial", self.source)
        self.assertEqual(self.actions(run), ["run", "run"])
        self.assertEqual(
            read_tsv(run / "tables/functional_matrices/kofam_ko_counts.tsv"),
            [dict(accession="B", K00001="0"), dict(accession="A", K00001="1")],
        )
        self.run_pipeline("initial", self.source, resume=True)
        searches = [
            row
            for row in read_tsv(run / "pipeline_info/trace.tsv")
            if ":ANNOTATION_SEARCH " in row["name"]
        ]
        self.assertEqual(len(searches), 2)
        self.assertEqual({row["status"] for row in searches}, {"CACHED"})
        shutil.rmtree(self.root / "work-initial")
        reused = self.run_pipeline("reused", run)
        self.assertEqual(self.actions(reused), ["reuse", "reuse"])
        parser = self.project / "bin/normalize_kofam.py"
        parser.write_text(
            parser.read_text() + "\n# Interpretive-code change control.\n"
        )
        normalized = self.run_pipeline("renormalized", reused)
        self.assertEqual(self.actions(normalized), ["renormalize", "renormalize"])
        self.resource("fixture2")
        changed = self.run_pipeline("changed", normalized)
        self.assertEqual(self.actions(changed), ["run", "run"])
        self.resource("fixture3")
        (self.project / "bin/fail_b").touch()
        failed = self.run_pipeline("failed", changed, succeeds=False)
        states = {
            row["accession"]: row["status"]
            for row in read_tsv(failed / "tables/annotation_status.tsv")
            if row["tool"] == "kofam"
        }
        self.assertEqual(states, dict(A="success", B="failed"))
        self.assertEqual(
            read_tsv(failed / "tables/functional_matrices/kofam_ko_counts.tsv"),
            [dict(accession="B", K00001="NA"), dict(accession="A", K00001="1")],
        )
        self.assertEqual(
            read_json(failed / "samples/B/annotation/kofam/result.json")["exit_code"], 4
        )
        disabled = self.run_pipeline("disabled", failed, enabled=False)
        self.assertEqual(self.actions(disabled), ["skip", "skip"])
        self.assertEqual(
            read_tsv(disabled / "tables/functional_matrices/kofam_ko_counts.tsv"),
            [dict(accession="B"), dict(accession="A")],
        )
        self.assertEqual(
            len(read_tsv(disabled / "tables/protein_function_summary.tsv")), 2
        )
        self.assertEqual(
            original,
            {
                str(path.relative_to(self.source)): digest(path)
                for path in self.source.rglob("*")
                if path.is_file()
            },
        )
