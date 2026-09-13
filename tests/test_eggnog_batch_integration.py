"""Exercise shared eggNOG orchestration with complete synthetic native evidence.

These tests run the actual planner, Nextflow processes, aggregation and portable
source validation. They do not qualify biological search sensitivity.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import unittest
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))

import test_annotation_bundle as bundle_fixture
import test_annotation_normalizers as normalizer_fixture
import test_functional_annotation_integration as workflow_fixture
from annotation_common import (
    bundle_proteins,
    identity,
    read_json,
    read_tsv,
    write_json,
    write_tsv,
)
from annotation_resources import file_records
from annotation_source import validate_source
from annotation_workflow_fixture import publish_disabled

NEXTFLOW = shutil.which("nextflow")

EMAPPER = r"""import csv
import json
import sys
from pathlib import Path
import normalize_eggnog as egg

if sys.argv[1:] == ['--version']:
    print('emapper-3.0.0-beta6 synthetic orchestration control')
    sys.exit(0)
source = Path(sys.argv[sys.argv.index('-i') + 1])
mode = sys.argv[sys.argv.index('-m') + 1]
proteins = []
for line in source.read_text().splitlines():
    if line.startswith('>'):
        proteins.append(dict(query=line[1:].split()[0], sequence=''))
    else:
        proteins[-1]['sequence'] += line
input_count = len(proteins)
with Path(__LEDGER__).open('a') as handle:
    handle.write(json.dumps(dict(mode=mode, queries=[protein['query'] for protein in proteins])) + '\n')
if Path(__FAIL__).exists() and Path(__FAIL__).read_text().strip() in ('', mode):
    print('Intentional synthetic whole-batch native failure', file=sys.stderr)
    sys.exit(9)
output = Path(sys.argv[sys.argv.index('--output_dir') + 1])
prefix = sys.argv[sys.argv.index('-o') + 1]
if mode == 'diamond':
    omitted = set(json.loads(Path(__NOHITS__).read_text())) if Path(__NOHITS__).exists() else set()
    proteins = [protein for protein in proteins if protein['query'] not in omitted]
    with (output / (prefix + '.emapper.seed_orthologs')).open('w') as handle:
        handle.write('#' + '\t'.join(egg.SEED_NATIVE_COLUMNS) + '\n')
        for protein in proteins:
            length = len(protein['sequence'])
            handle.write(f"{protein['query']}\t1\t0\t100\t1\t{length}\t1\t{length}\t90\t100\t100\n")
        handle.write(f'## {len(proteins)} queries scanned\n')
    (output / (prefix + '.emapper.hits')).write_text('synthetic native search diagnostics\n')
    if '--no_annot' in sys.argv and '--report_orthologs' not in sys.argv:
        sys.exit(0)
elif mode == 'no_search':
    source_seeds = Path(sys.argv[sys.argv.index('--annotate_hits_table') + 1])
    seeds = [line.split('\t') for line in source_seeds.read_text().splitlines()
             if line and not line.startswith('#')]
    seeds.sort(key=lambda row: (int(row[1]), row[0]))
    Path(str(source_seeds) + '.sorted').write_text(''.join('\t'.join(row) + '\n' for row in seeds))
    identifiers = {row[0] for row in seeds}
    assert identifiers <= {protein['query'] for protein in proteins}
    proteins = [protein for protein in proteins if protein['query'] in identifiers]
else:
    raise ValueError('Unexpected synthetic mapper mode')
rows, namespaces = [], []
for protein in proteins:
    query = protein['query']
    row = dict.fromkeys(egg.HEADER, '-')
    row.update(query=query, seed_ortholog='1234', evalue='0', score='100',
               COG_category='J', GOs='GO:0000001', EC='ec:1.2.3.4', KEGG_ko='K00001',
               Preferred_name='hflB')
    # Model the observed native medium/low name change under pooled context.
    row['annotation_confidence'] = ''.join(
        ('m' if input_count == 1 else 'l') if name == 'Preferred_name'
        else 'h' if name in ('GOs', 'EC', 'KEGG_ko') else '-'
        for name in egg.FIELDS
    )
    rows.append(row)
    namespace = dict.fromkeys(egg.GO_HEADER, '-')
    namespace.update(query=query, gos_mf='GO:0000001', gos_mf_confidence='high')
    namespaces.append(namespace)
with (output / (prefix + '.emapper.annotations')).open('w') as handle:
    handle.write('## confidence codes: h=high m=medium l=low -=not annotated\n')
    handle.write('## confidence field order: ' + ' '.join(egg.FIELDS) + '\n')
    handle.write('#' + '\t'.join(egg.HEADER) + '\n')
    writer = csv.DictWriter(handle, fieldnames=egg.HEADER, delimiter='\t', lineterminator='\n')
    writer.writerows(rows)
    handle.write(f'## {len(rows)} queries scanned\n')
with (output / (prefix + '.emapper.annotations.go_namespaces.tsv')).open('w') as handle:
    writer = csv.DictWriter(handle, fieldnames=egg.GO_HEADER, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(namespaces)
with (output / (prefix + '.emapper.orthologs')).open('w') as handle:
    handle.write('#query\torth_type\tspecies\torthologs\n')
    for protein in proteins:
        handle.write(f"{protein['query']}\tone2one\tSynthetic species(1)\t*synthetic_ortholog\n")
    handle.write(f'## {len(proteins)} queries scanned\n')
"""

CHANGE_WORKFLOW = """nextflow.enable.dsl = 2
include { ANNOTATION_RESOURCES; FUNCTIONAL_ANNOTATION } from './subworkflows/local/functional_annotation'
include { IMPORT_ANNOTATION_SOURCE } from './modules/local/import_annotation_source'
include { AGGREGATE_ANNOTATIONS; ANNOTATION_ACCEPTANCE } from './modules/local/functional_annotation'
workflow {
    current = file(params.fixture_source, checkIfExists: true)
    previous = file(params.annotation_from, checkIfExists: true)
    ANNOTATION_RESOURCES()
    IMPORT_ANNOTATION_SOURCE(Channel.value(current), ANNOTATION_RESOURCES.out.receipt)
    bundles = IMPORT_ANNOTATION_SOURCE.out.samples.splitCsv(header: true, sep: '\\t').map { sample ->
        tuple([accession: sample.accession], current.resolve("samples/${sample.accession}/annotation/bundle"))
    }
    FUNCTIONAL_ANNOTATION(IMPORT_ANNOTATION_SOURCE.out.samples, bundles,
        ANNOTATION_RESOURCES.out.receipt, Channel.value(previous))
    AGGREGATE_ANNOTATIONS(FUNCTIONAL_ANNOTATION.out.plan, IMPORT_ANNOTATION_SOURCE.out.master,
        IMPORT_ANNOTATION_SOURCE.out.sample_status, FUNCTIONAL_ANNOTATION.out.bundle_files,
        FUNCTIONAL_ANNOTATION.out.results.map { item -> item[1] }.toList(),
        FUNCTIONAL_ANNOTATION.out.native_batches.toList())
    ANNOTATION_ACCEPTANCE(AGGREGATE_ANNOTATIONS.out.acceptance, IMPORT_ANNOTATION_SOURCE.out.samples)
}
"""


@unittest.skipUnless(NEXTFLOW, "Nextflow is required for shared eggNOG integration")
class EggnogBatchIntegrationTests(unittest.TestCase):
    def setUp(self):
        # Compose the existing fixture without inheriting or rerunning its tests.
        self.fixture = workflow_fixture.FunctionalAnnotationIntegrationTests()
        self.addCleanup(self.fixture.doCleanups)
        self.fixture.setUp()
        self.root, self.project, self.source = (
            self.fixture.root,
            self.fixture.project,
            self.fixture.source,
        )
        native_fixture = normalizer_fixture.NormalizerTests()
        native_fixture.setUp()
        self.addCleanup(native_fixture.tearDown)
        native_fixture.write_eggnog(rows=[], go_rows=[])
        self.database = self.root / "eggnog"
        shutil.copytree(native_fixture.resource, self.database)
        for name in (
            "eggnog_proteins.dmnd",
            "eggnog.db.fieldpresence.bin",
            "eggnog.db.taxids.bin",
            "eggnog.taxa.db",
            "eggnog.taxa.db.traverse.pkl",
        ):
            (self.database / name).write_bytes(b"Synthetic orchestration resource\n")
        contract = dict(
            component="eggnog",
            version="synthetic-v7",
            settings={"fixture": True},
            files=file_records(self.database),
        )
        write_json(
            self.database / "annotation_resource.json",
            dict(
                schema_version=1,
                component="eggnog",
                resource_id=identity(contract),
                contract=contract,
            ),
        )
        self.ledger = self.root / "native-calls.jsonl"
        self.fail_native = self.root / "fail-native"
        self.no_hits = self.root / "no-hits.json"
        executable = self.project / "bin/emapper.py"
        executable.write_text(
            f"#!{sys.executable}\n"
            + EMAPPER.replace("__LEDGER__", repr(str(self.ledger)))
            .replace("__FAIL__", repr(str(self.fail_native)))
            .replace("__NOHITS__", repr(str(self.no_hits)))
        )
        executable.chmod(0o755)
        diamond = self.project / "bin/diamond"
        diamond.write_text("#!/bin/sh\nprintf 'Synthetic DIAMOND version control\\n'\n")
        diamond.chmod(0o755)
        (self.project / "batch_change_control.nf").write_text(CHANGE_WORKFLOW)

    def calls(self, mode="diamond"):
        if not self.ledger.exists():
            return []
        records = [json.loads(line) for line in self.ledger.read_text().splitlines()]
        return [record["queries"] for record in records if record["mode"] == mode]

    def run_pipeline(
        self,
        name,
        previous,
        *,
        current=None,
        succeeds=True,
        resume=False,
        incompatible_archives=False,
    ):
        configuration = self.root / f"{name}.config"
        configuration.write_text(
            "params {\n annotation_tools = 'eggnog'\n"
            + f" eggnog_db = '{self.database}'\n"
            + " eggnog_container = 'sha256:"
            + "1" * 64
            + "'\n"
            + " python_container = 'sha256:"
            + "2" * 64
            + "'\n"
            + " annotation_cpus = 1\n annotation_memory = 18.GB\n"
            + " max_cpus = 1\n max_memory = 18.GB\n}\n"
        )
        output = self.root / name
        entrypoint = "reannotate.nf" if current is None else "batch_change_control.nf"
        command = [
            NEXTFLOW,
            "run",
            str(self.project / entrypoint),
            "-profile",
            "test,docker",
            "-c",
            str(configuration),
            "--annotation_from",
            str(previous),
            "--outdir",
            str(output),
            "-work-dir",
            str(self.root / f"work-{name}"),
        ]
        if current is not None:
            command.extend(["--fixture_source", str(current)])
        if resume:
            command.append("-resume")
        result = subprocess.run(
            command,
            cwd=self.root,
            text=True,
            capture_output=True,
            timeout=120,
            env=self.fixture.environment,
        )
        (self.root / f"{name}.log").write_text(result.stdout + result.stderr)
        self.assertEqual(
            result.returncode == 0, succeeds, result.stdout + result.stderr
        )
        if incompatible_archives:
            self.assertIn(
                "annotation_batches are incompatible", result.stdout + result.stderr
            )
            self.assertIn("fresh --outdir", result.stdout + result.stderr)
            trace = read_tsv(output / "pipeline_info/trace.tsv")
            self.assertFalse(any(":ANNOTATION_SEARCH " in row["name"] for row in trace))
            return output
        self.assertTrue(
            (output / "annotation_results.json").is_file(),
            result.stdout + result.stderr,
        )
        manifest, _ = validate_source(output)
        self.assertEqual(manifest["complete"], succeeds)
        return output

    def assert_publication(
        self, output, accessions, action, *, succeeds=True, searches=1
    ):
        manifest, batches = validate_source(output)
        self.assertEqual(manifest["accessions"], accessions)
        self.assertEqual(len(batches), 1)
        native = next(iter(batches.values()))
        self.assertEqual(set(native.members), set(accessions))
        self.assertEqual(list((output / "annotation_batches").iterdir()), [native.root])
        self.assertEqual(
            list((output / "annotation_batches").rglob("raw")), [native.root / "raw"]
        )
        self.assertEqual(list((output / "samples").rglob("raw")), [])
        plan = [
            row
            for row in read_tsv(output / "pipeline_info/annotation_plan.tsv")
            if row["tool"] == "eggnog"
        ]
        self.assertEqual([row["accession"] for row in plan], accessions)
        self.assertEqual({row["action"] for row in plan}, {action})
        self.assertEqual(
            {row["batch_id"] for row in plan}, {native.record["batch"]["batch_id"]}
        )
        self.assertEqual(len({row["task_directory"] for row in plan}), 1)
        trace = read_tsv(output / "pipeline_info/trace.tsv")
        self.assertEqual(
            sum(":ANNOTATION_SEARCH " in row["name"] for row in trace), searches
        )
        self.assertEqual(
            sum(":COMPLETE_EGGNOG_BATCH " in row["name"] for row in trace), 1
        )
        results, original_ids, tool_ids = {}, [], []
        for accession in accessions:
            sample = output / "samples" / accession / "annotation"
            _, proteins = bundle_proteins(sample / "bundle")
            self.assertEqual(len(proteins), 1)
            protein = proteins[0]
            original_ids.append(protein["protein_id"])
            tool_ids.append(protein["tool_id"])
            result = read_json(sample / "eggnog/result.json")
            results[accession] = result
            self.assertEqual(result["status"], "success" if succeeds else "failed")
            self.assertEqual(result["action"], action)
            self.assertEqual(
                result["native_batch"]["batch_result_id"],
                native.record["batch_result_id"],
            )
            self.assertFalse((sample / "eggnog/raw").exists())
            if succeeds:
                for table in ("eggnog_annotations.tsv", "eggnog_seed_hits.tsv"):
                    rows = read_tsv(sample / "eggnog/normalized" / table)
                    self.assertEqual(
                        [(row["accession"], row["gene_id"]) for row in rows],
                        [(accession, protein["gene_id"])],
                    )
                self.assertEqual(
                    result["evidence"]["features"]["ko"],
                    {protein["gene_id"]: ["K00001"]},
                )
                self.assertEqual(
                    result["evidence"]["features"]["Preferred_name"],
                    {protein["gene_id"]: ["hflB"]},
                )
        self.assertEqual(set(original_ids), {"gene_1"})
        self.assertEqual(len(set(tool_ids)), len(accessions))
        if succeeds:
            self.assertEqual(
                read_tsv(output / "tables/functional_matrices/eggnog_ko_counts.tsv"),
                [dict(accession=accession, K00001="1") for accession in accessions],
            )
            search = native.root / "raw/search/eggnog.emapper.seed_orthologs"
            self.assertTrue(
                search.read_text().endswith(f"## {len(accessions)} queries scanned\n")
            )
            self.assertFalse(search.with_name("eggnog.emapper.annotations").exists())
            for index in range(len(accessions)):
                member = native.root / "raw/annotations" / f"member{index:08d}"
                for name in ("annotations", "orthologs"):
                    self.assertTrue(
                        (member / f"eggnog.emapper.{name}")
                        .read_text()
                        .endswith("## 1 queries scanned\n")
                    )
        return native.record, results

    def input_source(self, name, accessions, changed=None):
        source = self.root / name
        (source / "tables").mkdir(parents=True)
        bundles = []
        for accession in accessions:
            origin = (changed or {}).get(
                accession, self.source / "samples" / accession / "annotation/bundle"
            )
            destination = source / "samples" / accession / "annotation/bundle"
            shutil.copytree(origin, destination)
            bundles.append(destination)
        write_tsv(
            source / "tables/master_table.tsv",
            ("Accession", "Gcode"),
            [dict(Accession=acc, Gcode=4) for acc in accessions],
        )
        write_tsv(
            source / "tables/sample_status.tsv",
            ("accession", "gcode_status"),
            [dict(accession=acc, gcode_status="done") for acc in accessions],
        )
        write_tsv(
            source / "tables/validated_samples.tsv",
            ("accession",),
            [dict(accession=acc) for acc in accessions],
        )
        publish_disabled(source, accessions, bundles)
        return source

    def test_one_native_batch_reuses_and_renormalizes_portably(self):
        first = self.run_pipeline("initial", self.source)
        packet, results = self.assert_publication(first, ["B", "A"], "run")
        self.assertEqual(len(self.calls()), 1)
        self.assertEqual(len(set(self.calls()[0])), 2)
        annotations = self.calls("no_search")
        self.assertEqual(len(annotations), 2)
        self.assertEqual(
            sorted(query for call in annotations for query in call),
            sorted(self.calls()[0]),
        )
        self.assertTrue(all(len(call) == 1 for call in annotations))
        resumed = self.run_pipeline("initial", self.source, resume=True)
        self.assert_publication(resumed, ["B", "A"], "run")
        searches = [
            row
            for row in read_tsv(resumed / "pipeline_info/trace.tsv")
            if ":ANNOTATION_SEARCH " in row["name"]
        ]
        self.assertEqual({row["status"] for row in searches}, {"CACHED"})
        self.assertEqual(len(self.calls()), 1)
        shutil.rmtree(self.root / "work-initial")
        portable = self.root / "portable first results"
        first.rename(portable)
        reused = self.run_pipeline("reused", portable)
        reused_packet, reused_results = self.assert_publication(
            reused, ["B", "A"], "reuse", searches=0
        )
        self.assertEqual(reused_packet, packet)
        self.assertEqual(len(self.calls()), 1)
        for accession in results:
            self.assertEqual(
                reused_results[accession]["normalized_files"],
                results[accession]["normalized_files"],
            )
        shutil.rmtree(self.root / "work-reused")
        shutil.rmtree(portable)
        parser = self.project / "bin/normalize_eggnog.py"
        parser.write_text(
            parser.read_text() + "\n# Interpretation-only integration control.\n"
        )
        renormalized = self.run_pipeline("renormalized", reused)
        changed_packet, changed_results = self.assert_publication(
            renormalized, ["B", "A"], "renormalize", searches=0
        )
        self.assertEqual(changed_packet, packet)
        self.assertEqual(len(self.calls()), 1)
        self.assertEqual(self.calls("no_search"), annotations)
        for accession in results:
            self.assertEqual(
                changed_results[accession]["search_fingerprint"],
                results[accession]["search_fingerprint"],
            )
            self.assertNotEqual(
                changed_results[accession]["normalization_fingerprint"],
                results[accession]["normalization_fingerprint"],
            )
            self.assertEqual(
                changed_results[accession]["evidence"], results[accession]["evidence"]
            )
        resource_file = self.database / "annotation_resource.json"
        resource = read_json(resource_file)
        resource["contract"]["version"] = "synthetic-v7-changed-method"
        resource["resource_id"] = identity(resource["contract"])
        write_json(resource_file, resource)
        self.run_pipeline(
            "renormalized",
            reused,
            resume=True,
            succeeds=False,
            incompatible_archives=True,
        )
        self.assertEqual(len(self.calls()), 1)

    def test_membership_and_input_changes_execute_new_native_batches(self):
        first = self.run_pipeline("initial", self.source)
        packet, results = self.assert_publication(first, ["B", "A"], "run")
        subset_source = self.input_source("subset-source", ["A"])
        subset = self.run_pipeline("subset", first, current=subset_source)
        subset_packet, subset_results = self.assert_publication(subset, ["A"], "run")
        self.assertNotEqual(
            subset_packet["batch"]["batch_id"], packet["batch"]["batch_id"]
        )
        self.assertEqual(subset_results["A"]["input_id"], results["A"]["input_id"])
        self.assertNotEqual(
            subset_results["A"]["search_fingerprint"],
            results["A"]["search_fingerprint"],
        )
        fixture = bundle_fixture.AnnotationBundleTests()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        dna = "ATGTGATTTTAA"
        fixture.genome.write_text(f">original_contig\n{dna}\n")
        fixture.faa.write_text(
            ">gene_1 a description with pipes|and_underscores\nMWF\n"
        )
        fixture.gff.write_text(fixture.gff.read_text().replace("ATGTGAGCTTAA", dna))
        genbank = SeqIO.read(fixture.gbk, "genbank")
        genbank.seq = Seq(dna)
        genbank.features[0].qualifiers["translation"] = ["MWF"]
        SeqIO.write(genbank, fixture.gbk, "genbank")
        changed_bundle = fixture.build("A", name="changed-A")
        changed_source = self.input_source(
            "changed-source", ["B", "A"], {"A": changed_bundle}
        )
        changed = self.run_pipeline("changed", first, current=changed_source)
        changed_packet, changed_results = self.assert_publication(
            changed, ["B", "A"], "run"
        )
        self.assertNotEqual(
            changed_packet["batch"]["batch_id"], packet["batch"]["batch_id"]
        )
        self.assertNotEqual(changed_results["A"]["input_id"], results["A"]["input_id"])
        self.assertEqual(changed_results["B"]["input_id"], results["B"]["input_id"])
        self.assertNotEqual(
            changed_results["B"]["search_fingerprint"],
            results["B"]["search_fingerprint"],
        )
        self.assertEqual([len(call) for call in self.calls()], [2, 1, 2])

    def test_whole_batch_native_failure_retains_raw_and_fails_every_member(self):
        self.fail_native.touch()
        failed = self.run_pipeline("failed", self.source, succeeds=False)
        packet, results = self.assert_publication(
            failed, ["B", "A"], "run", succeeds=False
        )
        self.assertEqual(packet["exit_code"], 9)
        self.assertEqual({result["exit_code"] for result in results.values()}, {9})
        self.assertEqual(len({result["reason"] for result in results.values()}), 1)
        self.assertEqual(len(self.calls()), 1)
        acceptance = read_json(failed / "annotation_acceptance.json")
        self.assertFalse(acceptance["complete"])
        self.assertEqual(
            {row["accession"] for row in acceptance["failures"]}, {"A", "B"}
        )
        for row in read_tsv(failed / "tables/master_table.tsv"):
            self.assertEqual(row["eggnog_accepted_proteins"], "NA")
        for row in read_tsv(failed / "tables/protein_function_summary.tsv"):
            self.assertEqual(row["eggnog_ko"], "NA")
        native = next(iter(validate_source(failed)[1].values()))
        self.assertIn(
            "Intentional synthetic whole-batch native failure",
            (native.root / "raw/search/tool.log").read_text(),
        )
        retained = list(
            (failed / "pipeline_info/annotation_failures").rglob("raw/exit_code.txt")
        )
        self.assertEqual([path.read_text() for path in retained], ["9\n"])
        self.fail_native.unlink()
        recovered = self.run_pipeline("failed", self.source, resume=True)
        self.assert_publication(recovered, ["B", "A"], "run")
        self.assertEqual(len(self.calls()), 2)
        self.assertEqual(retained[0].read_text(), "9\n")

    def test_encoded_accessions_survive_published_source_reuse(self):
        fixture = bundle_fixture.AnnotationBundleTests()
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        source = self.root / "encoded-source"
        (source / "tables").mkdir(parents=True)
        accessions = ["Sample A", "sample:1", "sample+1", "sample%1"]
        bundles = []
        for index, accession in enumerate(accessions):
            target = source / "samples" / accession / "annotation/bundle"
            shutil.copytree(fixture.build(accession, f"encoded-{index}"), target)
            bundles.append(target)
        write_tsv(
            source / "tables/master_table.tsv",
            ("Accession", "Gcode"),
            [dict(Accession=acc, Gcode=4) for acc in accessions],
        )
        write_tsv(
            source / "tables/sample_status.tsv",
            ("accession", "gcode_status"),
            [dict(accession=acc, gcode_status="done") for acc in accessions],
        )
        write_tsv(
            source / "tables/validated_samples.tsv",
            ("accession",),
            [dict(accession=acc) for acc in accessions],
        )
        publish_disabled(source, accessions, bundles)
        first = self.run_pipeline("encoded-first", source)
        self.assert_publication(first, accessions, "run")
        shutil.rmtree(self.root / "work-encoded-first")
        reused = self.run_pipeline("encoded-reused", first)
        self.assert_publication(reused, accessions, "reuse", searches=0)
        self.assertEqual(len(self.calls()), 1)


if __name__ == "__main__":
    unittest.main()
