"""Small published v0.4 sources for workflow and cohort validation tests."""

import shutil
import sys
from pathlib import Path

from aggregate_annotations import aggregate
from annotation_common import TOOLS, digest, identity, read_json, write_json


def install_synthetic_container_engine(project: Path) -> None:
    """Mock Docker only inside a temporary integration harness, never production."""
    script = project / "bin/docker"
    script.write_text(
        f"#!{sys.executable}\n"
        + r"""
import json
import os
import sys
from pathlib import Path

args = sys.argv[1:]
with Path(__file__).with_name('synthetic-container-calls.jsonl').open('a') as handle:
    handle.write(json.dumps(args) + '\n')
if args[0] in ('rm', 'stop', 'kill', 'pull'):
    sys.exit(0)
if args[0] == 'inspect':
    print(args[-1])
    sys.exit(0)
if args[0] != 'run':
    sys.exit('Unexpected synthetic Docker command: ' + repr(args))
i = 1
entrypoint = None
while args[i].startswith('-'):
    flag = args[i]
    if flag in ('-i', '-t', '--rm', '--privileged') or '=' in flag:
        i += 1
        continue
    if flag not in ('-v', '-w', '-u', '-e', '--env', '--name', '--cpu-shares', '--cpus', '--memory', '--memory-swap', '--entrypoint', '--user', '--platform'):
        sys.exit('Unexpected synthetic Docker option: ' + flag)
    value = args[i + 1]
    if flag in ('-e', '--env') and '=' in value:
        key, value = value.split('=', 1)
        os.environ[key] = value
    if flag == '--entrypoint':
        entrypoint = value
    if flag == '-w':
        os.chdir(value)
    i += 2
command = args[i + 1:]
if entrypoint:
    command.insert(0, entrypoint)
os.execvpe(command[0], command, os.environ)
"""
    )
    script.chmod(0o755)


def hold_bundle_until_search(project: Path, marker: Path, accession: str) -> None:
    """Withhold a bundle until an independent native search writes its marker."""
    entrypoint = project / 'reannotate.nf'
    text = entrypoint.read_text()
    delay = '''process HOLD_SLOW_BUNDLE {
    label 'process_single'
    container params.python_container
    input:
    tuple val(meta), path(bundle)
    output:
    tuple val(meta), path('released')
    script:
    """
    for attempt in \\$(seq 1 60); do
        test -f 'MARKER' && break
        sleep 1
    done
    test -f 'MARKER'
    cp -R bundle released
    """
}

'''.replace('MARKER', str(marker))
    branch = f'''    HOLD_SLOW_BUNDLE(bundles.filter {{ item -> item[0].accession == {accession!r} }})
    ready_bundles = bundles.filter {{ item -> item[0].accession != {accession!r} }}.mix(HOLD_SLOW_BUNDLE.out)
    FUNCTIONAL_ANNOTATION('''
    if text.count('IMPORT_ANNOTATION_SOURCE.out.samples, bundles,') != 1:
        raise AssertionError('Streaming fixture entrypoint changed')
    text = text.replace('workflow {', delay + 'workflow {', 1)
    text = text.replace('    FUNCTIONAL_ANNOTATION(', branch, 1)
    text = text.replace('IMPORT_ANNOTATION_SOURCE.out.samples, bundles,',
                        'IMPORT_ANNOTATION_SOURCE.out.samples, ready_bundles,', 1)
    entrypoint.write_text(text)


def refresh_tables(source: Path) -> None:
    """Re-sign an intentionally edited fixture table before semantic validation."""
    path = source / "annotation_results.json"
    record = read_json(path)
    record["tables"] = {name: digest(source / name) for name in record["tables"]}
    write_json(path, record)


def publish_disabled(source: Path, accessions: list[str], bundles=()) -> None:
    """Build the actual output contract with all functional tools disabled."""
    plan = dict(
        schema_version=1,
        accessions=accessions,
        enabled_tools=[],
        tasks=[
            dict(
                accession=accession,
                tool=tool,
                action="skip",
                input_proteins=None,
                status="skipped_disabled",
                reason="tool_disabled",
            )
            for accession in accessions
            for tool in TOOLS
        ],
    )
    plan["plan_id"] = identity(plan)
    plan_path = source / "fixture_plan.json"
    write_json(plan_path, plan)
    output = source / "fixture_aggregation"
    aggregate(
        plan_path,
        source / "tables/master_table.tsv",
        source / "tables/sample_status.tsv",
        list(bundles),
        [],
        output,
    )
    shutil.copytree(output, source, dirs_exist_ok=True)
    shutil.rmtree(output)
    plan_path.unlink()
