"""Small published v0.4 sources for workflow and cohort validation tests."""

import shutil
from pathlib import Path

from aggregate_annotations import aggregate
from annotation_common import TOOLS, digest, identity, read_json, write_json


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
