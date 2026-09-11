#!/usr/bin/env python3
"""Apply a fail-closed output-only adaptation to one pinned eggNOG source."""

import hashlib
import json
import shutil
from pathlib import Path

import eggnogmapper


def replace_once(source: str, old: str, new: str) -> str:
    if source.count(old) != 1:
        raise ValueError(f"Expected one source anchor: {old!r}")
    return source.replace(old, new, 1)


def main() -> None:
    package = Path(eggnogmapper.__file__).parent
    provenance = Path("/opt/nf-annotation/provenance")
    provenance.mkdir(parents=True, exist_ok=True)
    paths = {
        "annotator/e7/annotate.py": "dd084c836d2751a71c5452853971b1002108a8b9b3f01b966edef498a74bb641",
        "annotation/output.py": "e6e88ed0c131685f6fd319d6a6740a70509f4c2689f38ad4ebd68425f52d6072",
    }
    sources = {}
    for relative, expected in paths.items():
        path = package / relative
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"Unexpected eggNOG source bytes: {relative}")
        sources[relative] = path.read_text()

    relative = "annotator/e7/annotate.py"
    source = replace_once(
        sources[relative],
        "import array\n",
        "from eggnogmapper.go_namespace_export import record_namespaces\n\nimport array\n",
    )
    source = replace_once(
        source,
        "        eng = self.engine\n        merged: Set[str] = set()",
        "        eng = self.engine\n"
        "        record_namespaces(self.annotations, self._go_values, self._go_tier, eng.TIER_CONFIDENCE)\n"
        "        merged: Set[str] = set()",
    )
    source = replace_once(
        source,
        "        merged_values: Dict[str, Set[str]] = {}",
        "        namespace_values: Dict[str, Set[str]] = {}\n"
        "        namespace_tiers: Dict[str, int] = {}\n"
        "        merged_values: Dict[str, Set[str]] = {}",
    )
    source = replace_once(
        source,
        "                    bucket = merged_values.setdefault(output_field, set())",
        "                    namespace_values.setdefault(field, set()).update(values)\n"
        "                    namespace_tiers[field] = min(namespace_tiers.get(field, prio_key[2]), prio_key[2])\n"
        "                    bucket = merged_values.setdefault(output_field, set())",
    )
    source = replace_once(
        source,
        "        # Materialise multi-source merged outputs (currently just GOs).",
        "        record_namespaces(annotations, namespace_values, namespace_tiers, self.TIER_CONFIDENCE)\n\n"
        "        # Materialise multi-source merged outputs (currently just GOs).",
    )
    source = replace_once(
        source,
        "            # Merge GO sub-namespaces in _GO_NS_FIELDS order (matches eager",
        "            record_namespaces(annotations, go_values, go_tier, self.TIER_CONFIDENCE)\n"
        "            # Merge GO sub-namespaces in _GO_NS_FIELDS order (matches eager",
    )
    changed = {relative: source}

    relative = "annotation/output.py"
    source = replace_once(
        sources[relative],
        "import os\n",
        "from eggnogmapper.go_namespace_export import HEADER as GO_HEADER, write_row as write_go_row\n\nimport os\n",
    )
    source = replace_once(
        source,
        "                       applied_filters=None):\n\n    if resume == True:",
        "                       applied_filters=None):\n\n"
        "    if resume:\n"
        "        raise ValueError('GO export requires a fresh tool run; use Nextflow task resume')\n\n"
        "    if resume == True:",
    )
    source = replace_once(
        source,
        '    with open(annot_file, file_mode, encoding="utf-8", newline="\\n") as ANNOTATIONS_OUT:',
        '    with open(annot_file + ".go_namespaces.tsv", "w", encoding="utf-8", newline="\\n") as GO_OUT, open(annot_file, file_mode, encoding="utf-8", newline="\\n") as ANNOTATIONS_OUT:\n'
        '        print("\\t".join(GO_HEADER), file=GO_OUT)',
    )
    source = replace_once(
        source,
        "                                       eggnog_db=eggnog_db)\n\n            yield",
        "                                       eggnog_db=eggnog_db)\n"
        "                write_go_row(GO_OUT, annotation)\n\n            yield",
    )
    changed[relative] = source
    for relative, source in changed.items():
        compile(source, relative, "exec")
    # Preserve exact upstream code for the build's independent equivalence checks.
    for relative, source in changed.items():
        (provenance / Path(relative).name).write_text(sources[relative])
        (package / relative).write_text(source)
    helper = Path(__file__).with_name("go_namespace_export.py")
    shutil.copyfile(helper, package / helper.name)
    record = {
        "upstream_commit": "b3757a6d226047527a729546e58ff530d76a5d7d",
        "export_schema": 1,
        "original_sha256": paths,
        "patched_sha256": {
            relative: hashlib.sha256(source.encode()).hexdigest()
            for relative, source in changed.items()
        },
        "export_helper_sha256": hashlib.sha256(helper.read_bytes()).hexdigest(),
        "installer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (provenance / "go_export.json").write_text(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    main()
