"""Container contract tests for the CRISPRCasFinder runtime image."""

from __future__ import annotations

import os
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCKER_DIR = ROOT / "docker" / "ccfinder"
DOCKERFILE = DOCKER_DIR / "Dockerfile"
PATCH_FILE = DOCKER_DIR / "CRISPRCasFinder.runtime.patch"
SMOKE_SCRIPT = DOCKER_DIR / "ccfinder-smoke-test.sh"
WRAPPER_SCRIPT = DOCKER_DIR / "macsyfinder-wrapper.sh"
MACSYFINDER_TARBALL = DOCKER_DIR / "macsyfinder-1.0.5.tar.gz"
IMAGE_TAG = "codex-ccfinder-test:runtime-contract"
RUN_DOCKER_TESTS = os.environ.get("RUN_DOCKER_TESTS") == "1"


def run_command(args: list[str]) -> subprocess.CompletedProcess[str]:
    """Run one subprocess and return the completed process."""
    return subprocess.run(
        args,
        check=False,
        capture_output=True,
        text=True,
        cwd=ROOT,
    )


class CCFinderContainerContractTestCase(unittest.TestCase):
    """Lock the Dockerfile to the inspected Singularity runtime contract."""

    def test_dockerfile_records_inspected_singularity_runtime(self) -> None:
        """Require the Dockerfile to pin the inspected Singularity image reference."""
        dockerfile_text = DOCKERFILE.read_text(encoding="utf-8")

        self.assertIn(
            "quay.io/singularity/singularity:v3.8.2-slim@sha256:2558c93c63aa8f721aae106ff1052dff5442171d7c00aed3e81d0fe62172b118",
            dockerfile_text,
        )
        self.assertIn("ARG CCFINDER_VERSION=4.2.30", dockerfile_text)
        self.assertIn("ARG CCFINDER_REF=release-4.2.20", dockerfile_text)
        self.assertIn("ARG MACSYFINDER_VERSION=1.0.5", dockerfile_text)
        self.assertIn(
            'org.opencontainers.image.base.name="ubuntu:18.04"', dockerfile_text
        )
        self.assertIn(
            'io.nf_myco_annotation.reference_runtime="CRISPRCasFinder.simg"',
            dockerfile_text,
        )
        self.assertIn(
            'io.nf_myco_annotation.reference_runtime_version="4.2.30"', dockerfile_text
        )
        self.assertIn(
            "COPY macsyfinder-1.0.5.tar.gz /tmp/macsyfinder-1.0.5.tar.gz",
            dockerfile_text,
        )
        self.assertIn(
            "cpanm --notest \\\n        Bio::DB::Fasta \\\n        Bio::Seq \\\n        Bio::SeqIO \\\n        Bio::AlignIO;",
            dockerfile_text,
        )
        self.assertNotIn(
            'ensure_perl_module "Bio::AlignIO" "Bio::AlignIO";', dockerfile_text
        )
        self.assertNotIn(
            'ensure_perl_module "Bio::SeqIO" "Bio::SeqIO";', dockerfile_text
        )
        self.assertNotIn(
            'ensure_perl_module "Bio::DB::Fasta" "Bio::DB::Fasta";', dockerfile_text
        )
        self.assertIn(
            "COPY --from=builder /usr/local/lib/x86_64-linux-gnu/perl/5.26.1 /usr/local/lib/x86_64-linux-gnu/perl/5.26.1",
            dockerfile_text,
        )
        self.assertIn(
            "COPY --from=builder /usr/local/share/macsyfinder /usr/local/share/macsyfinder",
            dockerfile_text,
        )
        self.assertIn(
            "COPY --from=builder /usr/local/share/perl/5.26.1 /usr/local/share/perl/5.26.1",
            dockerfile_text,
        )
        self.assertIn(
            "COPY --from=builder /etc/macsyfinder /etc/macsyfinder", dockerfile_text
        )

    def test_runtime_patch_matches_the_inspected_runtime_contract(self) -> None:
        """Require the runtime patch to encode the inspected legacy runtime behaviour."""
        patch_text = PATCH_FILE.read_text(encoding="utf-8")
        wrapper_script = WRAPPER_SCRIPT.read_text(encoding="utf-8")

        self.assertIn('my $version = "4.2.30";', patch_text)
        self.assertIn(
            'my $PathVmatch2 = "/usr/local/CRISPRCasFinder/bin/vmatch2";', patch_text
        )
        self.assertIn(
            'my $PathMkvtree2 = "/usr/local/CRISPRCasFinder/bin/mkvtree2";', patch_text
        )
        self.assertIn(
            'my $PathVsubseqselect2 = "/usr/local/CRISPRCasFinder/bin/vsubseqselect2";',
            patch_text,
        )
        self.assertIn(
            'my $PathMacsyfinder = "/usr/local/CRISPRCasFinder/macsyfinder";',
            patch_text,
        )
        self.assertIn('my $casfinder = "CasFinder-2.0.4";', patch_text)
        self.assertIn('$casdb = $casfinder."/DEF-Class-2.0.4/";', patch_text)
        self.assertIn('$casdb = $casfinder."/DEF-Typing-2.0.4/";', patch_text)
        self.assertIn('$casdb = $casfinder."/DEF-SubTyping-2.0.4/";', patch_text)
        self.assertIn('my $profiles = $casfinder."/CASprofiles-2.0.4/";', patch_text)
        self.assertIn(
            "makesystemcall(\"$PathVmatch2 \" . join(' ',@vmatchoptions));", patch_text
        )
        self.assertIn(
            "makesystemcall(\"$PathVsubseqselect2 \" . join(' ',@subselectoptions));",
            patch_text,
        )
        self.assertIn('makesystemcall("$PathMkvtree2 -db $inputfile " .', patch_text)
        self.assertIn('$macsyfinder = "$PathMacsyfinder -w $cpuMacSyFinder', patch_text)
        self.assertIn("$found= 1; # PA", patch_text)
        self.assertGreater(MACSYFINDER_TARBALL.stat().st_size, 10_000_000)
        self.assertIn('exec /usr/local/bin/macsyfinder "$@"', wrapper_script)

    def test_smoke_script_uses_module_style_explicit_resource_paths(self) -> None:
        """Require the smoke harness to exercise the clean module-style output contract."""
        smoke_script = SMOKE_SCRIPT.read_text(encoding="utf-8")

        self.assertIn(
            "/usr/local/CRISPRCasFinder/install_test/sequence.fasta", smoke_script
        )
        self.assertIn('run_root="${task_root}/ccfinder_run"', smoke_script)
        self.assertIn('tool_output_root="${task_root}/ccfinder_raw"', smoke_script)
        self.assertIn(": > index.html", smoke_script)
        self.assertIn('-outdir "${tool_output_root}"', smoke_script)
        self.assertIn("-soFile /usr/local/CRISPRCasFinder/sel392v2.so", smoke_script)
        self.assertIn(
            "-DBcrispr /usr/local/CRISPRCasFinder/supplementary_files/CRISPR_crisprdb.csv",
            smoke_script,
        )
        self.assertNotIn("set +e", smoke_script)
        self.assertNotIn("tool_exit_code=$?", smoke_script)
        self.assertNotIn("cannot move '.*' to a subdirectory of itself", smoke_script)
        self.assertIn(
            "find \"${tool_output_root}\" -type f -name 'Casfinder_summary_*.tsv'",
            smoke_script,
        )
        self.assertIn(
            "find \"${tool_output_root}\" -type f -name 'rawCas.fna'", smoke_script
        )
        self.assertIn(
            "find \"${tool_output_root}\" -type f -name 'result.json'", smoke_script
        )

    @unittest.skipUnless(
        RUN_DOCKER_TESTS, "Set RUN_DOCKER_TESTS=1 to build the container contract test."
    )
    def test_built_container_passes_the_smoke_contract(self) -> None:
        """Build the image and require the patched legacy runtime contract to pass the smoke test."""
        build_result = run_command(
            [
                "docker",
                "build",
                "--platform",
                "linux/amd64",
                "-t",
                IMAGE_TAG,
                str(DOCKER_DIR),
            ]
        )
        self.assertEqual(
            build_result.returncode,
            0,
            msg=f"Docker build failed.\nSTDOUT:\n{build_result.stdout}\nSTDERR:\n{build_result.stderr}",
        )

        run_result = run_command(
            [
                "docker",
                "run",
                "--rm",
                "--platform",
                "linux/amd64",
                IMAGE_TAG,
                "bash",
                "-lc",
                (
                    "set -euo pipefail; "
                    "grep -F 'my $version = \"4.2.30\";' /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl; "
                    "grep -F 'my $PathMacsyfinder = \"/usr/local/CRISPRCasFinder/macsyfinder\";' "
                    "/usr/local/CRISPRCasFinder/CRISPRCasFinder.pl; "
                    "test -f /usr/local/share/perl/5.26.1/Bio/DB/IndexedBase.pm; "
                    "test -f /usr/local/share/perl/5.26.1/Bio/DB/Fasta.pm; "
                    'test -z "$(perl -MBio::DB::IndexedBase -e 1 2>&1)"; '
                    "test -d /usr/local/share/macsyfinder/DEF; "
                    "test -f /etc/macsyfinder/macsyfinder.conf.new; "
                    "/usr/local/bin/ccfinder-smoke-test"
                ),
            ]
        )
        self.assertEqual(
            run_result.returncode,
            0,
            msg=f"Runtime contract failed.\nSTDOUT:\n{run_result.stdout}\nSTDERR:\n{run_result.stderr}",
        )
