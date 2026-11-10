"""Exercise rendered native wrapper commands, including soft failures."""

import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which("nextflow"), "Nextflow is required")
class CcfinderCommandTests(unittest.TestCase):
    def test_uniform_staging_and_failure_preservation(self):
        with tempfile.TemporaryDirectory(
            prefix="ccfinder-commands-", dir="/tmp"
        ) as name:
            root = Path(name)
            shutil.copytree(ROOT / "bin", root / "bin")
            module = (
                (ROOT / "modules/local/ccfinder.nf")
                .read_text()
                .replace("'/usr/local/CRISPRCasFinder'", repr(str(root / "native")))
            )
            (root / "ccfinder.nf").write_text(module)
            (root / "native").mkdir()
            (root / "native/CRISPRCasFinder.pl").write_text(r"""
use strict; use warnings; use JSON::PP;
if (@ARGV == 1 && $ARGV[0] eq '-v') { print "version fixture, test\n"; exit 0; }
my %args; for(my $i=0;$i<@ARGV;$i++){ $args{$ARGV[$i]}=$ARGV[$i+1] if $ARGV[$i] =~ /^-/; }
die 'genetic code' unless $args{'-gcode'} =~ /^(4|11)$/;
die 'threads' unless $args{'-cpuMacSyFinder'} eq '2' && $args{'-cpuProkka'} eq '2';
die 'scientific flags' unless exists($args{'-cas'}) && exists($args{'-ccvRep'}) && $args{'-levelMin'} eq '3';
open my $f,'<',$args{'-in'} or die $!; my @rows;
while(<$f>){ chomp; if(/^>(\S+)$/){ push @rows,{Id=>$1,Length=>0,Crisprs=>[],Cas=>[]}; }
elsif(/^>/){die 'unsafe input header';} else {$rows[-1]{Length}+=length($_);} }
die 'parent contigs' unless -f '../ctg00000001.fna';
mkdir $args{'-outdir'}; open my $out,'>',"$args{'-outdir'}/result.json" or die $!;
my $len=$rows[0]{Length};
if($len==10){ print {$out} encode_json({Sequences=>\@rows}); exit 7; }
if($len==20){ print {$out} '{invalid'; exit 0; }
if($len==30){pop @rows;}
print {$out} encode_json({Sequences=>\@rows});
""")
            samples = []
            for length, code in ((10, 4), (20, 11), (30, 4), (40, 11)):
                path = root / f"S{length}.fa"
                path.write_text(
                    f">same.1 description\twith tabs\n{'A' * length}\n>same.2 second\n{'C' * length}\n"
                )
                samples.append(f"tuple([accession:'S{length}'],file('{path}'),{code})")
            (root / "main.nf").write_text(
                "include { CCFINDER } from './ccfinder'\nworkflow { CCFINDER(Channel.of("
                + ",".join(samples)
                + ")) }\n"
            )
            (root / "nextflow.config").write_text(
                f"params {{ outdir='{root}/results'; soft_fail_attempts=1; ccfinder_extra_args=''; annotation_cleanup='success' }}\nprocess.cpus=2\n"
            )
            result = subprocess.run(
                [shutil.which("nextflow"), "run", "main.nf", "-ansi-log", "false"],
                cwd=root,
                env=dict(os.environ, NXF_DISABLE_CHECK_LATEST="true"),
                capture_output=True,
                text=True,
                timeout=180, check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            for length in (10, 20, 30):
                sample = root / f"results/samples/S{length}/ccfinder"
                self.assertEqual((sample / "result.json").stat().st_size, 0)
                self.assertTrue((sample / "native_result.json").stat().st_size)
            self.assertIn(
                "exit_code=7",
                (root / "results/samples/S10/ccfinder/ccfinder.log").read_text(),
            )
            sample = root / "results/samples/S40/ccfinder"
            result = json.loads((sample / "result.json").read_text())
            self.assertEqual(
                [row["Id"] for row in result["Sequences"]], ["same.1", "same.2"]
            )
            self.assertEqual(
                result["Sequences"][0]["Description"], "description\twith tabs"
            )
            self.assertTrue((sample / "ccfinder_storage.json").is_file())
