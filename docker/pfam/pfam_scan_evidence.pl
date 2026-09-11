#!/usr/bin/env perl
use strict;
use warnings;
use Bio::Pfam::Scan::PfamScan;
use Cwd qw(abs_path);
use Getopt::Long qw(GetOptions);
use JSON;

my ($fasta, $database, $outdir, $cpus);
GetOptions('fasta=s' => \$fasta, 'database=s' => \$database,
           'outdir=s' => \$outdir, 'cpus=i' => \$cpus) or die "Invalid options\n";
die "Require --fasta FILE --database DIR --outdir NEW_DIR --cpus N\n"
  unless defined($fasta) && -s $fasta && defined($database) && -d $database
     && defined($outdir) && !-e $outdir && defined($cpus) && $cpus > 0 && !@ARGV;
mkdir $outdir or die "Cannot create $outdir: $!\n";
$outdir = abs_path($outdir);
local $ENV{NF_PFAM_EVIDENCE_DIR} = $outdir;

# Run native gathering-threshold search once, retaining all clan alternatives.
my $scan = Bio::Pfam::Scan::PfamScan->new(
    -fasta => abs_path($fasta), -dir => abs_path($database),
    -hmmlib => ['Pfam-A.hmm'], -version => '1.6', -cpu => $cpus,
    -clan_overlap => 1,
);
$scan->search;
$scan->write_results("$outdir/pfam.raw.tsv");

# Apply the pinned upstream resolver to those same results. It uses domain
# coordinates, E-values, clans and nesting; significance flags do not enter it.
$scan->_resolve_clan_overlap;
for my $line (@{$scan->{_header}}) {
    $line =~ s/resolve clan overlaps: off/resolve clan overlaps: on/;
}
$scan->write_results("$outdir/pfam.resolved.tsv");
open(my $json, '>', "$outdir/pfam.resolved.json") or die "Cannot write JSON: $!\n";
print {$json} JSON->new->canonical->encode($scan->results), "\n";
close($json) or die "Cannot close JSON: $!\n";
