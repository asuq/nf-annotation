#!/usr/bin/env perl
use strict;
use warnings;
use JSON::PP;
use Digest::MD5 qw(md5_hex);
use Encode qw(decode FB_CROAK LEAVE_SRC);

sub slurp {
    my ($path) = @_;
    open my $f, '<:raw', $path or die "Cannot read $path: $!\n";
    local $/; return <$f>;
}
sub emit_json {
    my ($path, $value) = @_;
    open my $f, '>:raw', $path or die "Cannot write $path: $!\n";
    print {$f} JSON::PP->new->ascii->canonical->pretty->encode($value);
    close $f or die "Cannot close $path: $!\n";
}
sub records {
    my ($bytes) = @_;
    my (@records, %seen);
    for my $line (split /(?<=\n)/, $bytes) {
        if ($line =~ /^>/) {
            my $header = $line; $header =~ s/^>//; $header =~ s/\r?\n$//;
            my $text = decode('UTF-8', $header, FB_CROAK | LEAVE_SRC);
            my ($id, $description) = $text =~ /^(\S+)(?:[ \t]+(.*))?$/s;
            die "Invalid or duplicate FASTA identifier\n" if !defined($id) || $id =~ /[\x00-\x1f\x7f]/ || $seen{$id}++;
            push @records, {source_id=>$id, source_header=>$text, description=>$description // '', raw_sequence=>'', sequence=>''};
        } else {
            next if !@records && $line =~ /^\s*$/;
            die "Sequence before FASTA header\n" unless @records;
            my $sequence = $line; $sequence =~ s/\r?\n$//;
            die "Invalid DNA or sequence whitespace\n" unless $sequence =~ /^[ACGTRYSWKMBDHVNacgtryswkmbdhvn]*$/;
            $records[-1]{raw_sequence} .= $line;
            $records[-1]{sequence} .= $sequence;
        }
    }
    die "Empty FASTA or contig\n" if !@records || grep { !length($_->{sequence}) } @records;
    return \@records;
}
sub ledger {
    my ($bytes) = @_;
    my $records = records($bytes);
    my (@mapping, $native);
    $native = '';
    for my $i (0 .. $#$records) {
        my $r = $records->[$i];
        my $id = sprintf('ctg%08d', $i+1);
        my $entry = '>' . $id . "\n" . $r->{raw_sequence};
        # A final sequence without newline is retained exactly; a subsequent
        # header necessarily had a preceding newline in the source.
        $native .= $entry;
        push @mapping, {tool_id=>$id, source_id=>$r->{source_id}, source_header=>$r->{source_header},
            description=>$r->{description}, length=>length($r->{sequence}),
            sequence_md5=>md5_hex($r->{sequence}), sequence_lines_md5=>md5_hex($r->{raw_sequence})};
    }
    return ({schema=>'ccfinder-contig-map-v1', source_fasta_md5=>md5_hex($bytes),
             tool_fasta_md5=>md5_hex($native), contigs=>\@mapping}, $native, $records);
}
sub restore_values {
    my ($value, $aliases, $field) = @_;
    if (ref($value) eq 'ARRAY') { return [map { restore_values($_, $aliases, $field) } @$value]; }
    if (ref($value) eq 'HASH') { return {map { $_ => restore_values($value->{$_}, $aliases, $_) } keys %$value}; }
    # Numeric and biological sequence fields retain their native JSON types and
    # bytes. Only identifier-bearing strings can contain our generated aliases.
    if (!ref($value) && defined($value) && defined($field)
            && $field =~ /^(?:Id|Name|Version|Gene|Sequence|Seqid|seq_id|contig|Summary_CRISPR-Cas)$/) {
        $value =~ s/(?<![A-Za-z0-9])\b(ctg[0-9]{8})(?=_|\b)/
            exists($aliases->{$1}) ? $aliases->{$1} : die("Foreign native contig reference\n")/ge;
    }
    return $value;
}
my $mode = shift @ARGV // '';
if ($mode eq 'stage') {
    my ($input, $output, $mapping, $contig_dir, $list) = @ARGV;
    die "stage requires input output mapping contig-directory list\n" unless @ARGV == 5;
    die "Staging outputs already exist\n" if -e $output || -e $mapping || -e $list;
    my ($record, $native, $source) = ledger(slurp($input));
    open my $out, '>:raw', $output or die "$output: $!\n"; print {$out} $native; close $out or die "$!\n";
    open my $names, '>:raw', $list or die "$list: $!\n";
    for my $i (0 .. $#{$record->{contigs}}) {
        my $id = $record->{contigs}[$i]{tool_id};
        my $path = "$contig_dir/$id.fna";
        die "Contig staging collision\n" if -e $path;
        open my $f, '>:raw', $path or die "$path: $!\n";
        print {$f} ">$id\n", $source->[$i]{raw_sequence}; close $f or die "$!\n";
        print {$names} "$id.fna\n";
    }
    close $names or die "$!\n";
    emit_json($mapping, $record);
} elsif ($mode eq 'restore') {
    my ($input, $native_input, $mapping, $native_result, $output) = @ARGV;
    die "restore requires input tool-input mapping native-result output\n" unless @ARGV == 5;
    my ($expected, $native) = ledger(slurp($input));
    my $record = decode_json(slurp($mapping));
    my $encoder = JSON::PP->new->ascii->canonical;
    die "Contig mapping or staged sequence bytes changed\n" unless $encoder->encode($record) eq $encoder->encode($expected) && slurp($native_input) eq $native;
    my $doc = decode_json(slurp($native_result));
    die "Invalid native result structure\n" unless ref($doc) eq 'HASH' && ref($doc->{Sequences}) eq 'ARRAY';
    my %mapping = map { $_->{tool_id} => $_ } @{$record->{contigs}};
    my (%seen, %restored);
    for my $seq (@{$doc->{Sequences}}) {
        die "Malformed native sequence result\n" unless ref($seq) eq 'HASH';
        my $id = $seq->{Id} // '';
        my $entry = $mapping{$id};
        # Regex coercion on the original scalar can change JSON::PP's numeric
        # encoding. Inspect a copy so the retained native Length stays numeric.
        my $length = $seq->{Length};
        die "Unknown, duplicate or changed native contig\n" unless $entry && !$seen{$id}++ && defined($length) && $length =~ /^\d+$/ && $length == $entry->{length};
        die "Missing native CRISPR/Cas arrays\n" unless ref($seq->{Crisprs}) eq 'ARRAY' && ref($seq->{Cas}) eq 'ARRAY';
        my $restored = restore_values($seq, {$id=>$entry->{source_id}});
        $restored->{Id} = $entry->{source_id};
        $restored->{Version} = $entry->{source_id};
        $restored->{Description} = $entry->{description};
        $restored{$id} = $restored;
    }
    die "Native result omits a staged contig\n" unless keys(%seen) == keys(%mapping);
    $doc->{Sequences} = [map { $restored{$_->{tool_id}} } @{$record->{contigs}}];
    emit_json($output, $doc);
} else { die "Expected stage or restore\n"; }
