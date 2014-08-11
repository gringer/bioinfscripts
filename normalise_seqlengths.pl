#!/usr/bin/perl

use warnings;
use strict;

use Pod::Usage; ## uses pod documentation in usage code
use Getopt::Long qw(:config auto_version auto_help);
use POSIX; ## for ceil

our $VERSION = "0.1";
our $DEBUG = 0;

=head1 NAME

normalise_seqlengths.pl -- converts all sequences in a FASTA file to
the same length

=head1 SYNOPSIS

./normalise_seqlengths.pl <reads.fa> [options]

=head2 Options

=over 2

=item B<-help>

Only display this help message

=item B<-fraglength>

Target fragment length (in base-pairs, default 2000)

=item B<-overlap>

Minimum overlap length (in base-pairs, default 200)

=back

=head1 DESCRIPTION

Converts all sequences in the input FASTA file to the same
length. Sequences shorter than the target length are dropped, and
sequences longer than the target length are split into overlapping
subsequences covering the entire range. This prepares the sequences
for use in an overlap-consensus assembler requiring constant-length
sequences (such as edena).

=head1 METHODS

=cut

=head2 processSeq(id, seq, length, minoverlap)

Convert the sequence I<seq> into a set of subsequences (possibly none)
that are exactly I<length> bp, and overlap at least by I<minoverlap>
bp. The converted sequences are written to stdout in FASTA format.

=cut

sub processSeq{
  my ($id, $seq, $lf, $lo) = @_;
  if(!$id || (length($seq) < $lf)){
    return;
  }
  if(length($seq) == $lf){
    $seq =~ s/(.{70})/$1\n/g;
    $seq =~ s/\s+$//;
    print(">${id}\n${seq}\n");
  }
  my $ls = length($seq);
  # number of full or partial fragments in sequence (minimum number of splits)
  my $nf = ceil($ls / $lf);
  # total overlap if the minimum number of fragments were used
  my $dl = ($nf * $lf) - $ls;
  # additional fragments required, derived from...
  my $k = ceil(($lo * ($nf - 1) - $dl) / ($lf - $lo));
  # ... mean_overlap = ($dl + $k * $lf) / ($nf + $k - 1)
  my $mo = ($dl + $k * $lf) / ($nf + $k - 1);
  $nf += $k;
  for(my $i = 0; $i < $nf; $i++){
    my $ss = int($i * ($lf-$mo));
    my $se = $ss + $lf;
    my $subSeq = substr($seq, $ss, $lf);
    if($i == ($nf - 1)){
      $ss = $ls - $lf;
      $se = $ls;
      $subSeq = substr($seq, $ls - $lf);
    }
    my $subID = $id;
    $subID =~ s/^(.*?)(\s|$)/$1#$i$2/;
    if($i == 0){
      $subID =~ s/^(.*?)(\s|$)/$1 \[$ls bp, $nf fragments, $mo overlap\]$2/;
    } else {
      $subID =~ s/^(.*?)(\s|$)/$1 \[$ss..$se\]$2/;
    }
    $subSeq =~ s/(.{70})/$1\n/g;
    $subSeq =~ s/\s+$//;
    print(">${subID}\n${subSeq}\n");
  }
}

####################################################
# Command line parsing and verification starts here
####################################################

my $argLine = join(" ",@ARGV);

my $options =
  {
   "fraglength" => 2000,
   "overlap" => 200,
};

GetOptions($options,
           'fraglength|f=i',
           'overlap|o=i',
           'debug!' => \$DEBUG,
) or pod2usage(1);

my $seqID = "";
my $seq = "";
while(<>){
  chomp;
  chomp;
  if(/^>(.+)$/){
    processSeq($seqID, $seq, $options->{"fraglength"},
               $options->{"overlap"});
    $seqID = $1;
    $seq = "";
  } else {
    $seq .= $_;
  }
}
processSeq($seqID, $seq, $options->{"fraglength"},
           $options->{"overlap"});
