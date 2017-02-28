#!/usr/bin/perl

use warnings;
use strict;

use Pod::Usage; ## uses pod documentation in usage code
use Getopt::Long qw(:config auto_version auto_help);
use POSIX; ## for ceil

our $VERSION = "0.1";
our $DEBUG = 0;

=head1 NAME

fastx-repeatFilter.pl -- Repeat Match Algorithm for Long-read Evaluation

=head1 SYNOPSIS

./fastx-repeatFilter.pl <reads.fq> [options]

=head2 Options

=over 2

=item B<-help>

Only display this help message

=item B<-length>

Target repeat length

=item B<-skew>

Maximum skew between repeats for matching lines

=item B<-threshold>

Threshold score for identifying repeat reads

=item B<-fraction>

Fraction of read that must pass threshold score

=item B<-reverse>

Remove repeat-containing reads instead of extracting them

=back

=head1 DESCRIPTION

Identifies (filters) repeat-containing reads in a long-read dataset.

=head1 METHODS

something

=cut

=head2 processSeq(id, seq, qual, len, skew, threshold, fraction, reverse)

Analyses the sequence I<seq> to work out if it is likely to contain a
substantial proportion of repeats of length I<len>. Depending on the
value of I<reverse> and whether a repeat read was detected, either
return an empty string or the sequence.

=cut

sub processSeq{
  my ($id, $seq, $qual, $len, $skew, $threshold, $fraction, $reverse) = @_;
  if(!$id || (length($seq) < $len)){
    if($reverse){
      if($qual){
        print("@${id}\n${seq}\n+\n${qual}\n");
      } else {
        print(">${id}\n${seq}\n");
      }
    }
  }
}


####################################################
# Command line parsing and verification starts here
####################################################

my $argLine = join(" ",@ARGV);

my $options =
  {
   "length" => 100,
   "skew" => 3,
   "threshold" => 0.8,
   "fraction" => 0.5,
   "reverse" => 0
};

GetOptions($options,
           'length|l=i',
           'skew=i',
           'threshold=f',
           'fraction=i',
           'reverse|v!',
           'debug!' => \$DEBUG,
) or pod2usage(1);

my $inQual = 0; # false
my $seqID = "";
my $qualID = "";
my $seq = "";
my $qual = "";
while(<>){
  chomp;
  chomp;
  if(!$inQual){
    if(/^(>|@)((.+?)( .*?\s*)?)$/){
      my $newSeqID = $2;
      my $newShortID = $3;
      if($seqID){
        processSeq($seqID, $seq, $qual,
                   $options->{"length"}, $options->{"skew"},
                   $options->{"threshold"}, $options->{"fraction"},
                   $options->{"reverse"});
      }
      $seq = "";
      $qual = "";
      $seqID = $newSeqID;
    } elsif(/^\+(.*)$/) {
      $inQual = 1; # true
      $qualID = $1;
    } else {
      $seq .= $_;
    }
  } else {
    $qual .= $_;
    if(length($qual) >= length($seq)){
      $inQual = 0; # false
    }
  }
}

if($seqID){
  processSeq($seqID, $seq, $qual,
             $options->{"length"}, $options->{"skew"},
             $options->{"threshold"}, $options->{"fraction"},
             $options->{"reverse"});
}
