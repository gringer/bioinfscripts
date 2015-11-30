#!/usr/bin/perl

use warnings;
use strict;

use Pod::Usage; ## uses pod documentation in usage code
use Getopt::Long qw(:config auto_version auto_help);

our $VERSION = "0.3";

=head1 NAME

nanopore_barcode.pl -- generates a random barcode suitable for
the nanopore sequencer

=head1 SYNOPSIS

./nanopore_barcode.pl [options]

=head2 Options

=over 2

=item B<-length>

Length of barcode to generate, excluding prefix/suffix (default 38)

=item B<-count>

Number of barcodes to generate (default 1)

=item B<-threshold>

Difference threshold for barcode inclusion (default 0.2)

=item B<-prefix>

Prefix sequence for barcode (default I<GGTGCTG>)

=item B<-suffix>

Suffix sequence for barcode (default I<TTAACCT>)

=item B<-order>

Order of base selection (default I<ATCG>)

=item B<-baseprob>

Base probability, excluding previous base (commma-separated, default I<0.45,0.30,0.25>)

=back

=head1 DESCRIPTION

Note: when more than one barcode is generated, each new barcode is
checked to make sure it is sufficiently different from previously
added sequences.

=cut

my $length = 38;
my $count = 1;
my $threshold = 0.2;
my $prefix = "GGTGCTG";
my $suffix = "TTAACCT";
my $baseOrder = "ATCG";
my $baseProb = "0.25,0.50,0.75";

GetOptions('length=i' => \$length,
           'count=i'=> \$count,
           'threshold=f'=> \$threshold,
           'prefix=s' => \$prefix,
           'suffix=s' => \$suffix,
           'order=s' => \$baseOrder,
           'baseprob=s' => \$baseProb,
          ) or pod2usage(1);

my @baseProbs = split(/,/,$baseProb);
my @baseOrders = split(//, $baseOrder);

my $lastBase = "";

sub seqsDifferent {
  my ($seq1, $seq2, $diffThreshold) = @_;
  my $differences = 0;
  for(my $i = 0; $i < length($seq1); $i++){
    if(substr($seq1,$i,1) ne substr($seq2,$i,1)){
      $differences++;
    }
  }
  my $diffProp = $differences / length($seq1);
  return($diffProp >= $diffThreshold);
}

my @addedSeqs = ();
my $sequence = "";
my $seqNum = 0;

my $maxHPlength = 3;

while($seqNum < $count){
  $sequence = "";
  my $lastHP = (length($prefix) == 0) ? "C" : substr($prefix,-1,1);
  for(my $i = 0; $i < $length; $i++){
    my @nextBaseOptions = (length($lastHP) < $maxHPlength) ?
      @baseOrders :
        grep {$_ ne substr($lastHP,0,1)} @baseOrders;
    my $prob = rand();
    my $nextBase = $nextBaseOptions[$#nextBaseOptions];
    my $probMultiplier = (scalar(@nextBaseOptions) == 4) ? 1 : (4/3);
    if($prob < ($baseProbs[0] * $probMultiplier)){
      $nextBase = $nextBaseOptions[0];
    } elsif($prob < ($baseProbs[1] * $probMultiplier)){
      $nextBase = $nextBaseOptions[1];
    } elsif($prob < ($baseProbs[2] * $probMultiplier)){
      $nextBase = $nextBaseOptions[2];
    }
    $sequence .= $nextBase;
    if($nextBase eq substr($lastHP,0,1)){
      $lastHP .= $nextBase;
    } else {
      $lastHP = $nextBase;
    }
  }
  my $addedDiffCount =
    grep {seqsDifferent($_,$sequence,$threshold)} @addedSeqs;
  if(!@addedSeqs || ($addedDiffCount == scalar(@addedSeqs))){
    $seqNum++;
    printf(">Barcode_%03d\n%s\n", $seqNum, $prefix.$sequence.$suffix);
    push(@addedSeqs, $sequence);
  }
}
