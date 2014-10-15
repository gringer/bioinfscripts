#!/usr/bin/perl

use warnings;
use strict;

use Pod::Usage; ## uses pod documentation in usage code
use Getopt::Long qw(:config auto_version auto_help);

our $VERSION = "0.2";

=head1 NAME

nanopore_barcode.pl -- generates a random barcode suitable for
the nanopore sequencer

=head1 SYNOPSIS

./nanopore_barcode.pl [options]

=head2 Options

=over 2

=item B<-length>

Length of barcode to generate (default 15)

=item B<-count>

Number of barcodes to generate (default 1)

=item B<-threshold>

Difference threshold for barcode inclusion (default 0.2)

=back

=head1 DESCRIPTION

Note: when more than one barcode is generated, each new barcode is
checked to make sure it is sufficiently different from previously
added sequences.

=cut

my $length = 15;
my $count = 1;
my $threshold = 0.2;

GetOptions('length=i' => \$length,
           'count=i'=> \$count,
           'threshold=f'=> \$threshold,
) or pod2usage(1);

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

while($seqNum < $count){
    $sequence = "";
    my $lastBase = "C";
    for(my $i = 0; $i < $length; $i++){
        my @nextBaseOptions = grep {$_ ne $lastBase} ("A","T","C");
        my $nextBase = (rand() < 0.6) ?
            $nextBaseOptions[0] : $nextBaseOptions[1];
        $sequence .= $nextBase;
        $lastBase = $nextBase;
    }
    my $addedDiffCount =
      grep {seqsDifferent($_,$sequence,$threshold)} @addedSeqs;
    if(!@addedSeqs || ($addedDiffCount == scalar(@addedSeqs))){
      $seqNum++;
      printf(">Barcode_%03d\n%s\n", $seqNum, $sequence);
      push(@addedSeqs, $sequence);
    }
}
