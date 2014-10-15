#!/usr/bin/perl

use warnings;
use strict;

use Pod::Usage; ## uses pod documentation in usage code
use Getopt::Long qw(:config auto_version auto_help);

our $VERSION = "0.1";

=head1 NAME

nanopore_barcode.pl -- generates a random barcode suitable for
the nanopore sequencer

=head1 SYNOPSIS

./nanopore_barcode.pl -length <int> -count <int>

=head2 Options

=over 2

=item B<-length>

Length of barcode to generate (default 30)

=item B<-count>

Number of barcodes to generate (default 1)

=back

=head1 DESCRIPTION

You will need to verify independently that the chosen barcodes are at
least 20% different from each other.

=cut

my $length = 30;
my $count = 1;

GetOptions('length=i' => \$length,
           'count=i'=> \$count,
) or pod2usage(1);

my $lastBase = "";

for(my $seqNum = 1; $seqNum <= $count; $seqNum++){
    my $sequence = "";
    my $lastBase = "C";
    for(my $i = 0; $i < $length; $i++){
        my @nextBaseOptions = grep {$_ ne $lastBase} ("A","T","C");
        my $nextBase = (rand() < 0.6) ?
            $nextBaseOptions[0] : $nextBaseOptions[1];
        $sequence .= $nextBase;
        $lastBase = $nextBase;
    }
    printf(">Barcode_%03d\n%s\n", $seqNum, $sequence);
}

