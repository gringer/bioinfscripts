#!/usr/bin/perl

use warnings;
use strict;

use List::MoreUtils qw(pairwise);

my @rows = ();
my $delim = "\t";

while(<>){
    chomp;
    if(/(\s)/){
	$delim = $1;
    }
    my @F = split(/\s+/);
    if(!@rows){
	@rows = @F;
    } else {
	@rows = pairwise {$a . $delim . $b} @rows, @F;
    }
}

foreach my $row (@rows){
    print($row."\n");
}
