#!/usr/bin/perl

# vcf_windowSorter.pl -- sorts a VCF file by position within a
# specified window size. The effect of this is to convert a file
# containing overlapping gene regions (with a given window of overlap)
# into a sorted file that can be indexed.

use warnings;
use strict;

my $windowSize = 50000;
 my $currentChr = "";
my $minPos = -1;
my $maxPos = -1;

my %posCache = ();
my @storedPoss = ();

while(<>){
    my $line = $_;
    if(!/^#/){
	my ($chr, $pos) = (/^(.*?)\t([0-9]+)/);
	if($chr ne $currentChr){
	    # spit out cache
	    for my $oldPos (@storedPoss){
		print($posCache{$oldPos});
	    }
	    # reset min/max values and clear cache
	    $currentChr = $chr;
	    $minPos = -1;
	    $maxPos = -1;
	    %posCache = ();
	    @storedPoss = ();
	}
	#print(join(":",($chr, $pos))."\n");
	if($pos > ($minPos + $windowSize)){
	    #printf(STDERR "Pos: $pos, MinPos: $minPos, MaxPos: $maxPos, ArraySize: %d\n", scalar(@storedPoss));
	    while(@storedPoss && ($storedPoss[0] <= ($pos - $windowSize))){
		# spit out first bit of cache
		my $oldPos = shift(@storedPoss);
		# extract and remove position information from cache
		print(delete($posCache{$oldPos}));
	    }
	    if(!(@storedPoss)){ # no remaining values, so no maximum value
		$maxPos = -1;
	    }
	}
	if(!exists($posCache{$pos})){ # avoid inserting positions multiple times
	    # pos might need to slot into array, so check with max value
	    if($pos < $maxPos){
		# this is expensive, but hopefully doesn't happen too often
		@storedPoss = sort({$a <=> $b} ($pos, @storedPoss));
		# no change to maxPos
	    } else {
		push(@storedPoss, $pos);
		$maxPos = $pos;
	    }
	}
	# store line in cache
	$posCache{$pos} .= $line;
	# regenerate minimum value; storedPoss should have at least 1 element
	$minPos = $storedPoss[0];
    } else {
	print($_);
    }
}

# print out remaining cached locations
for my $oldPos (@storedPoss){
    print($posCache{$oldPos});
}
