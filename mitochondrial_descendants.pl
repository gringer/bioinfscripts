#!/usr/bin/perl

use warnings;
use strict;

my %foundingLine = ();
my %found = ();
my %maternalIDs = ();
my %sex = ();

while (<>) {
    my @F = split(/\s+/);
    $maternalIDs{$F[1]} = $F[3];
    $sex{$F[1]} = $F[4];
    if(($F[2] == 0) && ($F[3] == 0)){
	$foundingLine{$F[1]}{$F[1]} = 1;
    }
}

print(STDERR "Loaded all IDs\n");

my $changed = 1;

print(STDERR "Finding mitochondrial descendants...");
while($changed){
    $changed = 0;
    foreach my $line (keys(%foundingLine)){
	my %lineInds = %{$foundingLine{$line}};
	foreach my $ind (keys(%lineInds)){
	    if(!$found{$ind}){
		#printf(STDERR "Finding mitochondrial descendants of <%d>\n", $ind);
		foreach my $ind2 (keys(%maternalIDs)){
		    if($maternalIDs{$ind2} == $ind){
			$foundingLine{$line}{$ind2} = 1;
		    }
		}
		$found{$ind} = 1;
		$changed = 1;
	    }
	}
    }
}

print(STDERR " stored all descendants\n");

foreach my $line (sort {$a <=> $b} (keys(%foundingLine))){
    if(scalar(keys(%{$foundingLine{$line}})) > 1){
	printf("Mitochondrial descendants of <%s>:\n", $line);
	my %lineInds = %{$foundingLine{$line}};
	foreach my $ind (sort {$a <=> $b} (keys(%lineInds))){
	    printf(" $ind %d\n", $sex{$ind});
	}
    }
}
