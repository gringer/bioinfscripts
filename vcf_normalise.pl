#!/usr/bin/perl

## vcf_normalise.pl -- removes VCF information from all but the
## highest coverage position, considering SNPs/INDELs separately

## RSRS    12      .       TATCA   TCATCA,TCAATCA  0       .       INDEL;IDV=13;IMF=0.0634146;DP=92629;...
## RSRS    12      .       T       C,A,G   0       .       DP=92629;...
## RSRS    12      .       T       C,G,A   0       .       DP=7490;...
## RSRS    12      .       T       TT      0       .       INDEL;IDV=1;IMF=0.2;DP=7490;...
## =>
## RSRS    12      .       TATCA   TCATCA,TCAATCA  0       .       INDEL;IDV=13;IMF=0.0634146;DP=92629;...
## RSRS    12      .       T       C,A,G   0       .       DP=92629;...

use warnings;
use strict;

my $indelCache = "";
my $snpCache = "";
my $bestSnpCov = 0;
my $bestIndelCov = 0;
my $oldChr = "";
my $oldLoc = -1;
my $chr = "";
my $loc = 0;
my $indel = "";
my $cov = 0;

while(<>){
    if(/^#/){
	print;
	next;
    }
    my $line = $_;
    if($line =~ /^(.*?)\t([0-9]+)\t.*?(INDEL.*?;)?DP=([0-9]+)/){
	($chr, $loc, $indel, $cov) = ($1, $2, $3, $4);
    }
    if(($chr ne $oldChr) || ($loc != $oldLoc)){
	print($indelCache);
	print($snpCache);
	if($indel){
	    $indelCache = $line;
	    $bestIndelCov = $cov;
	    $snpCache = "";
	    $bestSnpCov = 0;
	} else {
	    $indelCache = "";
	    $bestIndelCov = 0;
	    $snpCache = $line;
	    $bestSnpCov = $cov;
	}
	$oldChr = $chr;
	$oldLoc = $loc;
    } elsif($indel){
	if($cov > $bestIndelCov){
	    $indelCache = $line;
	    $bestIndelCov = $cov;
	}
    } else{
	if($cov > $bestSnpCov){
	    $snpCache = $line;
	    $bestSnpCov = $cov;
	}
    }
}

print($indelCache);
print($snpCache);
