#!/usr/bin/perl

use warnings;
use strict;

my ($tfamFileName) = shift(@ARGV);

my @maleGTs = ();
my @femaleGTs = ();

#print(STDERR "Processing TFAM file...");
open(my $tfamFile, "<", $tfamFileName) or die("Cannot open TFAM file '$tfamFileName'");
my $tfamPos = 0;
my $lastGTPos = -1;
while(<$tfamFile>){
    chomp;
    my ($iid, $fid, $mat, $pat, $sex, $pheno) = split(/\s+/);
    if($sex == 1){
	push(@maleGTs, $tfamPos);
    } elsif($sex == 2){
	push(@femaleGTs, $tfamPos);
    }
    $tfamPos++;
    $lastGTPos++;
}
close($tfamFile);
#print(STDERR " done\n");

my @maleAlleles = map {($_*2, $_*2+1)} @maleGTs;
my @femaleAlleles = map {($_*2, $_*2+1)} @femaleGTs;

my $processed = 0;

print(STDERR ".");
#printf("Marker,Sex,Type,Bases,Count\n");
my %counts = ();
while(<>){
    chomp;
    my $line = $_;
    my ($chr, $marker, $mapPos, $genPos, @alleles) = split(/\s+/, $_);
    my @genotypes = map {@alleles[$_*2].@alleles[$_*2+1]} (0..$lastGTPos);
    %counts = ();
    grep {
	my $gt = $_;
	my $a1 = substr($gt,0,1);
	my $a2 = substr($gt,1,1);
	$counts{"a"}{$gt}++;
	$counts{"a"}{$a1}++;
	$counts{"a"}{$a2}++;
	$counts{"m"}{$gt}++;
	$counts{"m"}{$a1}++;
	$counts{"m"}{$a2}++;
    } @genotypes[@maleGTs];
    grep {
	my $gt = $_;
	my $a1 = substr($gt,0,1);
	my $a2 = substr($gt,1,1);
	$counts{"a"}{$gt}++;
	$counts{"a"}{$a1}++;
	$counts{"a"}{$a2}++;
	$counts{"f"}{$gt}++;
	$counts{"f"}{$a1}++;
	$counts{"f"}{$a2}++;
    } @genotypes[@femaleGTs];
    foreach my $sex ("a","m","f"){
	foreach my $bases (keys(%{$counts{$sex}})){
	    printf("%s,%s,%s,%s,%d\n", $marker, $sex,
		   (length($bases) == 1) ? "a" : "g",
		   $bases,
		   $counts{$sex}{$bases});
	}
    }
    $processed++;
    if($processed > 1000){
	print(STDERR ".");
	$processed = 0;
    }
}
#print(STDERR " done\n");
