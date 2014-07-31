#!/usr/bin/perl
use warnings;
use strict;

my $snpFileName = shift(@ARGV);
open(my $snpFile, "<", $snpFileName);

my %inMap = ();
while(<$snpFile>){
    chomp;
    $inMap{$_} = 1;
}
close($snpFile);

$inMap{"dbSNP RS ID"} = 1;

while(<>){
    if(/^#/){
        print;
    } else {
        my @fields = split(/\t/);
        if($inMap{$fields[14]}){
            print;
        }
    }
}
