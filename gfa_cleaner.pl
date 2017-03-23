#!/usr/bin/perl

## gfa_cleaner.pl -- removes duplicate links and links that include
## vertexes with high numbers of outgoing or incoming edges

use warnings;
use strict;

my %seen = ();
my %counts = ();
my %displayBuffer = ();

my $maxCount = 4;

while(<>){
  if(/^[^L]/){
    print;
    next;
  }
  my $line = $_;
  chomp;
  my @F = split(/\t/);
  my $matchFwd = join(";",@F[(1,2,3,4)]);
  my $strA = join(";",@F[(1,2)]);
  grep {tr/\-\+/\+\-/} @F[(4,2)];
  my $matchRev = join(";",@F[(3,4,1,2)]);
  my $strB = join(";",@F[(3,4)]); ## after sign flip
  if(!$seen{$matchFwd}){
    $counts{$strA}++;
    $counts{$strB}++;
    $seen{$matchFwd} = 1;
    $seen{$matchRev} = 1;
    $displayBuffer{$line}{$strA} = 1;
    $displayBuffer{$line}{$strB} = 1;
  }
}

foreach my $line (sort(keys(%displayBuffer))){
  my @strs = keys(%{$displayBuffer{$line}});
  if(grep {$counts{$_} > $maxCount} @strs){
    next;
  }
  print($line);
}
