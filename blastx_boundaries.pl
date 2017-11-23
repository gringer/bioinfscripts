#!/usr/bin/perl

use warnings;
use strict;

my %starts = ();
my %ends = ();

while(<>){
  chomp;
  my @F = split(/\s+/);
  my $qs=$F[6];
  my $qe=$F[7];
  if($qe < $qs){
    ($qs, $qe) = ($qe, $qs);
  }
  my $sig = $F[0]."-".$F[1];
  if(!$starts{$sig} || ($starts{$sig} > $qs)){
    $starts{$sig} = $qs
  }
  if(!$ends{$sig} || ($ends{$sig} < $qe)){
    $ends{$sig} = $qe
  }
}

foreach my $key (keys(%starts)){
  if($key =~ /(.*)-(.*)/){
    my $target = $1;
    my $query = $2;
    printf("%s %s:%d-%d\n", $query, $target, $starts{$key}, $ends{$key});
  }
}
