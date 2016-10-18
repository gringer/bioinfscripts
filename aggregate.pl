#!/usr/bin/perl

use warnings;
use strict;

my $keyName = "";
my @valueAggregate = ();

while(<>){
  chomp;
  my @F = split(",", $_);
  if($F[0] ne $keyName){
    if($keyName){
      printf("${keyName},%d\n", scalar(@valueAggregate));
    }
    $keyName = $F[0];
    @valueAggregate = ();
  }
  push(@valueAggregate, $F[1]);
}
