#!/usr/bin/perl
use warnings;
use strict;

my $seq = "";
my $seqID = "";
my $keep = 0;
while(<>){
  chomp;
  if(/^>((.+?)( .*?\s*)?)$/){
    my $newID = $1;
    my $newShortID = $2;
    if($seq){
      printf(">%s [%d bp]\n", $seqID, length($seq));
    }
    $seq = "";
    $seqID = $newID;
  } else {
    $seq .= $_;
  }
}
if($seq){
  printf(">%s [%d bp]\n", $seqID, length($seq));
}
