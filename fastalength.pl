#!/usr/bin/perl
use warnings;
use strict;

my $sortable = 0; # false

if($ARGV[0] eq "-s"){
  shift(@ARGV);
  $sortable = 1; # true
}

my $seq = "";
my $seqID = "";
my $keep = 0;
while(<>){
  chomp;
  if(/^>((.+?)( .*?\s*)?)$/){
    my $newID = $1;
    my $newShortID = $2;
    if($seq){
      if($sortable){
        printf("%d %s\n", length($seq), $seqID);
      } else {
        printf(">%s [%d bp]\n", $seqID, length($seq));
      }
    }
    $seq = "";
    $seqID = $newID;
  } else {
    $seq .= $_;
  }
}
if($seq){
  if($sortable){
    printf("%d %s\n", length($seq), $seqID);
  } else {
    printf(">%s [%d bp]\n", $seqID, length($seq));
  }
}
