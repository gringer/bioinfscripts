#!/usr/bin/perl
use warnings;
use strict;

my $inQual = 0; # false
my $seqID = "";
my $qualID = "";
my $seq = "";
my $qual = "";
while(<>){
  chomp;
  chomp;
  if(!$inQual){
    if(/^@(.+)$/){
      my $newSeqID = $1;
      if($seqID){
        printf("@%s\n%s\n+\n%s\n", $seqID, $seq, $qual);
      }
      $seq = "";
      $qual = "";
      $seqID = $newSeqID;
    } elsif(/^\+(.*)$/) {
      $inQual = 1; # true
      $qualID = $1;
    } else {
      $seq .= $_;
    }
  } else {
    $qual .= $_;
    if(length($qual) >= length($seq)){
      $inQual = 0; # false
    }
  }
}

if($seqID){
  printf("@%s\n%s\n+\n%s\n", $seqID, $seq, $qual);
}
