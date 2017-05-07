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
    if(/^(>|@)((.+?)( .*?\s*)?)$/){
      my $newSeqID = $2;
      if($seq){
	printf("#seq:%s\n", $seqID);
	printf(" %s\n", $seq);
	if($qual){
	  printf("#qual:%s\n", $seqID);
	  printf(" %s\n", $qual);
	}
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
  printf("#seq:%s\n", $seqID);
  printf(" %s\n", $seq);
  if($qual){
    printf("#qual:%s\n", $seqID);
    printf(" %s\n", $qual);
  }
}
