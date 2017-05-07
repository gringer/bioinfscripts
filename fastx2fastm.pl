#!/usr/bin/perl

use warnings;
use strict;

my $inQual = 0; # false
my $seqID = "";
my $qualID = "";
my $seq = "";
my $qual = "";
while(<>){
  if(!$inQual){
    if(/^(>|@)((.+?)( .*?\s*)?)$/){
      my $newSeqID = $2;
      if($seq){
	$seq =~ s/\n/\n /g;
	$seq =~ s/ $//;
	printf("#seq:%s\n", $seqID);
	printf(" %s", $seq);
	if($qual){
	  $qual =~ s/\n/\n /g;
	  $qual =~ s/ $//;
	  printf("#qual:%s\n", $seqID);
	  printf(" %s", $qual);
	}
      }
      $seq = "";
      $qual = "";
      $seqID = $newSeqID;
      chomp $seqID;
      chomp $seqID;
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
  printf(" %s", $seq);
  if($qual){
    printf("#qual:%s\n", $seqID);
    printf(" %s", $qual);
  }
}
