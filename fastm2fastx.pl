#!/usr/bin/perl

use warnings;
use strict;

my $fastqMode = 0; # false
my $inQual = 0;
my $inSeq = 0;

my $seqID = "";
my $seq = "";
my $qual = "";

while(<>){
  if(/^#seq(:(.*$))?/){
    my $newSeqID = $2;
    $inQual = 0; # false
    $inSeq = 1; # true
    if($seqID){
      if($fastqMode){
	printf("@%s\n%s+\n%s", $seqID, $seq, $qual);
      } else {
	printf(">%s\n%s", $seqID, $seq);
      }
    }
    $seqID = $newSeqID;
    $seq = "";
  } elsif(/^#qual(:(.*$))?/){
    my $newSeqID = $2;
    if($newSeqID && ($newSeqID ne $seqID)){
      warn(sprintf("Quality ID [%s] and sequence ID [%s] do not match.".
		   " Setting ID to sequence ID [%s]", 
		   $newSeqID, $seqID, $seqID));
    }
    $fastqMode = 1;
    $inSeq = 0; # false
    $inQual = 1; # true
  } elsif(/^#/){ # some other unknown tag
    $inSeq = 0; # false
    $inQual = 0; # false
  }
  if(/^ (.*$)/){
    if($inSeq){
      $seq .= $1 . "\n";
    }
    if($inQual){
      $qual .= $1 . "\n";
    }
  }
}

if($seqID){
  if($fastqMode){
    printf("@%s\n%s+\n%s", $seqID, $seq, $qual);
  } else {
    printf(">%s\n%s", $seqID, $seq);
  }
}
