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
      $seqID = $1;
      $seq = "";
    } elsif(/^\+(.*)$/) {
      $inQual = 1; # true
      $qualID = $1;
      $qual = "";
      my $printedSeq = $seq;
      $printedSeq =~ s/(.{60})/$1\n/g;
      printf(">%s\n%s\n", $seqID, $printedSeq);
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
