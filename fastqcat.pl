#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_version auto_help pass_through); # for option parsing

my $inQual = 0; # false
my $seqID = "";
my $qualID = "";
my $seq = "";
my $qual = "";
my $minLen = 0;

GetOptions('minlength=i' => \$minLen) or
  die("Error in command line arguments");

while(<>){
  chomp;
  chomp;
  if(!$inQual){
    if(/^@(.+)$/){
      my $newSeqID = $1;
      if($seqID && (length($seq) >= $minLen)){
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

if($seqID && (length($seq) >= $minLen)){
  printf("@%s\n%s\n+\n%s\n", $seqID, $seq, $qual);
}
