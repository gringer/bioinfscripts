#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);

my $minLength = 0;
my $singleLine = 0;

GetOptions("minLength=i" => \$minLength, "singleLine!" => \$singleLine) or
  die("Error in command line arguments");

my $inQual = 0; # false
my $seqID = "";
my $qualID = "";
my $seq = "";
my $qual = "";
while(<>){
  chomp;
  chomp;
  if(/^\s+$/){
    next;
  }
  if(!$inQual){
    if(/^@(.+)$/){
      $seqID = $1;
      $seq = "";
    } elsif(/^\+(.*)$/) {
      $inQual = 1; # true
      $qualID = $1;
      $qual = "";
      if(length($seq) > $minLength){
        my $printedSeq = $seq;
        if(!$singleLine){
          $printedSeq =~ s/(.{60})/$1\n/g;
          $printedSeq =~ s/\n$//g;
        }
        printf(">%s\n%s\n", $seqID, $printedSeq);
      }
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
