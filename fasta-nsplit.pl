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
      my $inc = 0;
      while($seq =~ s/NNNN+(.*)//){
        my $newSeq = $1;
        printf(">%s.%s\n%s\n", $seqID, $inc++, $seq);
        $seq = $newSeq;
      }
      printf(">%s\n%s\n", $seqID, $seq);
    }
    $seq = "";
    $seqID = $newID;
  } else {
    $seq .= $_;
  }
}
if($seq){
  my $inc = 0;
  while($seq =~ s/N+(.*)//){
    my $newSeq = $1;
    printf(">%s.%s\n%s\n", $seqID, $inc++, $seq);
  }
  printf(">%s\n%s\n", $seqID, $seq);
}
