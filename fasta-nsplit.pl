#!/usr/bin/perl
use warnings;
use strict;

my $seq = "";
my $shortSeqID = "";
my $seqID = "";
my $keep = 0;
my $cumLength = 0;
while(<>){
  chomp;
  if(/^>((.+?)( .*?\s*)?)$/){
    my $newID = $1;
    my $newShortID = $2;
    if($seq){
      my $inc = 0;
      while($seq =~ s/(NNNN+)(.*)//){
	my $nStretch = $1;
        my $newSeq = $2;
        printf(">%s.%s\n%s\n", $seqID, $inc++, $seq) if ($seq);
	$cumLength += length($seq);
	printf(STDERR "%s\t%d\t%d\n", $shortSeqID, $cumLength, 
	       $cumLength + length($nStretch));
	$cumLength += length($nStretch);
        $seq = $newSeq;
      }
      printf(">%s\n%s\n", $seqID, $seq) if ($seq);
    }
    $seq = "";
    $shortSeqID = $newShortID;
    $seqID = $newID;
    $cumLength = 0;
  } else {
    $seq .= $_;
  }
}
if($seq){
  my $inc = 0;
  while($seq =~ s/(NNNN+)(.*)//){
    my $nStretch = $1;
    my $newSeq = $2;
    printf(">%s.%s\n%s\n", $seqID, $inc++, $seq) if ($seq);
    $cumLength += length($seq);
    printf(STDERR "%s\t%d\t%d\n", $shortSeqID, $cumLength, 
	   $cumLength + length($nStretch));
    $cumLength += length($nStretch);
    $seq = $newSeq;
  }
  printf(">%s\n%s\n", $seqID, $seq) if ($seq);
}
