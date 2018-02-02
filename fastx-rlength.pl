#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);
#use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error);
#use IO::File;

sub rc {
  my ($seq) = @_;
  $seq =~ tr/ACGTUYRSWMKDVHBXN-/TGCAARYSWKMHBDVXN-/;
  # work on masked sequences as well
  $seq =~ tr/acgtuyrswmkdvhbxn/tgcaaryswkmhbdvxn/;
  return(scalar(reverse($seq)));
}

sub rev {
  my ($seq) = @_;
  return(scalar(reverse($seq)));
}

sub printStats {
  my ($seq, $seqID, $trim, $kmerLength) = @_;
  my $len = length($seq);
  my $sseq = "";
  if($seqID && (length($seq) > $trim) && (length($seq) > $kmerLength)){
    my $countTotal = 0;
    my $countMax = 0;
    my $maxKmer = "";
    my %rptPos = ();
    my %allGapCounts = ();
    my %minGaps = ();
    my $revCount = 1;
    my $rcCount = 1;
    for(my $p = 0; ($p + $kmerLength) <= $len; $p++){
      $sseq = substr($seq, $p, $kmerLength);
      if(exists($rptPos{$sseq})){
	my $gap = $p - $rptPos{$sseq};
	$allGapCounts{$gap}++;
	if(!exists($minGaps{$sseq}) || ($minGaps{$sseq} > $gap)){
	  $minGaps{$sseq} = $gap;
	}
      }
      if(exists($rptPos{rev($sseq)})){
	$revCount++;
      }
      if(exists($rptPos{rc($sseq)})){
	$rcCount++;
      }
      $rptPos{$sseq} = $p;
    }
    if($revCount == 1){
      $revCount = 0;
    }
    if($rcCount == 1){
      $rcCount = 0;
    }
    my $numKmers = scalar(keys(%rptPos));
    my $kmerRatio = $numKmers/($len - $kmerLength + 1);
    my $numRepeats = scalar(keys(%minGaps));
    my @gaps = sort {$a <=> $b} (values(%minGaps));
    my $medianGap = (@gaps) ? $gaps[$#gaps / 2] : 0;
    my $medianCount = 0;
    my $modalGap = 0;
    my $modalCount = 0;
    my $rangeCountMed = 0;
    my $rangeCountMod = 0;
    if($medianGap){
      my %gapCounts = ();
      foreach my $gap (@gaps){
	$gapCounts{$gap}++;
      }
      $medianCount = ${allGapCounts{$medianGap}};
      my @modalSort = sort {$allGapCounts{$b} <=> $allGapCounts{$a}} (@gaps);
      $modalGap = $modalSort[0];
      $modalCount = $allGapCounts{$modalGap};
      for(my $gP = int($medianGap * 0.99); ($gP <= ($medianGap / 0.99));
	  $gP++){
	$rangeCountMed += $allGapCounts{$gP} if($allGapCounts{$gP});
      }
      for(my $gP = int($modalGap * 0.99); ($gP <= ($modalGap / 0.99));
	  $gP++){
	$rangeCountMod += $allGapCounts{$gP} if($allGapCounts{$gP});
      }
    }
    printf("%8d %0.3f %6d %6d %5d %5d %6d %5d %5d %6d %6d %6d %s\n",
	   $len, $kmerRatio, 
	   $numRepeats,
	   $countTotal,
	   $medianCount,
	   $rangeCountMed,
	   $medianGap,
	   $modalCount,
	   $rangeCountMod,
	   $modalGap,
	   $revCount,
	   $rcCount,
	   $seqID);
  }
}

my $trim = 0;
my $kmerLength = 17; ## number of bases in hash keys

GetOptions("trim=s" => \$trim) or
    die("Error in command line arguments");

my $inQual = 0; # false
my $seqID = "";
my $qualID = "";
my $seq = "";
my $qual = "";
my $buffer = "";
printf("%8s %5s %6s %6s %5s %5s %6s %5s %5s %6s %6s %6s %s\n",
       "length", "kRat", "cntRep", "cntTot",
       "medCt", "RCMed", "medGap",  "modCt", "RCMod",
       "modGap", "rvCnt", "rcCnt", "SeqID");
while(<>){
  chomp;
  chomp;
  if(!$inQual){
    if(/^(>|@)((.+?)( .*?\s*)?)$/){
      my $newSeqID = $2;
      my $newShortID = $3;
      printStats($seq, $seqID, $trim, $kmerLength);
      $seq = "";
      $qual = "";
      $buffer = "";
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

printStats($seq, $seqID, $trim, $kmerLength);
