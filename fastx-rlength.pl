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

my $trim = 0;
my $kmerLength = 17; ## number of bases in hash keys

GetOptions("trim=s" => \$trim) or
    die("Error in command line arguments");

my @rlengths = ();
my $inQual = 0; # false
my $seqID = "";
my $qualID = "";
my $seq = "";
my $sseq = "";
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
      my $len = length($seq);
      if($seqID && (length($seq) > $trim) && (length($seq) > $kmerLength)){
        my $countTotal = 0;
        my $countMax = 0;
        my $maxKmer = "";
	my @rptCounts = ();
	my %rptHash = ();
        my %gapCounts = ();
	for(my $p = 0; ($p + $kmerLength) <= $len; $p++){
          $sseq = substr($seq, $p, $kmerLength);
          if($sseq !~ /^([ACGT][ACGT])\1{6,}$/){
            push(@{$rptHash{$sseq}}, $p);
          }
        }
	my $numKmers = scalar(keys(%rptHash));
	my $kmerRatio = $numKmers/($len - $kmerLength + 1);
	my @repeatedKmers = grep {scalar(@{$rptHash{$_}}) > 1} keys(%rptHash);
        my @revKmers = grep {defined $rptHash{rev($_)}} keys(%rptHash);
        my @rcKmers = grep {defined $rptHash{rc($_)}} keys(%rptHash);
	my @gaps = ();
	foreach my $kmer (@repeatedKmers){
	  my @posList = @{$rptHash{$kmer}};
	  my $posCount = scalar(@posList);
	  $countTotal += $posCount;
	  if($posCount > $countMax){
	    $countMax = $posCount;
	  }
	  if(@posList){
	    push(@rptCounts, scalar(@posList));
	    #print(join(" ",@posList)."\n");
	    my $lastPos = $posList[0];
	    for(my $p = 0; $p <= $#posList; $p++){
	      my $nextPos = $posList[$p];
              my $gap = $nextPos - $lastPos;
              if($gap > 0){
                push(@gaps, $gap);
                $gapCounts{$gap}++;
              }
	      $lastPos = $nextPos;
	    }
	  }
	}
	@gaps = sort {$a <=> $b} (@gaps);
	@rptCounts = sort {$a <=> $b} (@rptCounts);
	my $medianGap = (@gaps) ? $gaps[$#gaps / 2] : 0;
	my $medianCount = $medianGap ? ${gapCounts{$medianGap}} : 0;
        my $modalGap = 0;
        my $modalCount = 0;
        my $rangeCountMed = 0;
        my $rangeCountMod = 0;
        if($medianCount){
          my @modalSort = sort {$gapCounts{$b} <=> $gapCounts{$a}} (@gaps);
          $modalGap = $modalSort[0];
          $modalCount = $gapCounts{$modalGap};
          for(my $gP = int($medianGap * 0.99); ($gP <= ($medianGap / 0.99));
              $gP++){
            $rangeCountMed += $gapCounts{$gP} if($gapCounts{$gP});
          }
          for(my $gP = int($modalGap * 0.99); ($gP <= ($modalGap / 0.99));
              $gP++){
            $rangeCountMod += $gapCounts{$gP} if($gapCounts{$gP});
          }
        }
        my $numRepeats = scalar(@repeatedKmers);
        printf("%8d %0.3f %6d %6d %5d %5d %6d %5d %5d %6d %6d %6d %s\n",
               $len, $kmerRatio, scalar(@repeatedKmers),
               $countTotal,
               $medianCount,
               $rangeCountMed,
               $medianGap,
               $modalCount,
               $rangeCountMod,
               $modalGap,
               scalar(@revKmers),
               scalar(@rcKmers),
               $seqID);
        push(@rlengths, $medianGap) if $medianGap;
      }
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

my $len = length($seq);
if($seqID && (length($seq) > $trim) && (length($seq) > $kmerLength)){
  my $countTotal = 0;
  my $countMax = 0;
  my $maxKmer = "";
  my @rptCounts = ();
  my %rptHash = ();
  my %gapCounts = ();
  for(my $p = 0; ($p + $kmerLength) <= $len; $p++){
    $sseq = substr($seq, $p, $kmerLength);
    if($sseq !~ /^([ACGT][ACGT])\1{6,}$/){
      push(@{$rptHash{$sseq}}, $p);
    }
  }
  my $numKmers = scalar(keys(%rptHash));
  my $kmerRatio = $numKmers/($len - $kmerLength + 1);
  my @repeatedKmers = grep {scalar(@{$rptHash{$_}}) > 1} keys(%rptHash);
  my @revKmers = grep {defined $rptHash{rev($_)}} keys(%rptHash);
  my @rcKmers = grep {defined $rptHash{rc($_)}} keys(%rptHash);
  my @gaps = ();
  foreach my $kmer (@repeatedKmers) {
    my @posList = @{$rptHash{$kmer}};
    my $posCount = scalar(@posList);
    $countTotal += $posCount;
    if ($posCount > $countMax) {
      $countMax = $posCount;
    }
    if (@posList) {
      push(@rptCounts, scalar(@posList));
      #print(join(" ",@posList)."\n");
      my $lastPos = $posList[0];
      for (my $p = 0; $p <= $#posList; $p++) {
        my $nextPos = $posList[$p];
        my $gap = $nextPos - $lastPos;
        if ($gap > 0) {
          push(@gaps, $gap);
          $gapCounts{$gap}++;
        }
        $lastPos = $nextPos;
      }
    }
  }
  @gaps = sort {$a <=> $b} (@gaps);
  @rptCounts = sort {$a <=> $b} (@rptCounts);
  my $medianGap = (@gaps) ? $gaps[$#gaps / 2] : 0;
  my $medianCount = $medianGap ? ${gapCounts{$medianGap}} : 0;
  my $modalGap = 0;
  my $modalCount = 0;
  my $rangeCountMed = 0;
  my $rangeCountMod = 0;
  if ($medianCount) {
    my @modalSort = sort {$gapCounts{$b} <=> $gapCounts{$a}} (@gaps);
    $modalGap = $modalSort[0];
    $modalCount = $gapCounts{$modalGap};
    for (my $gP = int($medianGap * 0.99); ($gP <= ($medianGap / 0.99));
         $gP++) {
      $rangeCountMed += $gapCounts{$gP} if($gapCounts{$gP});
    }
    for (my $gP = int($modalGap * 0.99); ($gP <= ($modalGap / 0.99));
         $gP++) {
      $rangeCountMod += $gapCounts{$gP} if($gapCounts{$gP});
    }
  }
  my $numRepeats = scalar(@repeatedKmers);
  printf("%8d %0.3f %6d %6d %5d %5d %6d %5d %5d %6d %6d %6d %s\n",
         $len, $kmerRatio, scalar(@repeatedKmers),
         $countTotal,
         $medianCount,
         $rangeCountMed,
         $medianGap,
         $modalCount,
         $rangeCountMod,
         $modalGap,
         scalar(@revKmers),
         scalar(@rcKmers),
         $seqID);
  push(@rlengths, $medianGap) if $medianGap;
}

## calculate statistics
@rlengths = sort {$b <=> $a} (@rlengths);
my $sum = 0;
my @cumLengths = map {$sum += $_} (@rlengths);

printf(STDERR "Total sequences: %d\n", scalar(@rlengths));
printf(STDERR "Longest repeat: %0.3f\n", $rlengths[0]);
printf(STDERR "Shortest repeat: %0.3f\n", $rlengths[$#rlengths]);
printf(STDERR "Mean repeat length: %0.3f\n", $sum / scalar(@rlengths));
printf(STDERR "Median repeat length: %0.3f\n", $rlengths[$#rlengths / 2]);
