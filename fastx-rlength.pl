#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);
#use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error);
use IO::File;

my $trim = 0;
my $maxUnit = 1000;
my $sampleSize = 1000; ## number of positions to check
my $maxChunks = 10; ## maximum number of chunks to split a contig into

GetOptions("trim=s" => \$trim) or
  die("Error in command line arguments");

my @rlengths = ();
my $inQual = 0; # false
my $seqID = "";
my $qualID = "";
my $seq = "";
my $qual = "";
my $buffer = "";
while(<>){
  chomp;
  chomp;
  if(!$inQual){
    if(/^(>|@)((.+?)( .*?\s*)?)$/){
      my $newSeqID = $2;
      my $newShortID = $3;
      my $len = length($seq);
      if($seqID && (length($seq) > $trim)){
        my $scoreTotal = 0;
        my $scoreMax = 0;
        my $scoreMin = $sampleSize;
        my $scoreCount = 0;
        my $maxOffs = 0;
        my $minOffs = 0;
        if ($len > ($maxUnit * 2)) {
          my $localMU = (($len / 5) > $maxUnit) ? $maxUnit : int($len / 5);
          ## limit search region to 10*$maxUnit bp chunks, to allow for non-global-repeats
          my $chunkCount = ($len > ($localMU * 10)) ? int($len / ($localMU * 10) + 1) : 1;
          $chunkCount = $maxChunks if $chunkCount > $maxChunks;
          my $chunkSize = int($len / $chunkCount + 1);
          my @seqArr = split("",$seq);
          for (my $chunkOffs = 0; ($chunkOffs + $chunkSize) < $len; $chunkOffs += $chunkSize) {
            ## Scanning for repeats within the range [$chunkOffs .. $chunkOffs + $chunkSize]
            for (my $slipOffs = 4; $slipOffs < $localMU; $slipOffs++) {
              my $score = 0;
              for (my $i = 0; $i < ($sampleSize); $i++) {
                my $seqPos = int(rand($chunkSize - $slipOffs));
                $score++ if ($seqArr[$chunkOffs + $seqPos] eq
                             $seqArr[$chunkOffs + $seqPos + $slipOffs]);
              }
              if ($score > $scoreMax) {
                $maxOffs = $slipOffs;
                $scoreMax = $score;
              }
              if ($score < $scoreMin) {
                $minOffs = $slipOffs;
                $scoreMin = $score;
              }
              $scoreTotal += $score;
              $scoreCount++;
            }
          }
        }
        my $scoreMean = ($scoreCount) ? ($scoreTotal / $scoreCount) : 0;
        my $minFrac = $scoreMean ? ($scoreMean - $scoreMin) / ($scoreMax - $scoreMean) : 0;
        my $maxFrac = $scoreMean ? ($scoreMax - $scoreMean) / ($scoreMean - $scoreMin) : 0;
        printf("%8d %0.3f %0.3f %0.3f %3d %0.3f %3d %0.3f %s\n",
               $len, $scoreMin/$sampleSize, $scoreMean/$sampleSize, $scoreMax/$sampleSize,
               $minOffs, $minFrac, $maxOffs, $maxFrac, $seqID);
        push(@rlengths, $maxOffs);
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
if ($seqID && (length($seq) > $trim)) {
  my $scoreTotal = 0;
  my $scoreMax = 0;
  my $scoreMin = $sampleSize;
  my $scoreCount = 0;
  my $maxOffs = 0;
  my $minOffs = 0;
  if ($len > ($maxUnit * 2)) {
    my $localMU = (($len / 5) > $maxUnit) ? $maxUnit : int($len / 5);
    ## limit search region to 10*$maxUnit bp chunks, to allow for non-global-repeats
    my $chunkCount = ($len > ($localMU * 10)) ? int($len / ($localMU * 10) + 1) : 1;
    $chunkCount = $maxChunks if $chunkCount > $maxChunks;
    my $chunkSize = int($len / $chunkCount + 1);
    my @seqArr = split("",$seq);
    for (my $chunkOffs = 0; ($chunkOffs + $chunkSize) < $len; $chunkOffs += $chunkSize) {
      ## Scanning for repeats within the range [$chunkOffs .. $chunkOffs + $chunkSize]
      for (my $slipOffs = 4; $slipOffs < $localMU; $slipOffs++) {
        my $score = 0;
        for (my $i = 0; $i < ($sampleSize); $i++) {
          my $seqPos = int(rand($chunkSize - $slipOffs));
          $score++ if ($seqArr[$chunkOffs + $seqPos] eq
                       $seqArr[$chunkOffs + $seqPos + $slipOffs]);
        }
        if ($score > $scoreMax) {
          $maxOffs = $slipOffs;
          $scoreMax = $score;
        }
        if ($score < $scoreMin) {
          $minOffs = $slipOffs;
          $scoreMin = $score;
        }
        $scoreTotal += $score;
        $scoreCount++;
      }
    }
  }
  my $scoreMean = ($scoreCount) ? ($scoreTotal / $scoreCount) : 0;
  my $minFrac = $scoreMean ? ($scoreMean - $scoreMin) / ($scoreMax - $scoreMean) : 0;
  my $maxFrac = $scoreMean ? ($scoreMax - $scoreMean) / ($scoreMean - $scoreMin) : 0;
  printf("%8d %0.3f %0.3f %0.3f %3d %0.3f %3d %0.3f %s\n",
         $len, $scoreMin/$sampleSize, $scoreMean/$sampleSize, $scoreMax/$sampleSize,
         $minOffs, $minFrac, $maxOffs, $maxFrac, $seqID);
  push(@rlengths, $maxOffs);
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
