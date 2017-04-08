#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use File::Temp qw(:seekable);
use Encode qw(encode_utf8);

my $idFileName = "";
my $readLength = 100; # 1 Tbp
my $maxCount = 0;
my $trim = 0;

GetOptions("readLength|l=i" => \$readLength,
           "count=i" => \$maxCount, ) or
  die("Error in command line arguments");

if(!$maxCount){
  die("Error: Reservoir sampling needs a maximum number of reads (-c <count>)\n");
} else {
  printf(STDERR "Reservoir sampling reads to output at most %d reads:",
        $maxCount);
}

sub processReads{
  my ($seq, $qual, $readLength, $maxCount, $readsRead, $readsProcessed, $tempFile) = @_;
  if(length($seq) < $readLength){
    return $readsProcessed;
  }
  if(length($seq) > $readLength){
    $seq = substr($seq, 0, $readLength);
    $qual = substr($qual, 0, $readLength) if $qual;
  }
  $readsProcessed++;
  my $swapPos = $readsProcessed-1;
  if($readsProcessed > $maxCount){
    $swapPos = int(rand($readsProcessed));
  }
  if($swapPos < $maxCount){
    my $outLines = "";
    if($qual){
      $outLines = sprintf("@%012d\n%s\n+\n%s\n", $readsRead, $seq, $qual);
    } else {
      $outLines = sprintf(">%012d\n%s\n", $readsRead, $seq);
    }
    if($readsProcessed > $maxCount){
      seek($tempFile, $swapPos * length(encode_utf8($outLines)), 0);
    }
    print($tempFile $outLines);
  }
  return $readsProcessed;
}

my $inQual = 0; # false
my $seqID = "";
my $qualID = "";
my $seq = "";
my $qual = "";
my $readsRead = 0;
my $readsProcessed = 0;
my $dotsPrinted = 0;
my $tempFile = File::Temp->new();

while(<>){
  chomp;
  chomp;
  if(!$inQual){
    if(/^(>|@)((.+?)( .*?\s*)?)$/){
      my $newSeqID = $2;
      my $newShortID = $3;
      if($seqID){
        if($maxCount && ($readsRead % 10000 == 0)){
          if($dotsPrinted % 50 == 0){
            if($readsRead > 1000){
              printf(STDERR " (%d reads read)", $readsRead);
            }
            printf(STDERR "\n  ");
          }
          print(STDERR ".");
          $dotsPrinted++;
        }
        $readsProcessed =
          processReads($seq, $qual, $readLength, $maxCount, ++$readsRead, $readsProcessed, $tempFile);
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

if($seqID){
  $readsProcessed =
          processReads($seq, $qual, $readLength, $maxCount, ++$readsRead, $readsProcessed, $tempFile);
}

printf(STDERR "\ndone (%d reads processed of length >= %d from %d total reads)\n",
       $readsProcessed, $readLength, $readsRead);

seek($tempFile, 0, SEEK_SET);
while(<$tempFile>){
  print;
}
