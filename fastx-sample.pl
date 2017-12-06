#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use File::Temp qw(:seekable);
use Encode qw(encode_utf8);

my $idFileName = "";
my $maxCount = 0;
my $trim = 0;

GetOptions("count=i" => \$maxCount, ) or
  die("Error in command line arguments");

if(!$maxCount){
  die("Error: Reservoir sampling needs a maximum number of reads (-c <count>)\n");
} else {
  printf(STDERR "Reservoir sampling reads to output at most %d reads:",
        $maxCount);
}

sub processReads{
  my ($seq, $qual, $maxCount, $readsRead, $readsProcessed, $readStore) = @_;
  $readsProcessed++;
  my $swapPos = ($readsProcessed <= $maxCount) ? ($readsProcessed-1) : int(rand($readsProcessed));
  if($swapPos < $maxCount){
    my $outLines = "";
    if($qual){
      $outLines = sprintf("@%012d\n%s\n+\n%s\n", $readsRead, $seq, $qual);
    } else {
      $outLines = sprintf(">%012d\n%s\n", $readsRead, $seq);
    }
    ${$readStore}{$swapPos} = $outLines;
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
my %readStore = ();

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
              printf(STDERR " (%d reads processed)", $readsRead);
            }
            printf(STDERR "\n  ");
          }
          print(STDERR ".");
          $dotsPrinted++;
        }
        $readsProcessed =
          processReads($seq, $qual, $maxCount, ++$readsRead, $readsProcessed, \%readStore);
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
          processReads($seq, $qual, $maxCount, ++$readsRead, $readsProcessed, \%readStore);
}

printf(STDERR "\ndone (%d reads processed from %d total reads)\n",
       $readsProcessed, $readsRead);

foreach my $id (keys(%readStore)){
  print $readStore{$id};
}
