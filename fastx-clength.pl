#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);
use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error);
use IO::File;

my $trim = 0;

GetOptions("trim=s" => \$trim) or
  die("Error in command line arguments");

my @clengths = ();
my @lengths = ();
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
      if($seqID && ($len > 1) && ($len > $trim)){
	bzip2 \$seq => \$buffer;
	my $cProp = length($seq) / length($buffer);
        # repetition statistic, appears to correlate with repetitiveness
        # a non-repetitive sequence tends to have this number < 5
        my $cStat = exp($cProp) / log($len);
        # normalised version of the repetitiveness statistic
        # [approximately normal when repetitive sequences are excluded]
        my $ncStat = ($cProp) / log($len);
        my ($ltProp, $midProp, $rtProp) = (0, 0, 0);
        my $bTrim = ($trim < 1) ? ($len * $trim) : $trim;
        if($trim && (($bTrim * 3) < $len)){
          my $ltSeq = substr($seq, 0, $bTrim);
          my $midSeq = substr($seq, $bTrim, -$bTrim);
          my $rtSeq = substr($seq, -$bTrim);
          bzip2 \$ltSeq => \$buffer;
          $ltProp = $bTrim / length($buffer);
          bzip2 \$midSeq => \$buffer;
          $midProp = ($len - $bTrim * 2) / length($buffer);
          bzip2 \$rtSeq => \$buffer;
          $rtProp = $bTrim / length($buffer);
        }
        printf("%0.3f %d %0.3f %0.3f %0.3f %0.3f %0.3f %s\n", $cProp, $len, $cStat, $ncStat, $ltProp, $midProp, $rtProp, $seqID);
	push(@lengths, $cProp);
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

if($seqID && (length($seq) > $trim)){
  my $len = length($seq);
  bzip2 \$seq => \$buffer;
  my $cProp = length($seq) / length($buffer);
  my ($ltProp, $midProp, $rtProp) = (0, 0, 0);
  my $bTrim = ($trim < 1) ? ($len * $trim) : $trim;
  if ($trim && (($bTrim * 3) < $len)) {
    my $ltSeq = substr($seq, 0, $bTrim);
    my $midSeq = substr($seq, $bTrim, -$bTrim);
    my $rtSeq = substr($seq, -$bTrim);
    bzip2 \$ltSeq => \$buffer;
    $ltProp = $bTrim / length($buffer);
    bzip2 \$midSeq => \$buffer;
    $midProp = ($len - $bTrim * 2) / length($buffer);
    bzip2 \$rtSeq => \$buffer;
    $rtProp = $bTrim / length($buffer);
  }
  printf("%0.3f %d %0.3f %0.3f %0.3f %s\n", $cProp, $len,
         $ltProp, $midProp, $rtProp, $seqID);
  push(@lengths, $cProp);
}

## calculate statistics
@lengths = sort {$b <=> $a} (@lengths);
my $sum = 0;
my @cumLengths = map {$sum += $_} (@lengths);

printf(STDERR "Total sequences: %d\n", scalar(@lengths));
printf(STDERR "Highest compression: %0.3f\n", $lengths[0]);
printf(STDERR "Lowest compression: %0.3f\n", $lengths[$#lengths]);
printf(STDERR "Mean compression: %0.3f\n", $sum / scalar(@lengths));
printf(STDERR "Median compression: %0.3f\n", $lengths[$#lengths / 2]);
