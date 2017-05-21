#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);
use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error);
use IO::File;

my $trim = 0;

GetOptions("trim=i" => \$trim) or
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
      if($seqID && (length($seq) > $trim)){
	my $len = length($seq);
	bzip2 \$seq => \$buffer;
	my $cProp = length($seq) / length($buffer);
        my ($ltProp, $midProp, $rtProp) = (0, 0, 0);
        if(($trim * 3) < $len){
          my $ltSeq = substr($seq, 0, $trim);
          my $midSeq = substr($seq, $trim, -$trim);
          my $rtSeq = substr($seq, -$trim);
          bzip2 \$ltSeq => \$buffer;
          $ltProp = $trim / length($buffer);
          bzip2 \$midSeq => \$buffer;
          $midProp = ($len - $trim * 2) / length($buffer);
          bzip2 \$rtSeq => \$buffer;
          $rtProp = $trim / length($buffer);
        }
        printf("%0.3f %d %0.3f %0.3f %0.3f %s\n", $cProp, $len, $ltProp, $midProp, $rtProp, $seqID);
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
  if (($trim * 3) < $len) {
    my $ltSeq = substr($seq, 0, $trim);
    my $midSeq = substr($seq, $trim, -$trim);
    my $rtSeq = substr($seq, -$trim);
    bzip2 \$ltSeq => \$buffer;
    $ltProp = $trim / length($buffer);
    bzip2 \$midSeq => \$buffer;
    $midProp = ($len - $trim * 2) / length($buffer);
    bzip2 \$rtSeq => \$buffer;
    $rtProp = $trim / length($buffer);
  }
  printf("%0.3f %d %0.3f %0.3f %0.3f %s\n", $cProp, $len, $ltProp, $midProp, $rtProp, $seqID);
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
