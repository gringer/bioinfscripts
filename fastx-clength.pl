#!/usr/bin/perl
use warnings;
use strict;

use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error);
use IO::File;

sub SIConvert{
  my ($val) = @_;
  my @largePrefixes = ("k", "M", "G");
  my $pID = 0;
  my $prefix = "";
  my $changed=0;
  while(($val > 1000) && ($pID < $#largePrefixes)){
    $prefix = $largePrefixes[$pID];
    $val /= 1000;
    $pID++;
  }
  return(sprintf("%.12g %s", $val, $prefix));
}

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
      if($seqID){
	bzip2 \$seq => \$buffer;
	my $cProp = (length($buffer) * 1000) / length($seq);
        printf("%0.3f %s\n", (1000 / $cProp), $seqID);
	push(@lengths, 1000 / $cProp);
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

if($seqID){
  bzip2 \$seq => \$buffer;
  my $cProp = (length($buffer) * 1000) / length($seq);
  printf("%0.3f %s\n", (1000 / $cProp), $seqID);
  push(@lengths, 1000 / $cProp);
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
