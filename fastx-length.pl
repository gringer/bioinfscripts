#!/usr/bin/perl
use warnings;
use strict;

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

my @lengths = ();
my $inQual = 0; # false
my $seqID = "";
my $qualID = "";
my $seq = "";
my $qual = "";
while(<>){
  chomp;
  chomp;
  if(!$inQual){
    if(/^(>|@)((.+?)( .*?\s*)?)$/){
      my $newSeqID = $2;
      my $newShortID = $3;
      if($seqID){
        printf("%d %s\n", length($seq), $seqID);
	push(@lengths, length($seq));
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
    printf("%d %s\n", length($seq), $seqID);
    push(@lengths, length($seq));
}

## calculate statistics
@lengths = sort {$b <=> $a} (@lengths);
my $sum = 0;
my @cumLengths = map {$sum += $_} (@lengths);

my $L50LengthNum = 0;
while($cumLengths[$L50LengthNum] < ($sum * 0.5)){
  $L50LengthNum++;
}

my $L90LengthNum = 0;
while($cumLengths[$L90LengthNum] < ($sum * 0.9)){
  $L90LengthNum++;
}

printf(STDERR "Total sequences: %d\n", scalar(@lengths));
printf(STDERR "Total length: %sb\n", SIConvert($sum));
printf(STDERR "Longest sequence: %sb\n", SIConvert($lengths[0]));
printf(STDERR "Shortest sequence: %sb\n", SIConvert($lengths[$#lengths]));
printf(STDERR "Mean Length: %sb\n", SIConvert(sprintf("%d", ($sum) / scalar(@lengths))));
printf(STDERR "Median Length: %sb\n", SIConvert($lengths[$#lengths / 2]));
printf(STDERR "N50: %d sequences; L50: %sb\n",
     ($L50LengthNum+1), SIConvert($lengths[$L50LengthNum]));
printf(STDERR "N90: %d sequences; L90: %sb\n",
     ($L90LengthNum+1), SIConvert($lengths[$L90LengthNum]));
