#!/usr/bin/perl
use warnings;
use strict;

my $sortable = 0; # false

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


if($ARGV[0] eq "-s"){
  shift(@ARGV);
  $sortable = 1; # true
}

my $seq = "";
my $seqID = "";
my $keep = 0;
my @lengths = ();
while(<>){
  chomp;
  if(/^>((.+?)( .*?\s*)?)$/){
    my $newID = $1;
    my $newShortID = $2;
    if($seq){
      if($sortable){
        printf("%d %s\n", length($seq), $seqID);
      } else {
        printf(">%s [%d bp]\n", $seqID, length($seq));
      }
      push(@lengths, length($seq));
    }
    $seq = "";
    $seqID = $newID;
  } else {
    $seq .= $_;
  }
}
if($seq){
  if($sortable){
    printf("%d %s\n", length($seq), $seqID);
  } else {
    printf(">%s [%d bp]\n", $seqID, length($seq));
  }
}

## calculate statistics
@lengths = sort {$b <=> $a} (@lengths);
my $sum = 0;
my @cumLengths = map {$sum += $_} (@lengths);
my $L50Length = $sum * 0.5;
my @L50cumLengths = grep {$_ < $L50Length} @cumLengths;
my $L50LengthNum = $#L50cumLengths;
if($L50cumLengths[$L50LengthNum] < $L50Length){
  $L50LengthNum++;
}
my $L90Length = $sum * 0.9;
my @L90cumLengths = grep {$_ < $L90Length} @cumLengths;
my $L90LengthNum = $#L90cumLengths;
if($L90cumLengths[$L90LengthNum] < $L90Length){
  $L90LengthNum++;
}

printf(STDERR "Total sequences: %d\n", scalar(@lengths));
printf(STDERR "Total length: %sbp\n", SIConvert($sum));
printf(STDERR "Longest sequence: %sbp\n", SIConvert($lengths[0]));
printf(STDERR "Shortest sequence: %sbp\n", SIConvert($lengths[$#lengths]));
printf(STDERR "N50: %d sequences; L50: %sbp\n",
     $L50LengthNum, SIConvert($lengths[$L50LengthNum]));
printf(STDERR "N90: %d sequences; L90: %sbp\n",
     $L90LengthNum, SIConvert($lengths[$L90LengthNum]));
