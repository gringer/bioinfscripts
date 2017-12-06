#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Digest::MD5 qw(md5_hex);

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

my $showMD5 = 0; # false
my $showFName = 0; # false

GetOptions("md5!" => \$showMD5, "fname!" => \$showFName, ) or
  die("Error in command line arguments");

# use stdin if no files supplied
if(!@ARGV){
  @ARGV = '-' unless (-t STDIN);
}

my @lengths = ();
my $inQual = 0; # false
my $seqID = "";
my $qualID = "";
my $seq = "";
my $qual = "";
my $seqFName = "";
foreach my $file (@ARGV) {
  # This little gunzip dance makes sure the script can handle both
  # gzip-compressed and uncompressed input, regardless of whether
  # or not it is piped
  my $z = new IO::Uncompress::Gunzip($file, "transparent", 1)
    or die "gunzip failed: $GunzipError\n";
  while(<$z>){
    chomp;chomp;chomp;
    if (!$inQual) {
      if (/^(>|@)((.+?)( .*?\s*)?)$/) {
        my $newSeqID = $2;
        my $newShortID = $3;
        if ($seqID) {
          print(length($seq));
          if($showMD5){
            print(" ".md5_hex($seq));
          }
          print(" ${seqID}");
          if($showFName){
            print(" ".$seqFName);
          }
          print("\n");
          push(@lengths, length($seq));
        }
        $seq = "";
        $qual = "";
        $seqID = $newSeqID;
        $seqFName = $file;
      } elsif (/^\+(.*)$/) {
        $inQual = 1;            # true
        $qualID = $1;
      } else {
        $seq .= $_;
      }
    } else {
      $qual .= $_;
      if (length($qual) >= length($seq)) {
        $inQual = 0;            # false
      }
    }
  }
}

if ($seqID) {
  print(length($seq));
  if($showMD5){
    print(" ".md5_hex($seq));
  }
  print(" ${seqID}");
  if($showFName){
    print(" ".$seqFName);
  }
  print("\n");
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

my $L10LengthNum = 0;
while($cumLengths[$L10LengthNum] < ($sum * 0.1)){
  $L10LengthNum++;
}

printf(STDERR "Total sequences: %d\n", scalar(@lengths));
printf(STDERR "Total length: %sb\n", SIConvert($sum));
printf(STDERR "Longest sequence: %sb\n", SIConvert($lengths[0]));
printf(STDERR "Shortest sequence: %sb\n", SIConvert($lengths[$#lengths]));
printf(STDERR "Mean Length: %sb\n", SIConvert(sprintf("%d", ($sum) / scalar(@lengths))));
printf(STDERR "Median Length: %sb\n", SIConvert($lengths[$#lengths / 2]));
printf(STDERR "N10: %d sequences; L10: %sb\n",
     ($L10LengthNum+1), SIConvert($lengths[$L10LengthNum]));
printf(STDERR "N50: %d sequences; L50: %sb\n",
     ($L50LengthNum+1), SIConvert($lengths[$L50LengthNum]));
printf(STDERR "N90: %d sequences; L90: %sb\n",
     ($L90LengthNum+1), SIConvert($lengths[$L90LengthNum]));
