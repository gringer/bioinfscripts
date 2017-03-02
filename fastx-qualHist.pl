#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $idFileName = "";
my $quiet = 0;
my $base = 33;

GetOptions("idfile=s" => \$idFileName, "quiet!" => \$quiet,
           "base" => \$base ) or
  die("Error in command line arguments");

# unknown commands are treated as identifiers
my @files = ();
while(@ARGV){
  my $arg = shift(@ARGV);
  if(-f $arg){
    push(@files, $arg);
  }
}
@ARGV = @files;

my %qualCounts = ();

my $baseCount = 0;

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
    grep {$qualCounts{$_}++} split(//,$_);
    $baseCount += length($_);
    $qual .= $_;
    if(length($qual) >= length($seq)){
      $inQual = 0; # false
    }
  }
}

my $cumCount = 0;
foreach my $qualChar (sort(keys(%qualCounts))){
  $cumCount += $qualCounts{$qualChar};
  printf("%s [%2d]: %6d (%6.2f%% / %6.2f%% )\n",
         $qualChar, ord($qualChar) - $base,
         $qualCounts{$qualChar}, $qualCounts{$qualChar} * 100 / $baseCount,
        $cumCount * 100 / $baseCount);
}

printf("Total sequence length: %d\n", $baseCount);
