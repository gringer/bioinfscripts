#!/usr/bin/perl
use warnings;
use strict;

## fastx-isofilter.pl -- Extract the longest isoform (or ORF) from transcripts

use Getopt::Long qw(:config auto_help pass_through);

my $quiet = 0;
my $trimString = "_i[0-9]+\$";
my $orfMode = 0;

GetOptions("trim=s" => \$trimString, "orf!" => \$orfMode,
           "quiet!" => \$quiet) or
  die("Error in command line arguments");

if($orfMode){
  $trimString = "_[0-9]+\$";
}

# unknown commands are treated as identifiers
my @files = ();
while(@ARGV){
  my $arg = shift(@ARGV);
  if(-f $arg){
    push(@files, $arg);
  } else {
    $trimString .= "|$arg";
  }
}
@ARGV = @files;

# use stdin if no files supplied
if(!@ARGV){
  @ARGV = '-' unless (-t STDIN);
}

if($trimString){
  $trimString =~ s/^\|//;
  $trimString = "($trimString)";
}

my $inQual = 0; # false
my $seqID = "";
my $fullID = "";
my $qualID = "";
my $seq = "";
my $qual = "";
my %fastXStrs = ();
my %fastXLengths = ();
while(<>){
  chomp;
  chomp;
  if(!$inQual){
    if(/^(>|@)((.+?)( .*?\s*)?)$/){
      my $newFullID = $2;
      my $newSeqID = $3;
      $newSeqID =~ s/$trimString//;
      if($seqID && (!$fastXLengths{$seqID} || ($fastXLengths{$seqID} < length($seq)))){
        $fastXLengths{$seqID} = length($seq);
        if(!$qual){
          $seq =~ s/(.{100})/$1\n/g;
          $seq =~ s/\n$//;
        }
        $fastXStrs{$seqID} = ($qual) ?
          sprintf("@%s\n%s\n+\n%s\n", $fullID, $seq, $qual) :
          sprintf(">%s\n%s\n", $fullID, $seq);
      }
      $seq = "";
      $qual = "";
      $seqID = $newSeqID;
      $fullID = $newFullID;
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

if($seqID && (!$fastXLengths{$seqID} || ($fastXLengths{$seqID} < length($seq)))){
  $fastXLengths{$seqID} = length($seq);
  if(!$qual){
    $seq =~ s/(.{100})/$1\n/g;
    $seq =~ s/\n$//;
  }
  $fastXStrs{$seqID} = ($qual) ?
    sprintf("@%s\n%s\n+\n%s\n", $seqID, $seq, $qual) :
    sprintf(">%s\n%s\n", $seqID, $seq);
}

foreach my $pat (sort {$fastXLengths{$b} <=> $fastXLengths{$a}} (keys(%fastXStrs))){
  print($fastXStrs{$pat});
}
