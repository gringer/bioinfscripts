#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);

my $quiet = 0;
my $filterSeqs = "";
my $reverse = 0;

GetOptions("filter=s" => \$filterSeqs, "reverse|v!" => \$reverse, "quiet!" => \$quiet) or
  die("Error in command line arguments");

# unknown commands are treated as identifiers
my @files = ();
while(@ARGV){
  my $arg = shift(@ARGV);
  if(-f $arg){
    push(@files, $arg);
  } else {
    $filterSeqs .= "|$arg";
  }
}
@ARGV = @files;

$filterSeqs =~ s/^\|//;
$filterSeqs = "($filterSeqs)";

if($reverse){
    printf(STDERR "Filter sequence (excluded): $filterSeqs\n");
} else {
    printf(STDERR "Filter sequence: $filterSeqs\n");
}

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
      if($seqID && ($reverse xor ($seq =~ /$filterSeqs/))){
        if($qual){
          printf("@%s\n%s\n+\n%s\n", $seqID, $seq, $qual);
        } else {
          printf(">%s\n%s\n", $seqID, $seq);
        }
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

if($seqID && ($reverse xor ($seq =~ /$filterSeqs/))){
  if($qual){
    printf("@%s\n%s\n+\n%s\n", $seqID, $seq, $qual);
  } else {
    printf(">%s\n%s\n", $seqID, $seq);
  }
}
