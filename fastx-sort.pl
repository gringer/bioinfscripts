#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);

my $quiet = 0;
my $searchPattern = "(^.*\$)";

GetOptions("pattern=s" => \$searchPattern, "quiet!" => \$quiet) or
  die("Error in command line arguments");

# Complain about non-file command line argument
my @files = ();
while(@ARGV){
  my $arg = shift(@ARGV);
  if(-f $arg){
    push(@files, $arg);
  } else {
    die("Unknown argument: $arg");
  }
}
@ARGV = @files;

my %fastXStrs = ();

my $inQual = 0; # false
my $seqID = "";
my $qualID = "";
my $seq = "";
my $qual = "";
while(<>){
  chomp;
  chomp;
  if(!$inQual){
    if(/^(>|@)(.*)$/){
      my $newSeqID = $2;
      if($seqID){
        if($seqID =~ /$searchPattern/){
          $fastXStrs{$seqID} = ($qual) ?
            sprintf("@%s\n%s\n+\n%s\n", $seqID, $seq, $qual) :
            sprintf(">%s\n%s\n", $seqID, $seq);
        } else {
          printf(STDERR "Warning: No match for pattern '$searchPattern' for sequence '$seqID'\n");
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

if($seqID){
  if($seqID =~ /$searchPattern/){
    $fastXStrs{$seqID} = ($qual) ?
      sprintf("@%s\n%s\n+\n%s\n", $seqID, $seq, $qual) :
      sprintf(">%s\n%s\n", $seqID, $seq);
  } else {
    printf(STDERR "Warning: No match for pattern '$searchPattern' for sequence '$seqID'\n");
  }
}

foreach my $pat (sort(keys(%fastXStrs))){
  print($fastXStrs{$pat});
}
