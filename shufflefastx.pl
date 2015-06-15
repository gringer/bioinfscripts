#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);
use List::Util qw(shuffle);

my $idFileName = "";
my $quiet = 0;
my $wrap = 1;

GetOptions("idfile=s" => \$idFileName, "quiet!" => \$quiet,
          "wrap!" => \$wrap) or
  die("Error in command line arguments");

my %idsToGet = ();

# unknown commands are treated as identifiers
my @files = ();
while(@ARGV){
  my $arg = shift(@ARGV);
  if(-f $arg){
    push(@files, $arg);
  } else {
    $idsToGet{$arg} = 1;
  }
}
@ARGV = @files;

if($idFileName){
  # read sequence IDs from input file
  open(my $idFile, "<", $idFileName);
  while(<$idFile>){
    chomp;
    s/^[>@]//;
    $idsToGet{$_} = 1;
  }
  close($idFile);
}

my $numToFind = scalar(keys(%idsToGet));

if(!$quiet && $numToFind){
  printf(STDERR "Read %d identifiers\n", $numToFind);
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
      if($seqID){
        my $newSeq = "";
        foreach my $rp (shuffle (0 .. length($seq))){
          $newSeq .= substr($seq,$rp,1);
        }
        if($qual){
          printf("@"."shuffled_%s\n%s\n+\n%s\n", $seqID, $newSeq, $qual);
        } else {
          if($wrap){
            $newSeq =~ s/(.{70})/$1\n/g;
            $newSeq =~ s/\s+$//;
          }
          printf(">shuffled_%s\n%s\n", $seqID, $newSeq);
        }
      }
      $seq = "";
      $qual = "";
      if(!$numToFind ||
         exists($idsToGet{$newSeqID}) || exists($idsToGet{$newShortID})){
        $seqID = $newSeqID;
      } else {
        $seqID = "";
      }
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
  my $newSeq = "";
  foreach my $rp (shuffle (0 .. length($seq))){
    $newSeq .= substr($seq,$rp,1);
  }
  if($wrap){
    $newSeq =~ s/(.{70})/$1\n/g;
    $newSeq =~ s/\s+$//;
  }
  if($qual){
    printf("@%s\n%s\n+\n%s\n", $seqID, $newSeq, $qual);
  } else {
    printf(">%s\n%s\n", $seqID, $newSeq);
  }
}
