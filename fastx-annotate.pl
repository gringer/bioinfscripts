#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);

my $idFileName = "";
my $quiet = 0;

GetOptions("idfile=s" => \$idFileName, "quiet!" => \$quiet) or
  die("Error in command line arguments");

my %idsToAnnotate = ();

# unknown commands are treated as identifiers
my @files = ();
while(@ARGV){
  my $arg = shift(@ARGV);
  if(-f $arg){
    push(@files, $arg);
  } else {
    $idsToAnnotate{$arg} = shift(@ARGV);
  }
}
@ARGV = @files;

if($idFileName){
  # read sequence IDs from input file
  printf(STDERR "Reading from file '$idFileName'\n");
  open(my $idFile, "<", $idFileName);
  while(<$idFile>){
    chomp;
    s/^[>@]//;
    s/\"//g;
    my @F = split(/,/,$_,2);
    $idsToAnnotate{$F[0]} = $F[1];
  }
  close($idFile);
}

if(!$quiet){
  printf(STDERR "Read %d identifiers\n", scalar(keys(%idsToAnnotate)));
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
        if($qual){
          printf("@%s\n%s\n+\n%s\n", $seqID, $seq, $qual);
        } else {
          printf(">%s\n%s\n", $seqID, $seq);
        }
      }
      $seq = "";
      $qual = "";
      if(exists($idsToAnnotate{$newSeqID})){
        $seqID = $newSeqID . " " . $idsToAnnotate{$newSeqID};
      } elsif(exists($idsToAnnotate{$newShortID})){
        $seqID = $newShortID . " " . $idsToAnnotate{$newShortID};
      } else {
        $seqID = $newSeqID;
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
  if($qual){
    printf("@%s\n%s\n+\n%s\n", $seqID, $seq, $qual);
  } else {
    printf(">%s\n%s\n", $seqID, $seq);
  }
}
