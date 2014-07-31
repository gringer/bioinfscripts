#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);

my $idFileName = "";
my $quiet = 0;

GetOptions("idfile=s" => \$idFileName, "quiet!" => \$quiet) or
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
    s/^>//;
    $idsToGet{$_} = 1;
  }
  close($idFile);
}

if(!$quiet){
  printf(STDERR "Read %d identifiers\n", scalar(keys(%idsToGet)));
}

my $seq = "";
my $seqID = "";
my $keep = 0;
while(<>){
  chomp;
  if(/^>((.+?)( .*?\s*)?)$/){
    my $newID = $1;
    my $newShortID = $2;
    if($seq){
      printf(">%s\n%s\n", $seqID, $seq);
    }
    $seq = "";
    if(exists($idsToGet{$newID}) || exists($idsToGet{$newShortID})){
      $seqID = $1;
      $seq = "";
      $keep = 1;
    } else {
      #printf(STDERR ">%s [%s]\n", $newID, $newShortID);
      $seqID = "";
      $keep = 0;
    }
  } elsif($keep) {
    $seq .= $_;
  }
}
if($seq){
  printf(">%s\n%s\n", $seqID, $seq);
}
