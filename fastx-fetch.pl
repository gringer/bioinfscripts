#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $idFileName = "";
my $quiet = 0;
my $minLength = 0;
my $maxLength = 10 ** 12; # 1 Tbp
my $count = -1;
my $invert = 0; # invert logic
my $trim = 0;

GetOptions("idfile=s" => \$idFileName, "quiet!" => \$quiet,
           "reverse|v!" => \$invert, "trim=i" => \$trim,
           "minLength=i" => \$minLength, "maxLength=i" => \$maxLength,
           "count=i" => \$count, ) or
  die("Error in command line arguments");

my %idsToGet = ();
if($trim){
  $minLength = $minLength + $trim * 2;
  $maxLength = $maxLength + $trim * 2;
}

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
  if(!$quiet){
    printf(STDERR "Attempting to read from input file ($idFileName)\n");
  }
  my $idFile = new IO::Uncompress::Gunzip "$idFileName" or
    die "Unable to open $idFileName\n";
  while(<$idFile>){
    chomp;
    s/^[>@]//;
    $idsToGet{$_} = 1;
  }
  close($idFile);
}

if(!$quiet){
  printf(STDERR "Read %d identifiers\n", scalar(keys(%idsToGet)));
}

if(!$quiet && $invert){
  printf(STDERR "Excluding IDs, rather than selecting\n");
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
      if($seqID && (length($seq) >= $minLength) && (length($seq) <= $maxLength)){
        if($trim > 0){
          $seq = substr($seq, $trim, length($seq)-($trim * 2));
          if($qual){
            $qual = substr($qual, $trim, length($seq)-($trim * 2));
          }
        }
        if($qual){
          printf("@%s\n%s\n+\n%s\n", $seqID, $seq, $qual);
        } else {
          printf(">%s\n%s\n", $seqID, $seq);
        }
      }
      if($count-- == 0){
        $seqID = "";
        last;
      }
      $seq = "";
      $qual = "";
      if((!(keys(%idsToGet)) || exists($idsToGet{$newSeqID}) || exists($idsToGet{$newShortID})) xor $invert){
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

if($seqID && (length($seq) >= $minLength) && (length($seq) <= $maxLength)){
  if($qual){
    printf("@%s\n%s\n+\n%s\n", $seqID, $seq, $qual);
  } else {
    printf(">%s\n%s\n", $seqID, $seq);
  }
}
