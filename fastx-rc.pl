#!/usr/bin/perl
use warnings;
use strict;

## fastx-rc -- reverse complement a sequence

use Getopt::Long qw(:config auto_help pass_through);

my $quiet = 0;

sub rc {
  my ($seq) = @_;
  $seq =~ tr/ACGTUYRSWMKDVHBXN-/TGCAARYSWKMHBDVXN-/;
  # work on masked sequences as well
  $seq =~ tr/acgtuyrswmkdvhbxn/tgcaaryswkmhbdvxn/;
  return(scalar(reverse($seq)));
}

GetOptions("quiet!" => \$quiet) or
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
          printf("@%s [RC]\n%s\n+\n%s\n", $seqID, rc($seq), scalar(reverse($qual)));
        } else {
          printf(">%s [RC]\n%s\n", $seqID, rc($seq));
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
  if($qual){
    printf("@%s [rc]\n%s\n+\n%s\n", $seqID, rc($seq), scalar(reverse($qual)));
  } else {
    printf(">%s [rc]\n%s\n", $seqID, rc($seq));
  }
}
