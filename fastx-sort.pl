#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);

my $quiet = 0;
my $searchPattern = ""; ## "(^.*\$)";
my $numeric = 0;
my $length = 0;
my $r = 0;

GetOptions("reverse|r!" => \$r, "pattern=s" => \$searchPattern,
           "quiet!" => \$quiet,
           "numeric|n!" => \$numeric, "length!" => \$length) or
  die("Error in command line arguments");

if($r){
  print(STDERR "Reversing sort direction\n");
}

# Complain about non-file command line argument
my @files = ();
while(@ARGV){
  my $arg = shift(@ARGV);
  if(-e $arg){
    push(@files, $arg);
  } else {
    die("Unknown argument: $arg");
  }
}
@ARGV = @files;

my %fastXStrs = ();
my %fastXVals = ();

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
        if(!$searchPattern || ($seqID =~ /($searchPattern)/)){
          my $matchPattern = ($2) ? $2 : $1;
          $fastXVals{$seqID} = (!$searchPattern) ? length($seq) : $matchPattern;
          if(!$qual){
            $seq =~ s/(.{100})/$1\n/g;
            $seq =~ s/\n$//;
          }
          $fastXStrs{$seqID} .= ($qual) ?
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

if ($seqID) {
  if (!$searchPattern || ($seqID =~ /($searchPattern)/)) {
    my $matchPattern = ($2) ? $2 : $1;
    $fastXVals{$seqID} = (!$searchPattern) ? length($seq) : $matchPattern;
    if (!$qual) {
      $seq =~ s/(.{100})/$1\n/g;
      $seq =~ s/\n$//;
    }
    $fastXStrs{$seqID} .= ($qual) ?
      sprintf("@%s\n%s\n+\n%s\n", $seqID, $seq, $qual) :
      sprintf(">%s\n%s\n", $seqID, $seq);
  } else {
    printf(STDERR "Warning: No match for pattern '$searchPattern' for sequence '$seqID'\n");
  }
}

if($searchPattern){
  if($numeric){
    foreach my $pat (sort {$fastXVals{$r?$b:$a} <=> $fastXVals{$r?$a:$b}} (keys(%fastXStrs))){
      print($fastXStrs{$pat});
    }
  } else {
    foreach my $pat (sort {$fastXVals{$r?$b:$a} cmp $fastXVals{$r?$a:$b}} (keys(%fastXStrs))){
      print($fastXStrs{$pat});
    }
  }
} elsif($length){
  foreach my $pat (sort {$fastXVals{$r?$a:$b} <=> $fastXVals{$r?$b:$a}} (keys(%fastXStrs))){
    print($fastXStrs{$pat});
  }
} elsif($numeric){
  foreach my $pat (sort {$r?$b:$a <=> $r?$a:$b} (keys(%fastXStrs))){
    print($fastXStrs{$pat});
  }
} else {
  foreach my $pat (sort {$r?$b:$a cmp $r?$a:$b} (keys(%fastXStrs))){
    print($fastXStrs{$pat});
  }
}
