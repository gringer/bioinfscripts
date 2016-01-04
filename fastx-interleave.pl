#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $singleLine = 0;

my %seqFiles = ();
my @seqFileOrder = ();

GetOptions("singleLine!" => \$singleLine) or
  die("Error in command line arguments");

while(@ARGV){
  my $argument = shift(@ARGV);
  if(-f $argument){ # file existence check
    $seqFiles{$argument} = 1;
    push(@seqFileOrder, $argument);
  } else {
    die("Error in command line arguments: $argument");
  }
}

foreach my $seqFileName (keys(%seqFiles)){
  my $tmpFile = new IO::Uncompress::Gunzip "$seqFileName" or
    die "Unable to open $seqFileName\n";
  $seqFiles{$seqFileName} = $tmpFile;
}

my $numFiles = scalar(@seqFileOrder);

if($numFiles < 2){
  die("Too few input files");
} else {
  printf(STDERR "Interleaving %d input files\n", $numFiles);
}

my @inQual = (0) x $numFiles; # false
my @seqID = ("") x $numFiles;
my @qualID = ("") x $numFiles;
my @seq = ("") x $numFiles;
my @qual = ("") x $numFiles;
my @printable = ("") x $numFiles;
my $lineCount = 0;
my $inFile = 0;

for(my @lines = map {
  $inFile = $seqFiles{$_}; my $res = <$inFile>; $res} @seqFileOrder;
    grep {$_} @lines; # stop if all input is invalid
    @lines = map {
      $inFile = $seqFiles{$_}; my $res = <$inFile>; $res} @seqFileOrder){
  $lineCount++;
  for(my $i = 0; $i < $numFiles; $i++){
    if(!$lines[$i]){
      next;
    }
    my $line = $lines[$i];
    chomp $line; chomp $line;
    #printf(STDERR "Line $lineCount,$i: $line\n");
    if ($line =~ /^\s+$/) {
      next;
    }
    if (!$inQual[$i]) {
      if($line =~ /^(>|@)(.*)$/){
        my $newSeqID = $2;
        if($seqID[$i]){
          if($qual[$i]){
            $printable[$i] .=
              sprintf("@%s\n%s\n+\n%s\n", $seqID[$i], $seq[$i], $qual[$i]);
          } else {
            $printable[$i] .=
              sprintf(">%s\n%s\n", $seqID[$i], $seq[$i]);
          }
        }
        $qual[$i] = "";
        $seq[$i] = "";
#        printf(STDERR "Line $lineCount,$i: ".
#               "setting SeqID for $i to $newSeqID\n");
        $seqID[$i] = $newSeqID;
      } elsif ($line =~ /^\+(.*)$/) {
        if(!$seqID[$i]){
          die("[QID ] no sequence ID for $i on line $lineCount\nline: ".$line);
        }
        $inQual[$i] = 1;            # true
        $qualID[$i] = $1;
        $qual[$i] = "";
      } else {
        if(!$seqID[$i]){
          die("[SEQ ] no sequence ID for $i on line $lineCount\nline: ".$line);
        }
        $seq[$i] .= $line;
      }
    } else {
      if(!$seqID[$i]){
        die("[QUAL] no sequence ID for $i on line $lineCount\nline: ".$line);
      }
      $qual[$i] .= $line;
      if (length($qual[$i]) >= length($seq[$i])) {
        $inQual[$i] = 0;            # false
      }
    }
  } # end loop over lines from files
  if(scalar(grep {$_} @printable) == $numFiles){ # print if all can be printed
    #printf(STDERR "/----------\\\n");
    for(my $i = 0; $i < $numFiles; $i++){
      print($printable[$i]);
      $printable[$i] = "";
    }
    ## reset variables
    #printf(STDERR "\\----------/\n");
  }
} # end loop over files

foreach my $seqFileName (keys(%seqFiles)){
  close($seqFiles{$seqFileName});
}
