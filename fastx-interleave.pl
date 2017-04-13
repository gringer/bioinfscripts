#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $singleLine = 0;
my $minLength = 0;
my $maxCount = 0; # use reservoir sampling to randomise reads

my %seqFiles = ();
my @seqFileOrder = ();

GetOptions("singleLine!" => \$singleLine, "minLength=i" => \$minLength,
           "count=i" => \$maxCount) or
  die("Error in command line arguments");

while(@ARGV){
  my $argument = shift(@ARGV);
  if(-e $argument){ # file existence check
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

if($minLength){
  printf(STDERR "Only outputting read groups for reads of at least %d bases\n", $minLength);
}

if($maxCount){
  printf(STDERR "Reservoir sampling to output at most %d read groups:",
        $maxCount);
}

my @inQual = (0) x $numFiles; # false
my @seqID = ("") x $numFiles;
my @qualID = ("") x $numFiles;
my @seq = ("") x $numFiles;
my @qual = ("") x $numFiles;
my @printable = ("") x $numFiles;
my $lineCount = 0;
my $inFile = 0;
my $recordsRead = 0;
my $dotsPrinted = 0;

my @printReservoir = ();

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
          if($printable[$i]){
            print(STDERR "Warning: double print for file $i\n");
          }
          if($qual[$i]){
            $printable[$i] .=
              sprintf("@%s\n%s\n+\n%s\n", $seqID[$i], $seq[$i], $qual[$i]);
          } else {
            $printable[$i] .=
              sprintf(">%s\n%s\n", $seqID[$i], $seq[$i]);
          }
          if(length($seq[$i]) < $minLength){
            $printable[$i] = "#";
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
    if(grep {/^#/} @printable){
      ## reset / clear all if any are non-printable
      for(my $i = 0; $i < $numFiles; $i++){
        $printable[$i] = "";
      }
    } else {
      if($maxCount && ($recordsRead % 10000 == 0)){
        if($dotsPrinted % 50 == 0){
          if($recordsRead > 1000){
            printf(STDERR " (%d read groups processed)", $recordsRead);
          }
          printf(STDERR "\n  ");
        }
        print(STDERR ".");
        $dotsPrinted++;
      }
      $recordsRead++;
      my $linesToAdd = "";
      for(my $i = 0; $i < $numFiles; $i++){
        $linesToAdd .= $printable[$i];
        $printable[$i] = "";
      }
      if(!$maxCount){
        print($linesToAdd);
      } elsif($maxCount >= $recordsRead){
        push(@printReservoir, $linesToAdd);
      } else {
        my $swapPos = rand($recordsRead);
        if($swapPos < $maxCount){
          $printReservoir[$swapPos] = $linesToAdd;
        }
      }
    }
  }
} # end loop over files

foreach my $seqFileName (keys(%seqFiles)){
  close($seqFiles{$seqFileName});
}

if($maxCount){
  printf(STDERR "\ndone (%d read groups processed)\n", $recordsRead);
  print(join("",@printReservoir));
}
