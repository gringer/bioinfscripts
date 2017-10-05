#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);
#use IO::Compress::Bzip2 qw(bzip2 $Bzip2Error);
#use IO::File;

my $trim = 0;
my $maxUnit = 1000;
my $sampleSize = 1000; ## number of positions to check
my $maxChunks = 10; ## maximum number of chunks to split a contig into
my $kmerLength = 13; ## number of bases in hash keys

GetOptions("trim=s" => \$trim) or
    die("Error in command line arguments");

my @rlengths = ();
my $inQual = 0; # false
my $seqID = "";
my $qualID = "";
my $seq = "";
my $qual = "";
my $buffer = "";
while(<>){
  chomp;
  chomp;
  if(!$inQual){
    if(/^(>|@)((.+?)( .*?\s*)?)$/){
      my $newSeqID = $2;
      my $newShortID = $3;
      my $len = length($seq);
      if($seqID && (length($seq) > $trim) && (length($seq) > $kmerLength)){
        my $countTotal = 0;
        my $countMax = 0;
        my $maxKmer = "";
	my @rptCounts = ();
	my %rptHash = ();
	for(my $p = 0; ($p + $kmerLength) <= $len; $p++){
	  push(@{$rptHash{substr($seq, $p, $kmerLength)}}, $p);
	}
	my $numKmers = scalar(keys(%rptHash));
	my $kmerRatio = $numKmers/($len - $kmerLength + 1);
	my @repeatedKmers = grep {scalar(@{$rptHash{$_}}) > 1} keys(%rptHash);
	my @gaps = ();
	foreach my $kmer (@repeatedKmers){
	  my @posList = @{$rptHash{$kmer}};
	  my $posCount = scalar(@posList);
	  $countTotal += $posCount;
	  if($posCount > $countMax){
	    $countMax = $posCount;
	  }
	  if(@posList){
	    push(@rptCounts, scalar(@posList));
	    #print(join(" ",@posList)."\n");
	    my $lastPos = $posList[0];
	    for(my $p = 0; $p <= $#posList; $p++){
	      my $nextPos = $posList[$p];
	      push(@gaps, $nextPos - $lastPos) if ($nextPos > $lastPos);
	      $lastPos = $nextPos;
	    }
	  }
	}
	@gaps = sort {$a <=> $b} (@gaps);
	@rptCounts = sort {$a <=> $b} (@rptCounts);
	#print(join(" ", @gaps),"\n") if (@gaps);
	my $medianGap = @gaps ? $gaps[$#gaps / 2] : 0;
	my $medianCount = @rptCounts ? $rptCounts[$#rptCounts / 2] : 0;
	my $numRepeats = scalar(@repeatedKmers);
        printf("%8d %0.3f %5d %5d %5d %5d %s\n",
               $len, $kmerRatio, scalar(@repeatedKmers), 
	       $countTotal, 
	       $medianCount,
	       $medianGap,
	       $seqID);
        push(@rlengths, $medianGap) if $medianGap;
      }
      $seq = "";
      $qual = "";
      $buffer = "";
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

my $len = length($seq);
if($seqID && (length($seq) > $trim) && (length($seq) > $kmerLength)){
  my $countTotal = 0;
  my $countMax = 0;
  my $maxKmer = "";
  my @rptCounts = ();
  my %rptHash = ();
  for(my $p = 0; ($p + $kmerLength) <= $len; $p++){
    push(@{$rptHash{substr($seq, $p, $kmerLength)}}, $p);
  }
  my $numKmers = scalar(keys(%rptHash));
  my $kmerRatio = $numKmers/($len - $kmerLength + 1);
  my @repeatedKmers = grep {scalar(@{$rptHash{$_}}) > 1} keys(%rptHash);
  my @gaps = ();
  foreach my $kmer (@repeatedKmers){
    my @posList = @{$rptHash{$kmer}};
    my $posCount = scalar(@posList);
    $countTotal += $posCount;
    if($posCount > $countMax){
      $countMax = $posCount;
    }
    if(@posList){
      push(@rptCounts, scalar(@posList));
      #print(join(" ",@posList)."\n");
      my $lastPos = $posList[0];
      for(my $p = 0; $p <= $#posList; $p++){
	my $nextPos = $posList[$p];
	push(@gaps, $nextPos - $lastPos) if ($nextPos > $lastPos);
	$lastPos = $nextPos;
      }
    }
  }
  @gaps = sort {$a <=> $b} (@gaps);
  @rptCounts = sort {$a <=> $b} (@rptCounts);
  #print(join(" ", @gaps),"\n") if (@gaps);
  my $medianGap = @gaps ? $gaps[$#gaps / 2] : 0;
  my $medianCount = @rptCounts ? $rptCounts[$#rptCounts / 2] : 0;
  my $numRepeats = scalar(@repeatedKmers);
  printf("%8d %0.3f %5d %5d %5d %5d %s\n",
	 $len, $kmerRatio, scalar(@repeatedKmers), 
	 $countTotal, 
	 $medianCount,
	 $medianGap,
	 $seqID);
  push(@rlengths, $medianGap) if $medianGap;
}

## calculate statistics
@rlengths = sort {$b <=> $a} (@rlengths);
my $sum = 0;
my @cumLengths = map {$sum += $_} (@rlengths);

printf(STDERR "Total sequences: %d\n", scalar(@rlengths));
printf(STDERR "Longest repeat: %0.3f\n", $rlengths[0]);
printf(STDERR "Shortest repeat: %0.3f\n", $rlengths[$#rlengths]);
printf(STDERR "Mean repeat length: %0.3f\n", $sum / scalar(@rlengths));
printf(STDERR "Median repeat length: %0.3f\n", $rlengths[$#rlengths / 2]);
