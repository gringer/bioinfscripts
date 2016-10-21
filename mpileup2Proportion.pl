#!/usr/bin/perl
use warnings;
use strict;
use Data::Dumper::Simple;

# mpileup2Proportion.pl -- generate base/INDEL proportion statistics for
#   output from 'samtools mpileup'
#
# example use:
# samtools mpileup -C50 -Q0 -e 20 -o 40 -f ref.fasta input.bam |
#   mpileup2Proportion.pl > output.prop.csv
#
# example output:
# Assembly,Position,Ref,Coverage,cR,pR,A,C,G,T,d,i
# mmusMT_PCR1,15564,G,89,85,95.5,0.0,2.2,95.5,1.1,1.1,2.2
# mmusMT_PCR1,15565,A,89,85,95.5,95.5,1.1,3.4,0.0,0.0,7.9
# mmusMT_PCR1,15566,A,89,83,93.3,93.3,0.0,2.2,1.1,0.0,0.0
# mmusMT_PCR1,15567,T,89,73,82.0,4.5,0.0,5.6,82.0,7.9,1.1
# mmusMT_PCR1,15568,A,89,71,79.8,79.8,3.4,2.2,1.1,13.5,2.2


use Getopt::Long qw(:config auto_help pass_through);

my $sampleName = "";
my $minCoverage = 0;
my $writeCounts = 0;
my $writeConsensus = 0;
my $consThresholdCov = 5;
my $deletionSens = 0.15;

GetOptions("mincoverage=i" => \$minCoverage,
	   "samplename=s" => \$sampleName,
           "deletionsensitivity=s" => \$deletionSens,
           "fasta=s" => \$writeConsensus,
           "counts!" => \$writeCounts) or
  die("Error in command line arguments");

my $assembly = "";

if($writeConsensus && !$sampleName){
  $sampleName = "samplename";
}

if(!$writeConsensus){
  if($sampleName){
    printf("%s,", "Sample");
  }
  printf("%s,%s,%s,%s,%s,%s\n",
         "Assembly", "Position", "Coverage", "ref", "cR",
         "pR,A,C,G,T,d,i,InsMode");
} else {
  printf(STDERR "%s,%s,%s,%s,%s,%s\n",
         "Assembly", "Position", "Coverage", "ref", "cR",
         "pR,A,C,G,T,d,i,InsMode,Variant");
}

my %refSeqs = ();
if($writeConsensus){
  if(!(-f $writeConsensus)){
    die("ERROR: Reference fasta sequence '${writeConsensus}' does not exist");
  }
  open(my $refFile, "<", $writeConsensus);
  my $seqID = "";
  while(<$refFile>){ # parse FASTA file
    chomp;
    if(/^>(.+?)( .*)?$/){
      $seqID = $1;
      $refSeqs{$seqID} = "*"; # placeholder to simplify substr
    } else {
      $refSeqs{$seqID} .= $_;
    }
  }
  close($refFile);
}

my %deletions = ();
my $oldRefName = "";
my $lastBase = 0;
my $seqChanged = 0;

while(<>){
  chomp;
  my ($refName, $pos, $refAllele, $cov, $bases, $rest) = split(/\t/, $_, 6);
  my $observedDiff = 0; # false
  if($cov < $minCoverage){
    next;
  }
  if($oldRefName ne $refName){
    if($writeConsensus){  ## complete old sequence (if any)
      if(!$refSeqs{$oldRefName}){
        print(STDERR "Warning: reference '${oldRefName}' not found\n") unless ($oldRefName eq "");
      } else {
        if($oldRefName && (length($refSeqs{$oldRefName}) < $lastBase)){
          print(substr($refSeqs{$oldRefName}, ($lastBase+1))."\n");
        } elsif($oldRefName && $seqChanged){
          print("\n");
        }
        $seqChanged = 0;
      }
      print(">${refName}\n"); ## write new sequence header
    }
    $oldRefName = $refName;
    $lastBase = 0;
  }
  if($writeConsensus){
    if(++$lastBase < $pos){  ## print sequence from the intervening gap
      print(substr($refSeqs{$refName}, ($lastBase), ($pos - $lastBase)));
    }
    $lastBase = $pos;
  }
  $_ = uc($bases);
  ## process insertions
  my %insertCounts = ();
  my $ic = 0;
  my $maxInserts = 0;
  my $maxInsertSeq = "";
  while(s/\+[0-9]+([ACGTNacgtn]+)//){
    my $insertSeq = $1;
    $insertCounts{$insertSeq}++;
    if($insertCounts{$insertSeq} > $maxInserts){
      $maxInsertSeq = $insertSeq;
      $maxInserts = $insertCounts{$insertSeq};
    }
    $ic++;
  }
  s/\^.//g; # remove "start of read" + "read mapping quality" indicators
  ## process deletions
  while(s/-([0-9]+)[ACGTNacgtn]+//){
      ## deletions are a special case, because *all* deletions are replaced with a single '*'
      ## this is worked around by ignoring '*' and parsing the delete string for future positions
      my $delSize = $1;
      for(my $i = 1; $i <= $delSize; $i++){
	  $deletions{$pos+$i}++;
      }
  }
  ## remove stray insertions and deletions (which shouldn't exist...)
  s/(\+|-)[0-9]+[ACGTNacgtn]+//g;
  my $rc = tr/,.//;
  #my $dc = tr/*//;
  my $dc = $deletions{$pos}?$deletions{$pos}:0;
  delete($deletions{$pos});
  my $ac = tr/aA//;
  my $cc = tr/cC//;
  my $gc = tr/gG//;
  my $tc = tr/tT//;
  my ($pr, $pi, $pd, $pa, $pc, $pg, $pt) = (0, 0, 0, 0, 0, 0, 0);
  my ($p0, $p1, $p2, $p3) = (0, 0, 0, 0);
  # note: insertions don't count towards total coverage,
  #       they are additional features attached to a read base
  my $total = $rc+$dc+$ac+$cc+$gc+$tc;
  # if($refAllele eq "A"){
  #   $ac = $rc;
  # } elsif($refAllele eq "C"){
  #   $cc = $rc;
  # } elsif($refAllele eq "G"){
  #   $gc = $rc;
  # } elsif($refAllele eq "T"){
  #   $tc = $rc;
  # }
  # was previously $coverage, not $total
  if($total > 0){
    ($pr, $pi, $pd, $pa, $pc, $pg, $pt) = map {$_ / $total}
      ($rc, $ic, $dc, $ac, $cc, $gc, $tc);
  }
  if($writeConsensus){
    ## determine consensus allele
    my $consAllele = $refAllele;
    if(($total > $consThresholdCov) &&
       (($pr < 0.5) || ($pd > $deletionSens) || ($pi > 0.5))){
      $seqChanged = 0;
      my %consCounts =
        (r => $rc,
         A => $ac,
         C => $cc,
         G => $gc,
         T => $tc,
         d => $dc);
      my @sortedAlleles = sort {$consCounts{$b} <=> $consCounts{$a}} keys(%consCounts);
      if(($sortedAlleles[0] eq "d") || ($pd > $deletionSens)){
        $consAllele = "";
      } elsif($sortedAlleles[0] ne "r"){
        $consAllele = $sortedAlleles[0];
      }
      printf(STDERR "%s,%d,%d,%s,", $refName, $pos, $cov, $refAllele);
      printf(STDERR "%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f",
             $rc, $pr, $pa, $pc, $pg, $pt, $pd, $pi);
      if($maxInsertSeq){
        printf(STDERR ",%s;%0.2f", $maxInsertSeq, $maxInserts / $ic);
        printf(STDERR ",[Insert %s]", $maxInsertSeq);
      } else {
        printf(STDERR ",");
        printf(STDERR ",[%s -> %s]", $refAllele, $consAllele);
      }
      print(STDERR "\n");
    }
    #print(",$consAllele");
    print($consAllele); ## print current reference allele
    if($pi > 0.5){ ## more than 50% of reads suggest insertion
      print(uc($maxInsertSeq));
    }
  } else {
    if($sampleName){
      printf("%s,", $sampleName);
    }
    printf("%s,%d,%d,%s,", $refName, $pos, $cov, $refAllele);
    if($writeCounts){
      printf("%d,%0.2f,%d,%d,%d,%d,%d,%d",
             $rc, $pr, $ac, $cc, $gc, $tc, $dc, $ic);
      if($maxInsertSeq){
        printf(",%s;%d", $maxInsertSeq, $maxInserts);
      }
    } else {
      printf("%0.2f,%0.2f,%0.4f,%0.4f,%0.4f,%0.4f,%0.2f,%0.2f",
               $rc, $pr, $pa, $pc, $pg, $pt, $pd, $pi);
      if($maxInsertSeq){
        printf(",%s;%0.2f", $maxInsertSeq, $maxInserts / $ic);
      }
    }
    print("\n");
  }
}

if($writeConsensus){
  if($oldRefName){
    print(substr($refSeqs{$oldRefName}, ($lastBase+1))."\n");
  }
}
