#!/usr/bin/perl
use warnings;
use strict;

# readstomper.pl -- generate base-packed proportion statistics for
#   output from 'samtools mpileup'
#
# example use:
# samtools mpileup -B -Q 0 -f circ-Nb-ec3-mtDNA.fasta LAST_OL_132394_vs_mtDNA.bam | \
#   /bioinf/scripts/readstomper.pl > LAST_OLmpileupProp_132394_vs_mtDNA.txt.gz
#
# example output:
# Assembly,Position,Coverage,ref,cR,pR,A,C,G,T,d,i,InsMode
# Nb_mtDNA,11456,86,T,79.00,0.92,0.0000,0.0465,0.0000,0.0000,0.03,0.01,A;1.00
# Nb_mtDNA,11457,86,A,81.00,0.94,0.0000,0.0000,0.0233,0.0000,0.03,0.02,TA;0.50
# Nb_mtDNA,11458,86,G,79.00,0.92,0.0349,0.0233,0.0000,0.0116,0.01,0.00
# Nb_mtDNA,11459,86,A,78.00,0.92,0.0000,0.0000,0.0235,0.0471,0.01,0.00
# Nb_mtDNA,11460,86,T,5.00,0.06,0.0233,0.8140,0.0233,0.0000,0.08,0.00
# Nb_mtDNA,11461,86,G,57.00,0.66,0.1512,0.1163,0.0000,0.0465,0.02,0.01,AATAACAACACGTAACCGA;1.00
# Nb_mtDNA,11462,86,G,78.00,0.91,0.0349,0.0349,0.0000,0.0116,0.01,0.00
# Nb_mtDNA,11463,86,G,66.00,0.77,0.0465,0.0698,0.0000,0.0116,0.10,0.01,GCAATA;1.00
# Nb_mtDNA,11464,86,T,74.00,0.86,0.0000,0.0814,0.0233,0.0000,0.03,0.00

use Getopt::Long qw(:config auto_help pass_through);

my $sampleName = "";
my $minCoverage = 0;
my $maxCoverage = 10**9;
my $writeCounts = 0;
my $writeConsensus = 0;
my $consThresholdCov = 5;
my $insThresholdCov = 5;
my $deletionSens = 0.15;
my $edit = 1;

GetOptions("mincoverage=i" => \$minCoverage,
           "transposoncoverage=i" => \$maxCoverage,
	   "samplename=s" => \$sampleName,
           "deletionsensitivity=s" => \$deletionSens,
           "insertioncoverage=i" => \$insThresholdCov,
           "fasta=s" => \$writeConsensus,
           "counts!" => \$writeCounts,
           "edit!" => \$edit ) or
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
  if($sampleName){
    printf(STDERR "%s,", "Sample");
  }
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
my $highCovStart = 0;

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
        } elsif($oldRefName){
          print("\n");
        }
      }
      $seqChanged = 0;
      print(">${refName}\n"); ## write new sequence header
    }
    $oldRefName = $refName;
    $lastBase = 0;
  }
  if($writeConsensus){
    if(++$lastBase < $pos){
      ## print sequence from the intervening gap
      print(substr($refSeqs{$refName}, ($lastBase), ($pos - $lastBase)));
    }
    $lastBase = $pos;
    if($cov >= $maxCoverage){
      print("N"); ## mask out likely transposon sequences
      next;
    }
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
  my $dc = $deletions{$pos}?$deletions{$pos}:0;
  delete($deletions{$pos});
  my $ac = tr/aA//;
  my $cc = tr/cC//;
  my $gc = tr/gG//;
  my $tc = tr/tT//;
  my $nc = tr/nN//;
  my ($pr, $pi, $pd, $pa, $pc, $pg, $pt, $pn) = (0, 0, 0, 0, 0, 0, 0);
  # note: insertions don't count towards total coverage,
  #       they are additional features attached to a read base
  my $total = $rc+$dc+$ac+$cc+$gc+$tc+$nc;
  if($total > 0){
    ($pr, $pi, $pd, $pa, $pc, $pg, $pt, $pn) = map {$_ / $total}
      ($rc, $ic, $dc, $ac, $cc, $gc, $tc, $nc);
  }
  if($writeConsensus){
    ## determine consensus allele
    my $consAllele = $refAllele;
    if($edit && ($total > $consThresholdCov) &&
       (($pr < 0.5) || ($pd > $deletionSens) || (($pi > 0.5) && ($maxInserts > $insThresholdCov)))){
      $seqChanged = 0;
      my %consCounts =
        (r => $rc,
         A => $ac,
         C => $cc,
         G => $gc,
         T => $tc,
         N => $nc,
         d => $dc);
      my @sortedAlleles = sort {$consCounts{$b} <=> $consCounts{$a}} keys(%consCounts);
      if(($sortedAlleles[0] eq "d") || ($pd > $deletionSens)){
        $consAllele = "";
      } elsif($sortedAlleles[0] ne "r"){
        $consAllele = $sortedAlleles[0];
      }
      if($consAllele ne $refAllele){
        if($sampleName){
          printf(STDERR "%s,", $sampleName);
        }
        printf(STDERR "%s,%d,%d,%s,", $refName, $pos, $cov, $refAllele);
        printf(STDERR "%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f",
               $rc, $pr, $pa, $pc, $pg, $pt, $pd, $pi);
        if($maxInserts > $insThresholdCov){
          printf(STDERR ",%s;%0.2f", $maxInsertSeq, $maxInserts / $ic);
          printf(STDERR ",[Insert %s]", $maxInsertSeq);
        } else {
          printf(STDERR ",");
          printf(STDERR ",[%s -> %s]", $refAllele, $consAllele);
        }
        print(STDERR "\n");
      }
    }
    #print(",$consAllele");
    print($consAllele); ## print current reference allele
    if(($pi > 0.5) && ($maxInserts > $insThresholdCov)){ ## more than 50% of reads suggest insertion
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
