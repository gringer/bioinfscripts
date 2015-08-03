#!/usr/bin/perl
use warnings;
use strict;

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
my $colourChange = 0;
my $minCoverage = 0;
my $writeCounts = 0;

GetOptions("mincoverage=i" => \$minCoverage,
	   "samplename=s" => \$sampleName,
           "colour!" => \$colourChange,
           "counts!" => \$writeCounts) or
  die("Error in command line arguments");

my $assembly = "";

if($sampleName){
  printf("%-15s ", "Sample");
}
if($colourChange){
  warn("Warning: Colour change calculations are not yet properly implemented");
  printf("%s,%s,%s,%s,%s\n",
         "Assembly", "Position", "Coverage", "cR",
         "0,1,2,3,d,i");
} else {
  printf("%s,%s,%s,%s,%s\n",
         "Assembly", "Position", "Coverage", "cR",
         "pR,A,C,G,T,d,i");
}

while(<>){
#  print(STDERR $_);
  chomp;
  my ($refName, $pos, $refAllele, $cov, $bases, $rest) = split(/\t/, $_, 6);
  $_ = uc($bases);
  my $i = scalar(m/\+[0-9]+[ACGTNacgtn]+/g);
  s/\^.//g;
  s/(\+|-)[0-9]+[ACGTNacgtn]+//g;
  my $r = tr/,.//;
  my $d = tr/*//;
  my $a = tr/aA//;
  my $c = tr/cC//;
  my $g = tr/gG//;
  my $t = tr/tT//;
  my ($pr, $pi, $pd, $pa, $pc, $pg, $pt) = (0, 0, 0, 0, 0, 0, 0);
  my $total = $i+$r+$d+$a+$c+$g+$t;
  # if($refAllele eq "A"){
  #   $a = $r;
  # } elsif($refAllele eq "C"){
  #   $c = $r;
  # } elsif($refAllele eq "G"){
  #   $g = $r;
  # } elsif($refAllele eq "T"){
  #   $t = $r;
  # }
  # was previously $coverage, not $total
  if($total > 0){
    ($pr, $pi, $pd, $pa, $pc, $pg, $pt) = map {$_ / $total}
      ($r, $i, $d, $a, $c, $g, $t);
  }
  if($cov > $minCoverage){
      if($sampleName){
	  printf("%s,", $sampleName);
      }
      printf("%s,%d,%d,%s,", $refName, $pos, $cov, $refAllele);
      if($writeCounts){
        printf("%d,%d,%d,%d,%d,%d,%d\n",
               $r, $a, $c, $g, $t, $d, $i);
      } else {
        printf("%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f\n",
               $pr, $pa, $pc, $pg, $pt, $pd, $pi);
      }
  }
}
