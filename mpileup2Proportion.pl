#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);

my $sampleName = "";
my $colourChange = 0;

GetOptions("samplename=s" => \$sampleName,
          "colour!" => \$colourChange) or
  die("Error in command line arguments");

my $assembly = "";

if($sampleName){
  printf("%-15s ", "Sample");
}
if($colourChange){
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
  if($sampleName){
    printf("%s", $sampleName);
  }
  printf("%s,%d,%d,%s,", $refName, $pos, $cov, $refAllele);
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
  printf("%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f\n",
         $pr, $pa, $pc, $pg, $pt, $pd, $pi);
}
