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
  printf("%-20s %8s %8s %3s %s\n",
         "Assembly", "Position", "Coverage", "Ref",
         "   0    1    2    3    d    i");
} else {
  printf("%-20s %8s %8s %3s %s\n",
         "Assembly", "Position", "Coverage", "Ref",
         "   A    C    G    T    d    i");
}

while(<>){
#  print(STDERR $_);
  chomp;
  my ($refName, $pos, $refAllele, $cov, $bases, $rest) = split(/\t/, $_, 6);
  if($sampleName){
    printf("%-15s ", $sampleName);
  }
  printf("%-20s %8d %8d %3s", $refName, $pos, $cov, $refAllele);
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
  my $total = $i+$r+$d+$a+$c+$g+$t;
  if($refAllele eq "A"){
    $a = $r;
  } elsif($refAllele eq "C"){
    $c = $r;
  } elsif($refAllele eq "G"){
    $g = $r;
  } elsif($refAllele eq "T"){
    $t = $r;
  }
  # was previously $coverage, not $total
  ($r, $i, $d, $a, $c, $g, $t) = map {$_ / $total}
      ($r, $i, $d, $a, $c, $g, $t);
  printf(" %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f\n",
         $a, $c, $g, $t, $d, $i);
}
