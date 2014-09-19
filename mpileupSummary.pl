#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);

GetOptions() or
  die("Error in command line arguments");

my $assembly = "";

printf("%-15s %8s %3s %8s %s\n",
       "Assembly", "Position", "Ref", "Coverage",
       "   .   A   C   G   T   d   i");

while(<>){
  chomp;
  my @F = split(/\t/);
  if(scalar(@F) <= 4){
    next;
  }
  printf("%-15s %8d %3s %8d", $F[0], $F[1], $F[2], $F[3]);
  splice(@F,0,4);
  $_ = $F[0];
  my @inss = m/\+[0-9]+[ACGTNacgtn]+/g;
  my $i = scalar(@inss);
  s/(\+|-)[0-9]+[ACGTNacgtn]+//g;
  my $r = tr/,.//;
  my $d = tr/*//;
  my $a = tr/aA//;
  my $c = tr/cC//;
  my $g = tr/gG//;
  my $t = tr/tT//;
  printf("%4d %3d %3d %3d %3d %3d %3d\n", $r, $a, $c, $g, $t, $d, $i);
}
