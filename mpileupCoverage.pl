#!/usr/bin/perl
use warnings;
use strict;

while(<>){
  chomp;
  if(!$_){
    next;
  }
  my ($refName, $pos, $refAllele, $cov, $bases, $qual, $rest) =
    split(/\t/, $_, 7);
  my @covs = ($cov);
  my @skips = (($qual =~ tr/<>//));
  while($rest){
    ($cov, $bases, $qual, $rest) = split(/\t/, $rest, 4);
    push(@covs, $cov);
    push(@skips, $skip);
  }
}
