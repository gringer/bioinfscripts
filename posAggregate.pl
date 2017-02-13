#!/usr/bin/perl

## Aggregates adjacent LAST matches into a single result/match line

use warnings;
use strict;

my ($oldquery, $oldtarget, $olddir, $oldqs, $oldqe, $oldqml,
    $oldql, $oldts, $oldtml, $oldtl) =
  ("") x 10;

while(<>){
  chomp;
  my ($query, $target, $dir, $qs, $qe, $qml, $ql, $ts, $tml, $tl) =
    split(/\s+/);
  if(($query eq $oldquery) && ($target eq $oldtarget) &&
     ($olddir eq $dir) && ($oldts <= $ts)){
    $qs = $oldqs;
    $qml += $oldqml;
    $ts = $oldts;
    $tml += $oldtml;
  } elsif($oldquery ne "") {
    printf("%-55s %-15s %-3s %-4s %-4s %-4s %-4s %-4s %-4s %-4s\n",
           $oldquery, $oldtarget, $olddir, $oldqs, $oldqe, $oldqml,
           $oldql, $oldts, $oldtml, $oldtl);
  }
  ($oldquery, $oldtarget, $olddir, $oldqs, $oldqe, $oldqml,
   $oldql, $oldts, $oldtml, $oldtl) =
     ($query, $target, $dir, $qs, $qe, $qml, $ql, $ts, $tml, $tl);
}

printf("%-55s %-15s %-3s %-4s %-4s %-4s %-4s %-4s %-4s %-4s\n",
       $oldquery, $oldtarget, $olddir, $oldqs, $oldqe, $oldqml,
       $oldql, $oldts, $oldtml, $oldtl);
