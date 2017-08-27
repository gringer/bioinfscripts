#!/usr/bin/perl

## Aggregates adjacent LAST matches into a single result/match line

use warnings;
use strict;

my ($oldquery, $oldtarget, $olddir, $oldqs, $oldqe, $oldqml,
    $oldql, $oldqpct, $oldts, $oldte, $oldtml, $oldtl, $oldtpct) =
  ("") x 13;

while(<>){
  chomp;
  my ($query, $target, $dir, $qs, $qe, $qml,
    $ql, $qpct, $ts, $te, $tml, $tl, $tpct) =
      split(/,/);
  if(/^query/){
    print($_."\n");
    next;
  }
  if(($query eq $oldquery) && ($target eq $oldtarget) &&
     ($olddir eq $dir) &&
     ((($dir eq "+") && ($oldts <= $ts)) ||
      (($dir eq "-") && ($oldte >= $te)))){
    $qs = $oldqs;
    $qml += $oldqml;
    $qpct = (($oldqpct * $oldqml) + ($qpct * $qml)) / ($oldqml + $qml);
    $ts = $oldts if ($dir eq "+");
    $te = $oldte if ($dir eq "-");
    $tml += $oldtml;
    $tpct = (($oldtpct * $oldtml) + ($tpct * $tml)) / ($oldtml + $tml);
  } elsif($oldquery ne "") {
    $qpct = sprintf("%0.2f", $qpct);
    $tpct = sprintf("%0.2f", $tpct);
    print(join(",",($oldquery, $oldtarget, $olddir,
                    $oldqs, $oldqe, $oldqml, $oldql, $oldqpct,
                    $oldts, $oldte, $oldtml, $oldtl, $oldtpct))."\n");
  }
  ($oldquery, $oldtarget, $olddir, $oldqs, $oldqe, $oldqml,
   $oldql, $oldqpct, $oldts, $oldte, $oldtml, $oldtl, $oldtpct) =
     ($query, $target, $dir, $qs, $qe, $qml,
      $ql, $qpct, $ts, $te, $tml, $tl, $tpct);
}

$oldqpct = sprintf("%0.2f", $oldqpct);
$oldtpct = sprintf("%0.2f", $oldtpct);
print(join(",",($oldquery, $oldtarget, $olddir,
                $oldqs, $oldqe, $oldqml, $oldql, $oldqpct,
                $oldts, $oldte, $oldtml, $oldtl, $oldtpct))."\n");
