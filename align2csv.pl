#!/usr/bin/perl

use warnings;
use strict;

print("SeqID,RefID,Pos,Ref,Alt\n");

while(<>){
  if(/^#/){
    next;
  }
  chomp;
  my ($nameR, $nameQ, $score, $seqRef, $seqDiff, $seqQry) = split(/\t/, $_);
  my $refPos = 0;
  my $ins = ""; # false
  my $del = ""; # false
  my $indelPos = -1;
  for(my $i = 0; $i < length($seqRef); $i++){
    my $r = substr($seqRef, $i, 1);
    my $d = substr($seqDiff, $i, 1);
    my $q = substr($seqQry, $i, 1);
    if($r ne "-"){
      $refPos++;
      if(($q ne "-") && ($ins || $del)){
        printf("%s,%s,%05d,%s,%s\n", $nameQ, $nameR, $indelPos, $del, $ins);
        $ins = ""; # false
        $del = ""; # false
        $indelPos = -1;
      }
    }
    if($r eq "-"){
      $ins .= $q;
      if($indelPos == -1){
        $indelPos = $refPos+1;
      }
    }
    if($q eq "-"){
      $del .= $r;
      if($indelPos == -1){
        $indelPos = $refPos;
      }
    }
    if(!($ins || $del) && ($r ne $q)){
#    if($r ne $q){
      printf("%s,%s,%05d,%s,%s\n", $nameQ, $nameR, $refPos, $r, $q);
    }
  }
}
