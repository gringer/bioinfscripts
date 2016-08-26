#!/usr/bin/perl

use warnings;
use strict;

my $model = "";
my $cumLength = 0;
my $meanSig = 0;
my $bpPos = 0;
my $lastMove = 0;
my $lastStart = 0;

while(<>){
  chomp;
  my @F = split(",", $_);
  if($F[5] !~ /[0-9]/){
    print("channel,mux,read,model_state,move,bppos,mean,start,length\n");
    next;
  }
  ## convert start/length to samples, assuming 4kHz sampling
  @F[(5,7)] = map({$_ * 4000} @F[(5,7)]);
  if($F[8] ne $model){
    if($model){
      print(join(",",(@F[(1,2,3)],$model,$lastMove,$bpPos,$meanSig,$lastStart,$cumLength)),"\n");
    }
    $model = $F[8];
    $lastMove = $F[9];
    $lastStart = $F[5];
    $bpPos += $F[9];
    $cumLength = 0;
    $meanSig = 0;
  }
  $meanSig = ($meanSig * $cumLength + $F[4] * $F[7]) / ($cumLength + $F[7]);
  $cumLength = $cumLength + $F[7];
  ##print(join(";",@F[(1,2,3,8,9,4,5,7)]),"\n");
}
