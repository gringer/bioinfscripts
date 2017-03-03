#!/usr/bin/perl

use warnings;
use strict;

my $model = "";
my $cumLength = 0;
my $cumTotalLength = 0;
my $meanSig = 0;
my $bpPos = 0;
my $lastMove = 0;
my $lastStart = 0;

my %colnums = ();

while(<>){
  chomp;
  my @F = split(",", $_);
  if($F[5] !~ /[0-9]/){
    my $colNum = 0;
    grep {$colnums{$_} = $colNum++} @F;
    print("channel,mux,read,model_state,move,bppos,mean,start,startSecs,length\n");
    next;
  }
  ## convert start/length to samples, assuming 4kHz sampling
  @F[($colnums{start},$colnums{length})] =
    map({$_ * $F[$colnums{sampleRate}]} @F[($colnums{start},$colnums{length})]);
  if($F[$colnums{model_state}] ne $model){
    if($model){
      print(join(",",(@F[($colnums{runID},$colnums{channel},$colnums{mux})],
                      $model,$lastMove,$bpPos,sprintf('%0.2f',$meanSig),
                      sprintf('%.0f', $lastStart),
                      sprintf('%0.4f', $lastStart/$F[$colnums{sampleRate}]),
                      sprintf('%.0f',$cumLength)
                     )),"\n");
    }
    $model = $F[$colnums{model_state}];
    $lastMove = $F[$colnums{move}];
    $lastStart = $F[$colnums{start}];
    $bpPos += $F[$colnums{move}];
    $cumLength = 0;
    $meanSig = 0;
  }
  $meanSig = ($meanSig * $cumLength + $F[$colnums{mean}] * $F[$colnums{length}]) / ($cumLength + $F[$colnums{length}]);
  $cumLength = $cumLength + $F[$colnums{length}];
  $cumLength = int($cumLength);
  $cumTotalLength += $cumLength;
}
