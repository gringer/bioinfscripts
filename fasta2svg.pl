#!/usr/bin/perl

# fasta2svg.pl -- converts a fasta sequence to a colour-coded
# sequence.
#
# Homopolymer regions of a specified length are used to define colour
# region borders. The homopolymers themselves are coloured using the
# standard electropherogram colours for bases:
#
# A -- Green  C -- Blue
# G -- Yellow T -- Red
#
# Regions between homopolymer sequences are coded in the RGB space
# with each of the three components defined based on the proportion of
# the three possible binary divisions of bases:
#
# R component: M proportion (i.e. AC ratio) [alternate colour cyan]
# G component: Y proportion (i.e. TC ratio) [alternate colour magenta]
# B component: S proportion (i.e. GC ratio) [alternate colour yellow]
#
# These colours can be presented as-is, or scaled based on the
# observed range (or distribution) of proportions. The effect of
# proportion scaling would be an increase in the colour saturation of
# regions between homopolymer sequences.
#
# Copyright 2015, David Eccles (gringer) <bioinformatics@gringene.org>
#
# Permission to use, copy, modify, and/or distribute this software for
# any purpose with or without fee is hereby granted. The software is
# provided "as is" and the author disclaims all warranties with regard
# to this software including all implied warranties of merchantability
# and fitness. The parties responsible for running the code are solely
# liable for the consequences of code excecution.

use warnings;
use strict;

use Getopt::Long qw(:config auto_version auto_help pass_through); # for option parsing

my $hpLength = 1;
my @pScale = [0,1,0,1,0,1];

GetOptions('hplength=i' => \$hpLength, 'pscale=f{6}' => \@pScale) or
  die("Error in command line arguments");

my $seq = "";
my $seqCount = 1;
my $seqID = "";
my @newPS = (0.5,0.5,0.5,0.5,0.5,0.5);

print("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
print("<svg xmlns:svg=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n");


while(<>){
  chomp;
  if(/^>(.+)$/){
    my $newID = $1;
    if($seq){
      my @tPS = drawSeq($seq, $seqID, $seqCount, $hpLength);
      #print(STDERR join(":",@tPS)."\n");
      foreach my $i (0..2){
        $newPS[$i*2] = $tPS[$i*2] if($tPS[$i*2] < $newPS[$i*2]);
        $newPS[$i*2+1] = $tPS[$i*2+1] if($tPS[$i*2+1] > $newPS[$i*2+1]);
      }
      $seqCount++;
    }
    $seq = "";
    $seqID = $newID;
  } else {
    $seq .= $_;
  }
}

if($seq){
  my @tPS = drawSeq($seq, $seqID, $seqCount, $hpLength);
  foreach my $i (0..2){
    $newPS[$i*2] = $tPS[$i*2] if($tPS[$i*2] < $newPS[$i*2]);
    $newPS[$i*2+1] = $tPS[$i*2+1] if($tPS[$i*2+1] > $newPS[$i*2+1]);
  }
}

print("</svg>\n");


sub comp{
  my ($tSeq) = @_;
  $tSeq =~ tr/ACGT/TGCA/;
  return($tSeq);
}

sub updateRange{
  my ($cr, $newCr) = @_;
  #print(STDERR join(":",@{$cr})." ".join(":",@{$newCr})."\n");
  foreach my $i (0..2){
    ${$cr}[$i*2] = ${$newCr}[$i]
      if((${$newCr}[$i] > 0) && (${$newCr}[$i] < ${$cr}[$i*2]));
    ${$cr}[$i*2+1] = ${$newCr}[$i]
      if((${$newCr}[$i] < 1) && (${$newCr}[$i] > ${$cr}[$i*2+1]));
  }
}

sub contentColour{
  my ($tSeq) = @_;
  my $total = ($tSeq =~ tr/AaCcGgTt//);
  if($total == 0){
    return(("#000000", 0.5, 0.5, 0.5));
  }
  my $a = ($tSeq =~ tr/Aa//);
  my $c = ($tSeq =~ tr/Cc//);
  my $g = ($tSeq =~ tr/Gg//);
  my $t = ($tSeq =~ tr/Tt//);
  my ($S, $Y, $M) = map {$_ / $total} (($c+$g), ($c+$t), ($c+$a));
  my ($Sc, $Yc, $Mc) = map {$_ * 255} ($S, $Y, $M);
  return((sprintf("#%02X%02X%02X", $Mc, $Yc, $Sc), $S, $Y, $M));
}

sub drawSeq{
  my ($tSeq, $tSeqID, $tSeqCount, $hl) = @_;
  my $ccRange = [0.5,0.5,0.5,0.5,0.5,0.5];
  my %hpCols = ( A => "#00FF00", C => "#0000FF", G=> "#FFFF00", T=> "#FF0000");
  #$tSeq = substr($tSeq,0,200); # only show first 100 bases for testing purposes
  my $tSeqC = comp($tSeq);
  printf(" <g id=\"%s\" stroke-width=\"%s\">\n", $tSeqID, $hl*2);
  my $sf = $hl*4+2.5; ## size adjustment factor
  my $ofy = $hl+0.25; ## y offset for fwd/comp sequence
  my $xp = $sf;
  printf("  <g id=\"%s_FWD\">\n", $tSeqID);
  ## surrounding black rectangle so that white areas don't look odd
  printf("   <rect fill=\"black\" x=\"%s\" y=\"%s\" ".
         "width=\"%s\" height=\"%s\" />\n",
         $xp - 0.5, $tSeqCount*$sf-(($sf-1)/2), length($tSeq)+1, $sf-1);
  while($tSeq =~ s/^(.*?)(A{$hl,}|C{$hl,}|G{$hl,}|T{$hl,})//){
    my $preSeq = $1;
    my $hpSeq = $2;
    my $hpBase = substr($hpSeq,0,1);
    if($preSeq){
      ## display sequence prior to homopolymer
      my ($col, @CC) = contentColour($preSeq);
      if(length($preSeq) > 5){
        updateRange($ccRange, \@CC);
      }
      printf("   <path stroke=\"%s\" d=\"M%s,%s l%s,0\" />\n",
             $col, $xp, $tSeqCount*$sf-$ofy, length($preSeq));
    }
    $xp += length($preSeq);
    ## display homopolymer sequence
    printf("   <path stroke=\"%s\" d=\"M%s,%s l%s,0\" />\n",
           $hpCols{$hpBase}, $xp, $tSeqCount*$sf-$ofy, length($hpSeq));
    $xp += length($hpSeq);
  }
  if($tSeq){
    my ($col, @CC) = contentColour($tSeq);
    if(length($tSeq) > 5){
      updateRange($ccRange, \@CC);
    }
    ## display sequence after last homopolymer
    printf("   <path stroke=\"%s\" d=\"M%s,%s l%s,0\" />\n",
           $col, $xp, $tSeqCount*$sf-$ofy, length($tSeq));
  }
  printf("  </g>\n");
  $xp = $sf;
  ## complement sequence
  printf("  <g id=\"%s_REV\">\n", $tSeqID);
  while($tSeqC =~ s/^(.*?)(A{$hl,}|C{$hl,}|G{$hl,}|T{$hl,})//){
    my $preSeq = $1;
    my $hpSeq = $2;
    my $hpBase = substr($hpSeq,0,1);
    if($preSeq){
      ## display sequence prior to homopolymer
      my ($col, @CC) = contentColour($preSeq);
      if(length($preSeq) > 5){
        updateRange($ccRange, \@CC);
      }
      printf("   <path stroke=\"%s\" d=\"M%s,%s l%s,0\" />\n",
             $col, $xp, $tSeqCount*$sf+$ofy, length($preSeq));
    }
    $xp += length($preSeq);
    ## display homopolymer sequence
    printf("   <path stroke=\"%s\" d=\"M%s,%s l%s,0\" />\n",
           $hpCols{$hpBase}, $xp, $tSeqCount*$sf+$ofy, length($hpSeq));
    $xp += length($hpSeq);
  }
  if($tSeqC){
    my ($col, @CC) = contentColour($tSeqC);
    if(length($tSeqC) > 5){
      updateRange($ccRange, \@CC);
    }
    ## display sequence after last homopolymer
    printf("   <path stroke=\"%s\" d=\"M%s,%s l%s,0\" />\n",
           $col, $xp, $tSeqCount*$sf+$ofy, length($tSeqC));
  }
  printf("  </g>\n");
  printf(" </g>\n");
  return(@{$ccRange});
}
