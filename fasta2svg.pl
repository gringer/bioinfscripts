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

GetOptions('hplength=i' => \$hpLength) or
  die("Error in command line arguments");

my $seq = "";
my $seqCount = 1;
my $seqID = "";

print("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
print("<svg xmlns:svg=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n");


while(<>){
  chomp;
  if(/^>(.+)$/){
    my $newID = $1;
    if($seq){
      drawSeq($seq, $seqID, $seqCount, $hpLength);
      $seqCount++;
    }
    $seq = "";
    $seqID = $newID;
  } else {
    $seq .= $_;
  }
}

if($seq){
  drawSeq($seq, $seqID, $seqCount, $hpLength);
}

print("</svg>\n");

sub comp{
  my ($tSeq) = @_;
  $tSeq =~ tr/ACGT/TGCA/;
  return($tSeq);
}

sub contentColour{
  my ($tSeq) = @_;
  my $total = ($tSeq =~ tr/AaCcGgTt//);
  my $a = ($tSeq =~ tr/Aa//);
  my $c = ($tSeq =~ tr/Cc//);
  my $g = ($tSeq =~ tr/Gg//);
  my $t = ($tSeq =~ tr/Tt//);
  my $S = (($c+$g) * 255) / $total;
  my $Y = (($c+$t) * 255) / $total;
  my $M = (($c+$a) * 255) / $total;
  return(sprintf("#%02X%02X%02X", $M, $Y, $S))
}

sub drawSeq{
  my ($tSeq, $tSeqID, $tSeqCount, $hl) = @_;
  my %hpCols = ( A => "#00FF00", C => "#0000FF", G=> "#FFFF00", T=> "#FF0000");
  #$tSeq = substr($tSeq,0,200); # only show first 100 bases for testing purposes
  my $tSeqC = comp($tSeq);
  printf(" <g id=\"%s\" stroke-width=\"%s\">\n", $tSeqID, $hl);
  my $sf = 5*$hl; ## size adjustment factor
  my $sfy = $hl/2+0.5; ## size adjustment factor (y)
  my $xp = $sf;
  printf("  <g id=\"%s_FWD\">\n", $tSeqID);
  printf("   <rect fill=\"black\" x=\"%s\" y=\"%s\" width=\"%s\" height=\"%s\" />\n",
         $xp - 0.5, $tSeqCount*$sf-$hl-1, length($tSeq)+1, 2*$hl+2);
  while($tSeq =~ s/^(.*?)(A{$hl,}|C{$hl,}|G{$hl,}|T{$hl,})//){
    my $preSeq = $1;
    my $hpSeq = $2;
    my $hpBase = substr($hpSeq,0,1);
    if($preSeq){
      printf("   <path stroke=\"%s\" d=\"M%s,%s l%s,0\" />\n",
             contentColour($preSeq), $xp, $tSeqCount*$sf-$sfy, length($preSeq));
    }
    $xp += length($preSeq);
    printf("   <path stroke=\"%s\" d=\"M%s,%s l%s,0\" />\n",
           $hpCols{$hpBase}, $xp, $tSeqCount*$sf-$sfy, length($hpSeq));
    $xp += length($hpSeq);
  }
  if($tSeq){
    printf("   <path stroke=\"%s\" d=\"M%s,%s l%s,0\" />\n",
           contentColour($tSeq), $xp, $tSeqCount*$sf-$sfy, length($tSeq));
  }
  printf("  </g>\n");
  $xp = $sf;
  printf("  <g id=\"%s_REV\">\n", $tSeqID);
  #print(STDERR $tSeqC."\n");
  while($tSeqC =~ s/^(.*?)(A{$hl,}|C{$hl,}|G{$hl,}|T{$hl,})//){
    my $preSeq = $1;
    my $hpSeq = $2;
    my $hpBase = substr($hpSeq,0,1);
    if($preSeq){
      printf("   <path stroke=\"%s\" d=\"M%s,%s l%s,0\" />\n",
             contentColour($preSeq), $xp, $tSeqCount*$sf+$sfy, length($preSeq));
    }
    $xp += length($preSeq);
    printf("   <path stroke=\"%s\" d=\"M%s,%s l%s,0\" />\n",
           $hpCols{$hpBase}, $xp, $tSeqCount*$sf+$sfy, length($hpSeq));
    $xp += length($hpSeq);
  }
  if($tSeqC){
    printf("   <path stroke=\"%s\" d=\"M%s,%s l%s,0\" />\n",
           contentColour($tSeqC), $xp, $tSeqCount*$sf+$sfy, length($tSeqC));
  }
  printf("  </g>\n");
  printf(" </g>\n");
}
