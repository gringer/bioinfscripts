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
# These colours can be presented as-is, or increased to full
# saturation. With a consistent saturation setting, the same
# subsequence should appear identical regardless of its location
# (except possibly at the start and end of the sequence).
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

use Getopt::Long qw(:config auto_version auto_help pass_through);
use List::Util qw(min);

my $hpLength = 1;
my $saturate = 0;

GetOptions('hplength=i' => \$hpLength, 'saturate!' => \$saturate) or
  die("Error in command line arguments");

my $seq = "";
my $seqCount = 1;
my $seqID = "";

print("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
print("<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" ".
      "width=\"100%\" height=\"100%\" preserveAspectRatio=\"xMidYMid\">\n");


while(<>){
  chomp;
  if(/^>(.+)$/){
    my $newID = $1;
    if($seq){
      drawSeq($seq, $seqID, $seqCount, $hpLength, $saturate);
      $seqCount++;
    }
    $seq = "";
    $seqID = $newID;
  } else {
    $seq .= $_;
  }
}

if($seq){
  drawSeq($seq, $seqID, $seqCount, $hpLength, $saturate);
}

print("</svg>\n");


sub comp{
  my ($tSeq) = @_;
  $tSeq =~ tr/ACGT/TGCA/;
  return($tSeq);
}

sub contentColour{
  my ($tSeq, $doSat) = @_;
  my $total = ($tSeq =~ tr/AaCcGgTt//);
  if($total == 0){
    return(("#000000", 0.5, 0.5, 0.5));
  }
  my $a = ($tSeq =~ tr/Aa//);
  my $c = ($tSeq =~ tr/Cc//);
  my $g = ($tSeq =~ tr/Gg//);
  my $t = ($tSeq =~ tr/Tt//);
  my ($M, $Y, $S) = map {$_ / $total} (($c+$a), ($c+$t), ($c+$g));
  if($doSat){
    # fully saturate colours using HSP model, increasing R/G/B until
    # one colour reaches 1
    # code derived from Darel Rex Finley's function
    # see: http://alienryderflex.com/saturation.html
    my ($pR, $pG, $pB) = (0.299, 0.587, 0.114);
    my $p = sqrt($M*$M*$pR + $Y*$Y*$pG + $S*$S*$pB);
    my ($dR, $dG, $dB) = ($M-$p, $Y-$p, $S-$p);
    my ($kR, $kG, $kB) = map {$_ ? (1-$p) / $_ : 1} ($dR, $dG, $dB);
    my $minK = min(grep {$_ > 0} ($kR, $kG, $kB));
    ($M, $Y, $S) = map {$p + $_ * $minK} ($dR, $dG, $dB);
    # clip at 0 and 1
    ($M, $Y, $S) = map {($_<0) ? 0 : (($_>1) ? 1 : $_)} ($M, $Y, $S);
  }
  my ($Mc, $Yc, $Sc) = map {$_ * 255} ($M, $Y, $S);
  return(sprintf("#%02X%02X%02X", $Mc, $Yc, $Sc));
}

sub drawSeq{
  my ($tSeq, $tSeqID, $tSeqCount, $hl, $doSat) = @_;
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
      my $col = contentColour($preSeq, $doSat);
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
    my $col = contentColour($tSeq, $doSat);
    ## display sequence after last homopolymer
    printf("   <path stroke=\"%s\" d=\"M%s,%s l%s,0\" />\n",
           $col, $xp, $tSeqCount*$sf-$ofy, length($tSeq));
  }
  printf("  </g>\n");
  $xp = $sf;
  ## complement sequence
  printf("  <g id=\"%s_COMP\">\n", $tSeqID);
  while($tSeqC =~ s/^(.*?)(A{$hl,}|C{$hl,}|G{$hl,}|T{$hl,})//){
    my $preSeq = $1;
    my $hpSeq = $2;
    my $hpBase = substr($hpSeq,0,1);
    if($preSeq){
      ## display sequence prior to homopolymer
      my $col = contentColour($preSeq, $doSat);
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
    my $col = contentColour($tSeqC, $doSat);
    ## display sequence after last homopolymer
    printf("   <path stroke=\"%s\" d=\"M%s,%s l%s,0\" />\n",
           $col, $xp, $tSeqCount*$sf+$ofy, length($tSeqC));
  }
  printf("  </g>\n");
  printf(" </g>\n");
}
