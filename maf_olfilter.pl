#!/usr/bin/perl

## Filters out non-end and identical LAST results to simulate overlap alignments

use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $seqFileName = "";
my $ignoreSelf = 0;

my $endLeniency = 50; ## base pairs from end to allow
my $containFrac = 0.9; ## minimum containment fraction
my $containIdent = 0.9; ## minimum containment identity

my %quals = ();

my $qSeq = "";
my $qStart = 0;
my $qEnd = 0;
my $qLen = 0;
my $qStrand = "";
my $qMatchLen = 0;
my $qName = "";
my $tSeq = "";
my $tStart = 0;
my $tEnd = 0;
my $tMatchLen = 0;
my $tLen = 0;
my $tName = "";

my $lineBuffer = "";
my $GFAspec = 1;

my %matches = ();

if($GFAspec == 2){
  print("H\tVN:Z:2.0\n");
} else {
  print("H\tVN:Z:bogart/edges\n");
}

my $fastaFileName = shift(@ARGV);
if(!(-e "${fastaFileName}.fai")){
  print(STDERR "Generating index file... ");
  system(("samtools", "faidx", $fastaFileName));
  print("done\n");
}

open(my $fastaFile, "<", "${fastaFileName}.fai");
while(<$fastaFile>){
  chomp;
  my @F = split(/\s+/);
  if($GFAspec == 2){
    printf("S\t%s\t%d\t*\n", $F[0], $F[1]);
  } else {
    printf("S\t%s\t*\t%d\n", $F[0], $F[1]);
  }
}
close($fastaFile);

while(<>){
  if(/^$/){
    next;
  }
  if(!/^[as]/){
    next;
  }
  $lineBuffer .= $_;
  my @F = split(/\s+/);
  if($F[0] eq "a"){
    $qSeq = "";
    $tSeq = "";
  } elsif($F[0] eq "s"){
    if($tSeq){
      $qName = $F[1];
      $qStart = $F[2];
      $qMatchLen = $F[3];
      $qEnd = $qStart + $qMatchLen;
      $qStrand = $F[4];
      $qLen = $F[5];
      $qSeq = $F[6];
      if($GFAspec == 2){
        if(($qSeq ne $tSeq) &&
           (($qStart < $endLeniency) || ($tStart < $endLeniency) ||
            ($qEnd > ($qLen-$endLeniency)) || ($tEnd > ($tLen-$endLeniency)))){
          print(join("\t",("E","${tName}+","${qName}${qStrand}",$tStart,$tEnd,
                           $qStart,$qEnd))."\n");
        }
      } else {
        if($qSeq ne $tSeq){
          if(($tEnd > ($tLen-$endLeniency)) && ($qStart < $endLeniency)){
            print(join("\t",("L",$tName,"+",$qName,"+","*"))."\n");
          }
          if(($qEnd > ($qLen-$endLeniency)) && ($tStart < $endLeniency)){
            print(join("\t",("L",$qName,"+",$tName,"+","*"))."\n");
          }
          if(($tStart < $endLeniency) && ($qStart < $endLeniency)){
            print(join("\t",("L",$tName,"-",$qName,"+","*"))."\n");
          }
          if(($tEnd > ($tLen-$endLeniency)) && ($qEnd > ($qLen-$endLeniency))){
            print(join("\t",("L",$tName,"+",$qName,"-","*"))."\n");
          }
        }
      }
      $lineBuffer = "";
   } else {
      $tName = $F[1];
      $tStart = $F[2];
      $tMatchLen = $F[3];
      $tEnd = $tStart + $tMatchLen;
      $tLen = $F[5];
      $tSeq = $F[6];
    }
  }
}

