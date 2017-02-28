#!/usr/bin/perl

## Splits up a FASTA/FASTQ file based on the location of identified sequences (e.g. adapter/barcode sequences)
## [currently only identifies the location of the sequences]

use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $seqFileName = "";

GetOptions("seqfile=s" => \$seqFileName) or
  die("Error in command line arguments");

my %seqs = ();
my %quals = ();

# if(!$seqFileName){
#   die("Error: read sequence file must be specified, '-seqfile <file.fa>'");
# }

# ## read in sequences
# my $inQual = 0; # false
# my $seqID = "";
# my $qualID = "";

# my $seqFile = new IO::Uncompress::Gunzip "$seqFileName" or
#   die "Unable to open $seqFileName\n";
# while(<$seqFile>){
#   chomp;
#   chomp;
#   if(!$inQual){
#     if(/^(>|@)((.+?)( .*?\s*)?)$/){
#       my $newSeqID = $2;
#       my $newShortID = $3;
#       $seqID = $newShortID;
#       if($seqID){
#         $seqs{$seqID} = "";
#         $quals{$seqID} = "";
#       }
#     } elsif(/^\+(.*)$/) {
#       $inQual = 1; # true
#       $qualID = $1;
#       if(($qualID ne "") && ($qualID ne $seqID)){
#         warn("Sequence ID and Qual ID do not match: $seqID; $qualID");
#       }
#     } else {
#       $seqs{$seqID} .= $_;
#     }
#   } else {
#     $quals{$seqID} .= $_;
#     my $lq = length($quals{$seqID});
#     my $ls = length($seqs{$seqID});
#     if($lq >= $ls){
#       $inQual = 0; # false
#       if($lq != $ls){
#         warn(sprintf("Sequence and Qual length do not match: $seqID (%d; %d)",
#                     $ls, $lq));
#       }
#     }
#   }
# }
# close($seqFile);

printf(STDERR "Read in %d sequences\n", scalar(keys(%seqs)));

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

my %matches = ();

printf("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n",
       "query", "target", "dir", "qS", "qE", "qML", "qL",
       "tS", "tE", "tML", "tL");
while(<>){
  if(!/^[as]/){
    next;
  }
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
      if($qStrand eq "-"){ ## correct for reverse complement
        $qEnd = $qLen - $qStart;
        $qStart = $qEnd - $qMatchLen;
        $tEnd = $tLen - $tStart;
        $tStart = $tEnd - $tMatchLen;
      }
      my $matchLine =
        sprintf("%s,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d",
                $qName, $tName, $qStrand, $qStart, $qEnd, $qMatchLen, $qLen,
                $tStart, $tEnd, $tMatchLen, $tLen);
      print("$matchLine\n");
      $matches{$qName}{$qStart} .= ":matchLine";
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

