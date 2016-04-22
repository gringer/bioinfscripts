#!/usr/bin/perl

use warnings;
use strict;
use English;

use Pod::Usage; ## uses pod documentation in usage code
use Getopt::Long qw(:config auto_version auto_help pass_through);
use List::Util qw(max min); ## for max/min

our $VERSION = "1.00";

=head1 NAME

psl_scaffolder.pl - use self-mapped PSL file to scaffold a genome

=head1 SYNOPSIS

./psl_scaffolder.pl -query <file> [options] <mapping.psl>

=cut

sub rc {
  my ($seq) = @_;
  $seq =~ tr/ACGTUYRSWMKDVHBXN-/TGCAARYSWKMHBDVXN-/;
  # work on masked sequences as well
  $seq =~ tr/acgtuyrswmkdvhbxn/tgcaaryswkmhbdvxn/;
  return(scalar(reverse($seq)));
}

sub getConsensus {
  my ($b1, $b2) = @_;
  if(($b1 eq $b2) || ($b1 eq " ") || ($b2 eq " ")){
    ## equal bases, or absent bases, so consensus is easy
    return($b1);
  }
  # if different, convert to upper case to simplify lookup
  my $bc = uc(($b1 cmp $b2) ? $b1.$b2 : $b2.$b1);
  my %consensusLookup =
    (AC => "M", AM => "A", CM => "C",
     GT => "K", GK => "G", KT => "T",
     AG => "R", AR => "A", GR => "G",
     CT => "Y", CY => "C", TY => "T",
     AT => "W", AW => "A", TW => "T",
    );
  # if "simple" ambiguity can be found, return that, otherwise return N
  # (i.e. GT => K, -A => N, YM -> N)
  return( ($consensusLookup{$bc}) ? $consensusLookup{$bc} : "N");
}

sub getMatch {
  my ($b1, $b2) = @_;
  return((($b1 eq $b2) || ($b1 eq " ") || ($b2 eq " ") ||
         ($b1 eq "N") || ($b2 eq "N")) ? " " : "*");
}

############### Program starts here

# set default options
my @pslFiles = ();
my $projOpts =
  {
   "query" => 0, # contig file for query sequences
   "pid" => 90, # percent ID threshold
   "trimlimit" => 50, # max number of overlapping bases outside match region
  };

GetOptions($projOpts, 'query=s', 'pid=i', 'trimlimit=i');

# process remaining command line arguments (hopefully only PSL files)
while (@ARGV) {
  my $argument = shift @ARGV;
  if(-f $argument){
    push (@pslFiles, $argument);
  } else {
  pod2usage({-exitVal => 1,
               -message => "Error: Unknown command-line option or ".
             "non-existent file, '$argument'\n", -verbose => 0});
  }
}

@ARGV = @pslFiles;

if(!$projOpts->{"query"}){
  pod2usage({-exitVal => 1,
             -message => "Error: No query assembly file provided",
             -verbose => 0});
}

if(!(-f $projOpts->{"query"})){
  pod2usage({-exitVal => 1,
             -message => sprintf("Error: query file '%s' doesn't exist",
                                 $projOpts->{"query"}),
             -verbose => 0});
}

print(STDERR "Loading query sequences into memory...");
open(my $queryFile, "<", $projOpts->{"query"});
my $seqID = "";
my %querySeqs = ();
while(<$queryFile>){
  chomp;
  if(/^>((.+?)( .*?\s*)?)$/){
    ## line is sequence header
    $seqID = $2;
    $querySeqs{$seqID}{fullName} = $1;
    $querySeqs{$seqID}{sequence} = "";
  } else {
    if(!$seqID){
      pod2usage({-exitVal => 1,
                 -message => sprintf(" Error: query file '%s' doesn't look ".
                                     "like a FASTA file (no initial ID header)",
                                     $projOpts->{"query"}),
                 -verbose => 0});
    }
    ## line is sequence
    $querySeqs{$seqID}{"sequence"} .= $_;
  }
}
close($queryFile);

my %targetSeqs = %querySeqs;
my $nextScaffoldID = 1;

my %replacementSeqs = ();

printf(STDERR " loaded in %d sequences\n", scalar(keys(%querySeqs)));

print(STDERR "Processing results...");
while(<>){
  chomp;
  my @fields = split(/\t/);
  my ($matches, $misMatches, $repMatches, $nCount, $qNumInsert,
      $qBaseInsert, $tNumInsert, $tBaseInsert, $strand, $qName,
      $qSize, $qStart, $qEnd, $tName, $tSize,
      $tStart, $tEnd, $blockCount, $blockSizes, $qStarts,
      $tStarts, @rest) = @fields;
  if(!$tStarts){
    pod2usage({-exitVal => 1,
               -message => sprintf(" Error: mapping file doesn't look like a ".
                                   "PSL file ".
                                   "(expecting >=21 tab-separated values)",
                                   $projOpts->{"query"}),
               -verbose => 0});
  }
  ## calculate percent identity
  my $qAliSize = $qEnd - $qStart;
  my $tAliSize = $tEnd - $tStart;
  my $sizeDif = abs($qAliSize - $tAliSize);
  my $pid = 100 * ($matches + $repMatches -
                   ($qNumInsert + $tNumInsert + 3*log(1+$sizeDif))) /
                     ($matches + $repMatches + $misMatches);
  if(($pid >= $projOpts->{"pid"}) &&
     $querySeqs{$qName} && $targetSeqs{$tName}){
    my %meta = ();
    my $shortTarget = ($tSize < $qSize) ? 1 : 0;
    my $longTarget = (1 - $shortTarget);
    my $sName = $fields[9 + ($shortTarget * 4)];
    my $lName = $fields[9 + ($longTarget * 4)];
    my $sLen = $fields[10 + ($shortTarget * 4)];
    my $lLen = $fields[10 + ($longTarget * 4)];
    my $sStart = $fields[11 + ($shortTarget * 4)];
    my $lStart = $fields[11 + ($longTarget * 4)];
    my $sEnd = $fields[12 + ($shortTarget * 4)];
    my $lEnd = $fields[12 + ($longTarget * 4)];
    my @sBlStarts = split(/,/, $fields[19 + $shortTarget]);
    my @lBlStarts = split(/,/, $fields[19 + $longTarget]);
    my @blSizes = split(/,/, $fields[18]);
    my ($sSeq, $lSeq) = ($querySeqs{$qName}{sequence},
                         $querySeqs{$tName}{sequence});
    if($shortTarget){
      ($sSeq, $lSeq) = ($lSeq, $sSeq);
    }
    my $doRC = ($strand eq "-");
    my $preTrim = min($sStart, $lStart);
    my $postTrim = min($sLen - $sEnd, $lLen - $lEnd);
    ## Only continue on if there's a good likelihood that this will work
    ## i.e. trimLength * (1-%id) < threshold
    my $sPre = substr($sSeq, 0, $sStart);
    my $lPre = substr($lSeq, 0, $lStart);
    my $sMid = substr($sSeq, $sStart, $sEnd-$sStart);
    my $lMid = substr($lSeq, $lStart, $lEnd-$lStart);
    my $sPost = substr($sSeq, $sEnd);
    my $lPost = substr($lSeq, $lEnd);
    if ($doRC) {
      if ($shortTarget) { # target sequence is assumed to be forward strand
        ($lSeq, $lPre, $lMid, $lPost) = (rc($lSeq), rc($lPost), rc($lMid), rc($lPre));
        $preTrim = min($sStart, $lLen - $lEnd);
        $postTrim = min($sLen - $sEnd, $lStart);
      } else {
        ($sSeq, $sPre, $sMid, $sPost) = (rc($sSeq), rc($sPost), rc($sMid), rc($sPre));
        $preTrim = min($sLen - $sEnd, $lStart);
        $postTrim = min($sStart, $lLen - $lEnd);
      }
    }
    my $trimTotal = ($preTrim + $postTrim);
    if($trimTotal <= $projOpts->{"trimlimit"}){
      my $sPreTrim = substr($sPre, length($sPre)-$preTrim);
      my $sPostTrim = substr($sPost, 0, $postTrim);
      my $lPreTrim = substr($lPre, length($lPre)-$preTrim);
      my $lPostTrim = substr($lPost, 0, $postTrim);
      my $preLen = max(length($sPre), length($lPre));
      my $postLen = max(length($sPost), length($lPost));
      my $lastS = $sBlStarts[0];
      my $lastL = $lBlStarts[0];
      my $alSeqS = "";
      my $alSeqL = "";
      for (my $i = 0; $i <= $#blSizes; $i++) {
        my $gapS = $sBlStarts[$i] - $lastS;
        my $gapL = $lBlStarts[$i] - $lastL;
        my $gapLength = max($gapS, $gapL);
        my $fillS = $gapLength - $gapS;
        my $fillL = $gapLength - $gapL;
        $alSeqS .= ("-" x $fillS) . substr($sSeq, $sBlStarts[$i]-$gapS, $gapS);
        $alSeqL .= ("-" x $fillL) . substr($lSeq, $lBlStarts[$i]-$gapL, $gapL);
        $alSeqS .= substr($sSeq, $sBlStarts[$i], $blSizes[$i]);
        $alSeqL .= substr($lSeq, $lBlStarts[$i], $blSizes[$i]);
        $lastS = $sBlStarts[$i] + $blSizes[$i];
        $lastL = $lBlStarts[$i] + $blSizes[$i];
      }
      $alSeqS = $sPreTrim . $alSeqS . $sPostTrim;
      $alSeqL = $lPreTrim . $alSeqL . $lPostTrim;
      my $alConsensus = "";
      for (my $i = 0; $i < length($alSeqS); $i++) {
        $alConsensus .= getConsensus(substr($alSeqS,$i,1),substr($alSeqL,$i,1));
      }
      my $consensusLength = length($alConsensus);
      $alConsensus =
        substr($sPre, 0, length($sPre) - $preTrim).
          substr($lPre, 0, length($lPre) - $preTrim).
            $alConsensus.
              substr($sPost, $postTrim).substr($lPost, $postTrim);
      my $newSeqID = sprintf("psl_scaffold_%d", $nextScaffoldID++);
      if(!exists($replacementSeqs{$sName}{score}) ||
         ($trimTotal < $replacementSeqs{$sName}{score}) ||
         (($trimTotal == $replacementSeqs{$sName}{score}) &&
          ($consensusLength > $replacementSeqs{$sName}{clength}))){
        $replacementSeqs{$sName}{score} = $trimTotal;
        $replacementSeqs{$sName}{clength} = $consensusLength;
        $replacementSeqs{$sName}{fullName} =
          sprintf("%s [%s %s]", $newSeqID, $sName, $lName);
        $replacementSeqs{$sName}{sequence} = $alConsensus;
      }
      if(!exists($replacementSeqs{$lName}{score}) ||
         ($trimTotal < $replacementSeqs{$lName}{score}) ||
         (($trimTotal == $replacementSeqs{$lName}{score}) &&
          ($consensusLength > $replacementSeqs{$lName}{clength}))){
        $replacementSeqs{$lName}{score} = $trimTotal;
        $replacementSeqs{$lName}{clength} = $consensusLength;
        $replacementSeqs{$lName}{fullName} =
          sprintf("%s [%s %s]", $newSeqID, $sName, $lName);
        $replacementSeqs{$lName}{sequence} = $alConsensus;
      }
    }
  }
}
printf(STDERR " done\n");

my %displayed = ();

foreach my $seqID (sort(keys(%targetSeqs))){
  my $fullName = $targetSeqs{$seqID}{fullName};
  my $sequence = $targetSeqs{$seqID}{sequence};
  if(exists($replacementSeqs{$seqID})){
    print(STDERR "Found match for $seqID\n");
    $fullName = $replacementSeqs{$seqID}{fullName};
    $sequence = $replacementSeqs{$seqID}{sequence};
  }
  if(!$displayed{$fullName}){
    printf(">%s\n%s\n", $fullName, $sequence);
    $displayed{$fullName} = 1;
  }
}

foreach my $seqID (sort(keys(%replacementSeqs))){
  if(!$displayed{$seqID}){
    print(STDERR "No match for $seqID\n");
  }
}
