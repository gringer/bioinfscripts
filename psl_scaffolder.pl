#!/usr/bin/perl

use warnings;
use strict;

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
  $seq =~ tr/ACGTUYRSWMKDVHBXN-/TGCARYSWKMHBDVXN-/;
  # work on masked sequences as well
  $seq =~ tr/acgtuyrswmkdvhbxn/tgcaryswkmhbdvxn/;
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
    (AC => "M",
     GT => "K",
     AG => "R",
     CT => "Y",
     AT => "W");
  # if "simple" ambiguity can be found, return that, otherwise return N
  # (i.e. GT => K, -A => N, YM -> N)
  return( ($consensusLookup{$bc}) ? $consensusLookup{$bc} : "N");
}

sub getMatch {
  my ($b1, $b2) = @_;
  return(($b1 eq $b2) ? "*" : " ");
}


############### Program starts here

# set default options
my @pslFiles = ();
my $projOpts =
  {
   "query" => 0, # contig file for query sequences
   "pid" => 90, # percent ID threshold
  };

GetOptions($projOpts, 'query=s', 'pid=i');

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
    $querySeqs{$seqID}{"fullName"} = $1;
    $querySeqs{$seqID}{"sequence"} = "";
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
    my $shortTarget = ($qAliSize >= $tAliSize) ? 1 : 0;
    my $longTarget = ($qAliSize < $tAliSize) ? 1 : 0;
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
    my $sPre = substr($sSeq, 0, $sStart);
    my $lPre = substr($lSeq, 0, $lStart);
    my $sMid = substr($sSeq, $sStart, $sEnd-$sStart);
    my $lMid = substr($lSeq, $lStart, $lEnd-$lStart);
    my $sPost = substr($sSeq, $sEnd);
    my $lPost = substr($lSeq, $lEnd);
    my $doRC = ($strand eq "-");
    if($doRC){
      if($shortTarget){  # target sequence is assumed to be forward strand
        ($lPre, $lMid, $lPost) = (rc($lPost), rc($lMid), rc($lPre));
        $lSeq = rc($lSeq);
      } else {
        ($sPre, $sMid, $sPost) = (rc($sPost), rc($sMid), rc($sPre));
        $sSeq = rc($sSeq);
      }
    }
    my $preLength = max(length($sPre), length($lPre));
    printf("%s %s\n", $sName, $lName);
    printf("%d/%d-%d %d/%d-%d\n", $sLen, $sStart, $sEnd, $lLen, $lStart, $lEnd);
    printf("%d %d\n", length($sSeq), length($lSeq));
    #printf("%s [%d bp]\n",
    #      $sSeq, length($sSeq));
    #printf("%${preLength}s %s %s [%d bp]\n",
    #      $sPre, $sMid, $sPost, length($sPre.$sMid.$sPost));
    printf("%s [%d bp]\n",
           $sMid, length($sMid));
    #printf("%s [%d bp]\n",
    #       $lSeq, length($lSeq));
    #printf("%${preLength}s %s %s [%d bp]\n",
    #       $lPre, $lMid, $lPost, length($lPre.$lMid.$lPost));
    printf("%s [%d bp]\n",
           $lMid, length($lMid));
    print("Blocks: \n");
    print(join(" ",@blSizes)."\n");
    print(join(" ",@sBlStarts)."\n");
    print(join(" ",@lBlStarts)."\n");
    my $lastS = $sBlStarts[0];
    my $lastL = $lBlStarts[0];
    my $alSeqS = "";
    my $alSeqL = "";
    for(my $i = 0; $i <= $#blSizes; $i++){
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
    my $alConsensus = "";
    my $alMatch = "";
    for(my $i = 0; $i < length($alSeqS); $i++){
      $alConsensus .= getConsensus(substr($alSeqS,$i,1),substr($alSeqL,$i,1));
      $alMatch .= getMatch(substr($alSeqS,$i,1),substr($alSeqL,$i,1));
    }
    print($alSeqS."\n");
    print($alConsensus."\n");
    print($alMatch."\n");
    print($alSeqL."\n");
  }
  #printf("%0.1f\n", $pid);
}
print(STDERR " done\n");
