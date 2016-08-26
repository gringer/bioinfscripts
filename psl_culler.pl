#!/usr/bin/perl

use warnings;
use strict;
use English;

use Pod::Usage; ## uses pod documentation in usage code
use Getopt::Long qw(:config auto_version auto_help pass_through);

our $VERSION = "1.00";

=head1 NAME

psl_culler.pl - remove contigs that are substantially covered by other contigs

=head1 SYNOPSIS

./psl_culler.pl -query <file> [options] <self_mapping.psl>

=cut

sub min {
  ($a, $b) = @_;
  return( ($a < $b) ? $a : $b);
}

sub max {
  ($a, $b) = @_;
  return( ($a > $b) ? $a : $b);
}

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
   "prefix" => "psl_scaffold_", # prefix for contig names
   "pid" => 90, # percent ID threshold
   "capfrac" => 0.9, # fraction captured to exclude sequence
  };

GetOptions($projOpts, 'query=s', 'pid=i', 'capfrac=f', 'prefix=s');

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
    $querySeqs{$seqID}{length} = 0;
    $querySeqs{$seqID}{captured} = 0;
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
    $querySeqs{$seqID}{"length"} += length($_);
  }
}
close($queryFile);

my %targetSeqs = %querySeqs;
my $nextScaffoldID = 1;

my %replacementSeqs = ();

printf(STDERR " loaded in %d sequences (last ID %s)\n",
       scalar(keys(%querySeqs)), $seqID);

print(STDERR "Processing results...");
while(<>){
  chomp;
  if(!$_){
    next;
  }
  my @fields = split(/\t/);
  my ($matches, $misMatches, $repMatches, $nCount, $qNumInsert,
      $qBaseInsert, $tNumInsert, $tBaseInsert, $strand, $qName,
      $qSize, $qStart, $qEnd, $tName, $tSize,
      $tStart, $tEnd, $blockCount, $blockSizes, $qStarts,
      $tStarts, @rest) = @fields;
  if(!$tStarts){
    pod2usage({-exitVal => 1,
               -message => sprintf(" Error: mapping file doesn't look ".
                                   "like a PSL file (expecting".
                                   ">=21 tab-separated values, got %d)",
                                  scalar(@fields)),
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
    if($querySeqs{$qName}{length} > $targetSeqs{$tName}{length}){
      ## only remove larger sequences; this prevents a larger sequence
      ## from being removed, and also removing its sub-sequences
      $querySeqs{$qName}{captured} += $qEnd - $qStart;
    }
  } elsif($pid < $projOpts->{"pid"}){
    # printf(STDERR "Rejecting match '%s' vs '%s': identity (%f) too low\n",
    #      $qName, $tName, $pid);
  }
}
printf(STDERR " done\n");

my %displayed = ();

foreach my $seqID (sort(keys(%querySeqs))){
  if(!$seqID){
    next;
  }
  if($querySeqs{$seqID}{length} == 0){
    printf(STDERR "Length of %s is 0\n", $seqID);
  }
  my $fullName = $querySeqs{$seqID}{fullName};
  my $sequence = $querySeqs{$seqID}{sequence};
  my $capProp = $querySeqs{$seqID}{captured} ?
    ($querySeqs{$seqID}{captured} / $querySeqs{$seqID}{length}) : 0;
  if($capProp && ($capProp > $projOpts->{capfrac})){
    printf(STDERR "culling $seqID; %0.1f%% captured\n", $capProp * 100);
  } else {
    printf(">%s\n%s\n", $fullName, $sequence);
    $displayed{$fullName} = 1;
  }
}
