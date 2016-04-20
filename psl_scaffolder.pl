#!/usr/bin/perl

use warnings;
use strict;

use Pod::Usage; ## uses pod documentation in usage code
use Getopt::Long qw(:config auto_version auto_help pass_through);

our $VERSION = "1.00";

=head1 NAME

psl_scaffolder.pl - use self-mapped PSL file to scaffold a genome

=head1 SYNOPSIS

./psl_scaffolder.pl -query <file> [options] <mapping.psl>

=cut

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

printf(STDERR " loaded in %d sequences\n", scalar(keys(%querySeqs)));

print(STDERR "Processing results...");
while(<>){
  chomp;
  my ($matches, $misMatches, $repMatches, $nCount, $qNumInsert,
      $qBaseInsert, $tNumInsert, $tBaseInsert, $strand, $qName, $qSize,
      $qStart, $qEnd, $tName, $tSize, $tStart, $tEnd, $blockCount,
      $blockSizes, $qStarts, $tStarts, @rest) = split(/\t/);
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
  if($pid >= $projOpts->{"pid"}){
  }
  #printf("%0.1f\n", $pid);
}
print(STDERR " done\n");
