#!/usr/bin/perl

# addKraken2Silva.pl -- modifies SILVA FASTA files to add in
# embl taxIDs for the Kraken database and translate U to T

# Author: David Eccles (gringer), 2015 <bioinformatics@gringene.org>

use strict;
use warnings;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

sub usage {
  print(STDERR "usage: ./addKraken2Silva.pl <Map File> <input file>\n");
  print(STDERR "\nmodifies SILVA FASTA files to add in embl taxIDs for the Kraken database\n");
  print(STDERR "\nOther Options:\n");
  print(STDERR "-help        : Only display this help message\n");
  print(STDERR "\n");
}

my $mapFileName = 0; # false

my @files = ();

# extract command line arguments
while(@ARGV){
    my $argument = shift @ARGV;
    if(-f $argument){ # file existence check
        if(!$mapFileName){
            $mapFileName = $argument;
        } else {
	    push(@files, $argument);
        }
    } else {
        if($argument eq "-help"){
            usage();
            exit(0);
        } else {
	    print(STDERR "Error: command line parameter '$argument' not understood\n");
	    usage();
	    exit(1);
	}
    }
}

@ARGV = @files;

if(!$mapFileName){
    print(STDERR "Error: No valid map file given\n");
    usage();
    exit(1);
}

my $mapFile = 0;
$mapFile = new IO::Uncompress::Gunzip "$mapFileName" or
    die "Unable to open $mapFileName\n";

my %seqMap = ();

print(STDERR "Reading from map file...");

my $mapLinesCounter = 0;

while(<$mapFile>){
    chomp;
    my $seqID = 0;
    my $taxID = 0;
    if($_ =~ /^(.*?)\s/){
	$seqID = $1;
    }
    if($_ =~ /\s([0-9]+)$/){
	$taxID = $1;
    }
    if($seqID && $taxID){
	$seqMap{$seqID} = $taxID;
    }
    if($mapLinesCounter++ > 100000){
	$mapLinesCounter = 0;
	print(STDERR ".");
	last;
    }
}

printf(STDERR " found %d mappings from sequence IDs to taxa\n", scalar(keys(%seqMap)));

my $seqID = "";
my $seq = "";
while(<>){
  chomp;
  if(/^>(([^\.]+).*?)(( |$).*)$/){
      my $seqBase = $2;
      my $seqID = $1;
      my $rest = $3;
      ##print(STDERR "base: $seqBase, id:$seqID, rest:$rest\n");
      my $taxID = $seqMap{$seqBase};
      if($taxID){
	  $seqID.= "|kraken:taxid|".$taxID;
      }
      printf(">%s%s\n", $seqID, $rest);
      $seq = "";
  } else {
      $_ =~ tr/U/T/;
      print($_."\n");
  }
}
