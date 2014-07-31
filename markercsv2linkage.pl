#!/usr/bin/perl

# csv2merlin.pl -- Copy genotype data and pedigree information
# to create MERLIN-compatible QTDT files

# Author: David Eccles (gringer), 2013 <bioinformatics@gringer.org>

use warnings;
use strict;

use Text::CSV;

sub usage {
  print("usage: ./csv2merlin.pl <PED file> <genotype file> [options]\n");
  print("\nConvert CSV format to MERLIN-QTDT\n");
  print("\nOther Options:\n");
  print("-help        : Only display this help message\n");
  print("\n");
}

my @files = ();

my $pedFileName = "";
my $gtFileName = "";

# extract command line arguments
while(@ARGV){
    my $argument = shift @ARGV;
    if(-f $argument){ # file existence check
        if(!$pedFileName){
            $pedFileName = $argument;
            #printf(STDERR "Setting PED file name to '%s'\n", $pedFileName);
        } elsif (!$gtFileName) {
            $gtFileName = $argument;
            #printf(STDERR "Setting genotype file ".
            #       "name to '%s'\n", $gtFileName);
        } else {
            print(STDERR "Error: only two files can be specified\n");
            usage();
            exit(1);
        }
    } else {
        if($argument eq "-help"){
            usage();
            exit(0);
        }
    }
}

if(!$pedFileName || !$gtFileName){
    print(STDERR "Error: *two* files must be specified, cannot continue\n");
    usage();
    exit(1);
}

my $csv = Text::CSV->new ( { binary => 1 } )  # should set binary attribute.
   or die "Cannot use CSV: ".Text::CSV->error_diag ();

my $maxGTCount = 0;
my %gtLines = ();
open(my $gtFile, "<", $gtFileName) or die("Cannot open genotype file");
while(my $row = $csv->getline($gtFile)){
  my @F = @{$row};
  my $id = shift(@F);
  if(scalar(@F) > $maxGTCount){
    $maxGTCount = scalar(@F);
  }
  map{$_ = ($_)?$_:"0/0"} @F;
  $gtLines{$id} = join(" ",@F)."\n";
}
close($gtFile);

my $blankLine = " 0/0" x $maxGTCount . "\n";

open(my $pedFile, "<", $pedFileName) or die("Cannot open PED file");
while(<$pedFile>){
  chomp;
  print $_." ";
  my @F = split(/\s+/, $_, 3);
  if($gtLines{$F[1]}){
    print($gtLines{$F[1]});
  } else {
    print($blankLine);
  }
}
