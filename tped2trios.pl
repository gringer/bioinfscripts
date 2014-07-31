#!/usr/bin/perl

# tped2trios.pl -- Convert plink tped/tfam format file to a BEAGLE
# unphased trio data file. If a second file is provided indicating the
# genotyped individuals, then only trios containing two (or three)
# genotyped individuals will be output.

# Author: David Eccles (gringer), 2012 <bioinformatics@gringer.org>

use strict;
use warnings;

sub usage {
  print("usage: ./tped2trios.pl <tped file> [id list file] [options]\n");
  print("\nConvert plink tped/tfam file to BEAGLE trio file\n");
  print("\nOther Options:\n");
  print("-3      : ensure all trio members are genotyped (instead of 2)\n");
  print("-help   : Only display this help message\n");
  print("\n");
}

my @files = ();

my $tpedFileName = "";
my $listFileName = "";
my $trioMin = 2;

# extract command line arguments
while(@ARGV){
    my $argument = shift @ARGV;
    if(-f $argument){ # file existence check
        if(!$tpedFileName){
            $tpedFileName = $argument;
            printf(STDERR "Setting TPED file name to '%s'\n", $tpedFileName);
        } elsif (!$listFileName) {
            $listFileName = $argument;
            printf(STDERR "Setting genotyped individual list file ".
                   "name to '%s'\n", $listFileName);
        } else {
            print(STDERR "Error: only two files can be specified\n");
            usage();
            exit(1);
        }
    } else {
        if($argument eq "-3"){
          $trioMin = 3;
            print(STDERR "All trio individuals output will be in list file\n");
        }
        if($argument eq "-help"){
            usage();
            exit(0);
        }
    }
}

if(!$tpedFileName){
    print(STDERR "Error: no files specified, cannot continue\n");
    usage();
    exit(1);
}

my $tfamFileName = $tpedFileName;
$tfamFileName =~ s/\.tped$/\.tfam/;

if(!(-f $tfamFileName)){
    print(STDERR "Error: TFAM file associated with TPED file does not exist\n");
    printf(STDERR " - expecting '%s'\n", $tfamFileName);
    usage();
    exit(1);
}

if(!$listFileName){
    $listFileName = $tfamFileName;
}

my %genotypedIndivs = ();
# get IDs for people who have been genotyped
my $lineNum = 0;
open(my $listFile, "<", $listFileName) or die("Cannot open $listFileName");
while(<$listFile>){
    chomp;
    my @F = split(/\s+/, $_);
    # assume field 0 is family, field 1 is individual
    $genotypedIndivs{$F[1]."@".$F[0]} = $lineNum++;
}
close($listFile);
printf(STDERR "Added %d genotyped IDs\n", scalar(keys(%genotypedIndivs)));

my %tfamIndivs = ();
my @trioIDs = ();
# get locations/IDs and parents for indviduals in TFAM file
$lineNum = 0;
open(my $tfamFile, "<", $tfamFileName) or die("Cannot open $tfamFileName");
while(<$tfamFile>){
    chomp;
    my @F = split(/\s+/, $_);
    # assume field 0 is family, field 1 is individual
    # beagle output should be parent, parent, offspring
    my @IDs = ($F[3]."@".$F[0], $F[2]."@".$F[0], $F[1]."@".$F[0]);
#    print(join(";",@IDs)."\n");
    $tfamIndivs{$IDs[2]} = $lineNum++;
    if(scalar(grep {$genotypedIndivs{$_}} @IDs) >= $trioMin){
        push(@trioIDs, @IDs);
    }
}
close($tfamFile);

foreach my $id (@trioIDs){
    if(!defined($tfamIndivs{$id})){
        printf("column ID for $id is undefined, bailing out\n");
        usage();
        exit(1);
        # $genotypedIndivs{$id} = -1;
    }
}

# convert trio IDs into column numbers
my @trioColumns = 
    map {($tfamIndivs{$_} * 2, $tfamIndivs{$_} * 2 + 1)} @trioIDs;

# print first line (list of IDs)
print("I ind\@fam ".join(" ", map {($_,$_)} @trioIDs)."\n");

open(my $tpedFile, "<", $tpedFileName) or die("Cannot open $tpedFileName");
print(STDERR "Writing output (one '.' per 1,000 lines)...");
my $lineCount = 0;
while(<$tpedFile>){
    chomp;
    my @F = split(/\s+/, $_);
    push(@F, 0, 0);
    shift(@F); # remove chromosome number
    my $marker = shift(@F); # extract marker name
    shift(@F);shift(@F); # remove map position / base location
    # print out marker line with trio genotypes
    printf("M $marker %s\n", join(" ", map {$F[$_]} @trioColumns));
    if($lineCount++ % 1000 == 0){
        print(STDERR ".");
    }
}
print(STDERR " done!\n");
