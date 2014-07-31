#!/usr/bin/perl

# simplifyPed.pl -- Alter family IDs in a plink PED format file to
# simplify the pedigree structure. If a second file is provided
# indicating the genotyped individuals, then only trios with at least
# two genotyped individuals will be output and considered part of the
# same family.

# Author: David Eccles (gringer), 2012 <bioinformatics@gringer.org>

use strict;
use warnings;

sub usage {
  print("usage: ./simplifyPed.pl <PED file> [id list file] [options]\n");
  print("\nSimplify pedigree structure in plink PED file\n");
  print("\nOther Options:\n");
  print("-help        : Only display this help message\n");
  print("-phenogeno   : Set phenotype column to indicate genotyped status\n");
  print("\n");
}

my @files = ();

my $pedFileName = "";
my $listFileName = "";
my $phenoGenotype = 0; # false

# extract command line arguments
while(@ARGV){
    my $argument = shift @ARGV;
    if(-f $argument){ # file existence check
        if(!$pedFileName){
            $pedFileName = $argument;
            printf(STDERR "Setting PED file name to '%s'\n", $pedFileName);
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
        if($argument eq "-help"){
            usage();
            exit(0);
        } elsif($argument eq "-phenogeno"){
            $phenoGenotype = 1; # true
        }
    }
}

if(!$pedFileName){
    print(STDERR "Error: no files specified, cannot continue\n");
    usage();
    exit(1);
}

if(!$listFileName){
    $listFileName = $pedFileName;
    print(STDERR "Setting list file name to '$listFileName'\n");
}

my %genotypedIndivs = ();
# get IDs for people who have been genotyped
my $lineNum = 1;
open(my $listFile, "<", $listFileName) or die("Cannot open $listFileName");
while(<$listFile>){
    chomp;
    my @F = split(/\s+/, $_, 10);
    # assume field 0 is family, field 1 is individual
    $genotypedIndivs{$F[1]."@".$F[0]} = $lineNum++;
}
close($listFile);

my $maxFamID = 0;
my %indivFams = ();
my %pedIndivs = ();
my %trioIncluded = ();
# get IDs and parents for indviduals in PED file
$lineNum = 1;
open(my $pedFile, "<", $pedFileName) or die("Cannot open $pedFileName");
print(STDERR "Reading PED file...");
while(<$pedFile>){
    chomp;
    my @F = split(/\s+/, $_, 10);
    # assume field 0 is family, field 1 is individual
    if($maxFamID < $F[0]){
        $maxFamID = $F[0];
    }
    my @IDs = ($F[1]."@".$F[0], $F[2]."@".$F[0], $F[3]."@".$F[0]);
    print(STDERR ".");
    my $family = $lineNum++;
    if(scalar(grep {defined($genotypedIndivs{$_})} @IDs) >= 2){
        #printf(STDERR "trio included: %s [%d genotyped]\n",
        #    join(";",@IDs), scalar(grep {defined($genotypedIndivs{$_})} @IDs));
        $trioIncluded{$IDs[0]} = 1;
        for(my $i = 0; $i < 3; $i++){
            if(defined($indivFams{$IDs[$i]})){
                $family = $indivFams{$IDs[$i]};
                #printf(STDERR "individual %s has been seen before in family %s\n", $IDs[$i], $family);
            }
        }
        # sort out disagreements
        for(my $i = 0; $i < 3; $i++){
            if(!defined($indivFams{$IDs[$i]})){
                if($IDs[$i] !~ /^0@/){
                    $indivFams{$IDs[$i]} = $family;
                }
            } elsif($indivFams{$IDs[$i]} != $family) {
                my $disagreeFamily = $indivFams{$IDs[$i]};
                #printf(STDERR "individual %s has a disagreement with family %s\n", $IDs[$i], $family);
                my $minFamily = $family;
                my $otherFamily = $disagreeFamily;
                if($disagreeFamily < $minFamily){
                    $minFamily = $disagreeFamily;
                    $otherFamily = $family;
                }
                $family = $minFamily;
                for my $id (keys(%indivFams)){
                    if($indivFams{$id} == $otherFamily){
                        $indivFams{$id} = $minFamily;
                    }
                }
            }
        }
        # patch up remaining conflicts
        for(my $i = 0; $i < 3; $i++){
            if($IDs[$i] !~ /^0@/){
                $indivFams{$IDs[$i]} = $family;
            }
        }
    }
}
print(STDERR "\n");
close($pedFile);

my $nextUnknownFam = $lineNum + $maxFamID + 1;

open($pedFile, "<", $pedFileName) or die("Cannot open $pedFileName");
$lineNum = 1;
while(<$pedFile>){
    chomp;
    my @F = split(/\s+/, $_, 10);
    my $family = $lineNum++;
    my $indCode = $F[1]."@".$F[0];
    if(defined($indivFams{$indCode})){
        if(!defined($trioIncluded{$indCode})){
            # a genotyped parent who should be considered a founder
            $F[2] = 0;
            $F[3] = 0;
        }
        #printf(STDERR "found individual on line %d (family %d)\n", $lineNum - 1, $indivFams{$indCode});
        $F[0] = $indivFams{$indCode} + $maxFamID;
        if($phenoGenotype){ # phenotypes get repurposed as indicating genotyped status
            $F[5] = defined($genotypedIndivs{$indCode}) ? 2 : 1;
        }
    } else {
       # remaining unlinked individuals get their own family. This reduces tree complexity
      $F[0] = $nextUnknownFam++;
    }
    print(join(" ", @F)."\n");
}
close($pedFile);
