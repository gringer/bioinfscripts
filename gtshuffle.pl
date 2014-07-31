#!/usr/bin/perl

# gtshuffle.pl -- Shuffles a simplegt input file, producing output
# files containing random splits of the popultion

# Author: David Eccles (gringer), 2008 <programming@gringer.org>

# derived from gt2plink.pl

use strict;
use warnings;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use FileHandle;

sub usage {
  print("usage: ./gtshuffle.pl <input file> [options]\n");
  print("\nRandomly splits a simple genotype file into groups of individuals\n");
  print("\nOther Options:\n");
  print("-help        : Only display this help message\n");
  print("-output      : output base file name\n");
  print("-r <float>   : Proportion of individuals in the first group (Default: 0.5)\n");
  print("-n <integer> : Number of individuals in the first group (overrides ratio)\n");
  print("-s <integer> : Number of groups to split file into (overrides number, ratio)");
  print("\n");
}

# array shuffling function, thanks to BrowserUk on perlmonks.org
# taken from http://www.perlmonks.org/?node_id=626155
# modified slightly to turn it into something that works with 'use strict'
sub shuffle_bukNew  {
    # store a reference to the input array
    my $ref = \@_;
    # generate a list of numbers from 0 .. <length of input array - 1>
    my @x = 0 .. (scalar(@{$ref}) - 1);
    # create a new array generated from pieces spliced out of random
    # locations of the original array
    @{ $ref }[ map splice( @x, rand @x, 1 ), @x ];
}

my %markers = ();
my $complement = 0; # false # determines if complement, rather
                            # than intersect, should be chosen

my $genotypesInFilename = 0; # false
my $mapInFilename = 0; # false
my $ratioFirst = 0.5; # equal proportions in each group
my $countFirst = ""; # empty if not specified

my $numSplits = 2;

my $outputBaseName = 0; # false

# extract command line arguments
while(@ARGV){
    my $argument = shift @ARGV;
    if(-f $argument){ # file existence check
        if($genotypesInFilename){
            $mapInFilename = $argument;
        } else {
            $genotypesInFilename = $argument;
        }
    } else {
        if($argument eq "-help"){
            usage();
            exit(0);
        }
        if($argument eq "-output"){
            $outputBaseName = shift @ARGV;
        }
        if($argument eq "-r"){
            $ratioFirst = shift @ARGV;
        }
        if($argument eq "-s"){
            $numSplits = shift @ARGV;
        }
        if($argument eq "-n"){
            $countFirst = shift @ARGV;
        }
    }
}

if(!$genotypesInFilename){
    print(STDERR "Error: No valid genotype input file given\n");
    usage();
    exit(1);
}

if($numSplits < 1){
    print(STDERR "Error: Too few splits (please split into at least 1 group)\n");
    usage();
    exit(2);
}

# if no base name is provided, derive one from input file name
if(!$outputBaseName){
    if($genotypesInFilename =~ m/^(.*?)\./){
        # truncate at location of first '.'
        $outputBaseName = $1;
    } else {
        # just in case the genotype input file doesn't contain '.'
        $outputBaseName = $genotypesInFilename;
    }
}

my $foundFiles = "";
for (my $i = 1; $i <= $numSplits; $i++){
    if(-f $outputBaseName."_".$i.".txt"){
        $foundFiles = $foundFiles . ($outputBaseName."_".$i.".txt\n");
    }
}

if($foundFiles){
    print(STDERR "Error: Split output files already exist\n".
          "Please delete these before running this program:\n".
          $foundFiles."\n");
    usage();
    exit(2);
}


my $genotypesInFile = 0;

$genotypesInFile = new IO::Uncompress::Gunzip "$genotypesInFilename" or
    die "Unable to open $genotypesInFilename\n";

my $gtDataLength = 0;

my @outFiles = ();
my $filename;
my $fh;

print(STDERR "Preparing output files: \n");
for (my $i = 1; $i < $numSplits; $i++){
    printf(STDERR "%s_%d.txt\n", $outputBaseName, $i);
    $filename = sprintf("%s_%d.txt", $outputBaseName, $i);
    $fh = new FileHandle;
    $fh->open("> $filename");
    push(@outFiles, $fh);
}
printf(STDERR "%s_%d.txt\n", $outputBaseName, $numSplits);
$filename = sprintf("%s_%d.txt", $outputBaseName, $numSplits);
$fh = new FileHandle;
$fh->open("> $filename");
push(@outFiles, $fh);

my $lineCount = 0;
my $writtenLines = 0;
my $marker = "";
my @indLabels = ();
my @randLabels = ();
my %groupColumns = ();
my %seenMarker = ();
my $printedLabels = 0; # false
while (<$genotypesInFile>){
    $lineCount++;
    my $line = $_;
    if($line =~ /^##/){
        # Determine individual labels
        # This works even if more than one <ID> region is present in the
        # header line, as might be the case in a 'join'ed file
        if($line =~ /IDs:\s+(.*?)\s*>/){
            @indLabels = ();
        }
        while($line =~ /IDs:\s+(.*?)\s*>/){
            my @lineData = split(/\s+/, $1);
            push(@indLabels, @lineData);
            # Remove first defined group, onto the next one
            $line =~ s/^.*?>//;
        }
    } else{
        if($line =~ /^(.*?)[\s,]/){
            $marker = $1;
            if($seenMarker{$marker}){
                printf(STDERR "\nWarning: Already seen data for marker %s, found again on line %d. ".
                       "This repeated line will be ignored\n", $marker, $lineCount);
            } else {
                $line =~ s/^.*?[\s,]+//;
                if(!$gtDataLength){
                    $gtDataLength = length($line);
                }
                # checks for data corruption
                if(length($line) != $gtDataLength){
                    # carry out a basic check for input data oddness
                    printf(STDERR "\nWarning: Length of line %d of input (%d) doesn't match the expected length (%d). ".
                           "This line will be ignored\n", $lineCount, length($line), $gtDataLength);
                } else {
                    if(!$printedLabels){
                        # if individual labels haven't been defined, make them up as 1 .. <n of individuals>
                        # this makes it easier to identify what individuals came from where in the split files
                        if(!@indLabels){
                            my @lineData = split(/\s+/, $line);
                            @indLabels = (1 .. (@lineData));
                        }
                        # determine random groupings of individuals
                        # first: generate internal column reference IDs
                        @randLabels = (0 .. (@indLabels - 1));
                        # second: shuffle IDs into a random order
                        @randLabels = shuffle_bukNew(@randLabels);
                        if($numSplits == 2){
                            if($countFirst eq ""){
                                $countFirst = scalar(@indLabels) * $ratioFirst;
                            }
                            if($countFirst > scalar(@indLabels)){
                                printf(STDERR "Error: Too many individuals requested in first split (%d). This is greater than \n".
                                       "       the total number of individuals (%d).\n", $countFirst, scalar(@indLabels));
                                usage();
                                exit(3);
                            }
                            # store two groups with desired split
                            $groupColumns{1} = [ @randLabels[(0 .. ($countFirst - 1))] ];
                            $groupColumns{2} = [ @randLabels[($countFirst .. (@indLabels - 1))] ];
                            printf(STDERR "COUNTFIRST = %d, NUMLABELS = %d...", $countFirst, scalar(@indLabels));
                        } else {
                            # store $numSplits groups with equal split
                            # the reverse order should make sure that every individual is included
                            my $currentPos = 0;
                            for(my $groupsRemaining = $numSplits; $groupsRemaining > 0; $groupsRemaining--){
                                my $nextCount = int((@indLabels - $currentPos) / $groupsRemaining);
                                $groupColumns{$groupsRemaining} = [ @randLabels[($currentPos .. ($currentPos + $nextCount - 1))] ];
                                $currentPos = $currentPos + $nextCount;
                            }
                        }
                        printf(STDERR "Splitting into %s groups of (", $numSplits);
                        for(my $i = 1; $i < $numSplits; $i++){
                            printf(STDERR "%d, ", scalar(@{ $groupColumns{$i} }));
                        }
                        printf(STDERR "%d) individuals\n", scalar(@{ $groupColumns{$numSplits} }));
                        # also, print labels into the output files
                        for(my $i = 0; $i < $numSplits; $i++){
                            print({$outFiles[$i]} "## <Individual/Column IDs: ");
                            print({$outFiles[$i]} join(" ", @indLabels[@{ $groupColumns{($i + 1)} }]));
                            print({$outFiles[$i]} " > ##\n");
                        }
                        $printedLabels = 1; # true
                        printf(STDERR "Writing output (one '.' per 1000 lines)");
                    }
                    $seenMarker{$marker} = 1; # true
                    $writtenLines++;
                    if($writtenLines % 1000 == 0){
                        print (STDERR ".");
                    }
                    my @lineData = split(/[\s,]+/, $line);
                    for(my $i = 0; $i < $numSplits; $i++){
                        # output marker to both files
                        print({$outFiles[$i]} $marker." ");
                        # output genotypes for respective individuals to both files
                        print({$outFiles[$i]} join(" ",@lineData[@{ $groupColumns{($ i +1)} }])."\n");
                    }
                }
            }
        }
    }
}

foreach(@outFiles){
    close($_);
}
printf (STDERR " done (%d lines written)!\n", $writtenLines);
