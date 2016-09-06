#!/usr/bin/perl

# delta_stats.pl -- Calulate descriptive statistics and proportions
# for marker genotypes (assumes simplegt-formatted input file)

# Author: David Eccles (gringer), 2008 <programming@gringer.org>

# modified Oct 2010: genotype proportions, use expected proportions
# rather than observed proportions to reduce impact of missing data.
# modified Nov 2010: add count data output, allele frequency output

use warnings;
use strict;

sub usage {
    print(STDERR "usage: /delta_stats.pl < <input.file>\n");
    print(STDERR "\nOther Options:\n");
    print(STDERR "-help                   : display this information\n");
    print(STDERR "p<string> <int1> <int2> : define population as (<int1> .. <int2>)\n");
    print(STDERR "p<string> +<int>        : define population as (last+1 .. last+<int>)\n");
    print(STDERR "-nofreq                 : don't show frequencies\n");
    print(STDERR "-nonull                 : don't show missing proportion\n");
    print(STDERR "-count                  : show counts instead of frequencies\n");
    print(STDERR "-allelestats            : show allele (instead of genotype) statistics\n");
    print(STDERR "-csv                    : output comma-separated values\n");
    print(STDERR "\n");
}

sub max {
    my ($val1, $val2) = @_;
    if($val2 > $val1){
        return $val2;
    }
    else{
        return $val1;
    }
}

my %indivs = ();

# remove command line arguments (so while(<>) works)
my $lastPop = -1;
my $progname = $0;
my @pops = ();
my $showFreqs = 1; # true
my $showNull = 1; # true
my $showDelta = 1; # true
my $showCounts = 0; # false
my $alleleStats = 0; # false
my $csvFormat = 0; # false

my @fileNames = ();

while (@ARGV){
    my $argument = shift(@ARGV);
    if(-f $argument){
        push(@fileNames, $argument);
        printf(STDERR "Found input file '%s'\n", $argument);
    } else {
        if($argument eq "-help"){
            usage();
            exit(0);
        } elsif ($argument =~ /^p(.*)$/i){
            my $popname = $1;
            my $popFrom = shift(@ARGV);
            my $popTo = "";
            if($popFrom =~ /^\+([0-9]+)/){
                $popTo = $lastPop + $1;
                $popFrom = $lastPop + 1;
                $lastPop = $popTo;
            } else {
                $popTo = shift(@ARGV);
            }
            if($popFrom.$popTo =~ /^([0-9]+)$/){
                # Doing it this way allows additional columns to be
                # added to a pre-existing group
                if(!defined($indivs{$popname})){
                    push(@pops,$popname);
                    $indivs{$popname} = ();
                    printf(STDERR "Setting up population '%s'", $popname);
                } else {
                    printf(STDERR "Adding to population '%s'", $popname);
                }
                printf(STDERR " with range (%s .. %s)\n", $popFrom, $popTo);
                push (@{$indivs{$popname}}, ($popFrom .. $popTo));
            } else {
                printf(STDERR "Error: Population ranges for population '%s' (%s, %s) don't make sense\n",
                    $popname, $popFrom, $popTo);
            }
        } elsif ($argument =~ /^-count/){
            $showCounts = 1;
            printf(STDERR "Genotype counts will be displayed\n");
        } elsif ($argument =~ /^-nofreq/){
            $showFreqs = 0; #false
            printf(STDERR "Frequency/count statistics won't be displayed\n");
        } elsif ($argument =~ /^-allelestats/){
            $alleleStats = 1; #true
            printf(STDERR "Calculating allele (rather than genotype) frequencies/counts\n");
        } elsif ($argument =~ /^-csv/){
            $csvFormat = 1; #true
            printf(STDERR "Output will be CSV format\n");
        } elsif ($argument =~ /^-nonull/){
            $showNull = 0; #false
        } else {
            print(STDERR "Unknown argument '$argument'");
            usage();
            exit(1);
        }
    }
}

@ARGV = @fileNames;

my $nw = '%-0.7f '; #number width float formatting
my $nd = '%-9d '; #number width integer formatting
my $sw = '%-9s '; #number width string formatting
my $mw = '%-17s '; #marker width string formatting

if($csvFormat){
    $nw = ',%0.7f';
    $nd = ',%d';
    $sw = ',"%s"';
    $mw = '"%s"';
}

my $ncounts = 0;

# print out Column Titles

if (scalar(@pops) < 2){
    print(STDERR "Error: Fewer than two populations have been defined\n");
    usage();
    exit(1);
}

if (scalar(@pops) > 2){
    $showDelta = 0; # false
}

printf($mw, "Marker");
foreach my $pop1 (@pops){
    if($showFreqs){
        if($alleleStats){
            printf($sw, "p(A,$pop1)");
            printf($sw, "p(C,$pop1)");
        } else {
            printf($sw, "p(AA,$pop1)");
            printf($sw, "p(AC,$pop1)");
            printf($sw, "p(CC,$pop1)");
        }
    }
    if($showNull){
        printf($sw, "p(NN,$pop1)");
    }
}

if($showDelta){
    printf($sw, "delta");
}
print "\n";


my $pAA = 0;
my $pAC = 0;
my $pCC = 0;
my $pNN = 0;
my @popgts = ();

my @linePops = ();
my %gtCounts = ();

my %popgts = ();

while(<>){
    if(/^##/){
        next;
    }
    chomp;
    my ($marker, $genotypeLine) = split(/\s+/, $_, 2);
    printf($mw, $marker);
    if(!@linePops){
        # splits are slow, so only do this once
        @linePops = (0) x scalar(@_ = split(/\s+/, $genotypeLine));
        foreach my $pop1 (@pops){
            my @popCols = @{$indivs{$pop1}};
            @linePops[@popCols] = ($pop1) x scalar(@popCols);
        }
    }
    # replace complementary alleles and numeric alleles
    $genotypeLine =~ tr/1234aAcCgGtT/ACCAAACCCCAA/;
    # replace non-conforming alleles with N
    $genotypeLine =~ s/[^\sAC]/N/g;
    # make sure heterozygotes have the correct ordering
    $genotypeLine =~ s/CA/AC/g;
    %gtCounts = (); # reset count array
    foreach my $pop1 (@pops){
        # set individual counts to 0
        $gtCounts{$pop1}{"AA"} = 0;
        $gtCounts{$pop1}{"AC"} = 0;
        $gtCounts{$pop1}{"CC"} = 0;
        $gtCounts{$pop1}{"NN"} = 0;
    }
    # spin through line, storing counts for genotypes
    # this is the slowest thing in the code, so it should be heavily optimised
    my $col = 0;
    grep {
        ($gtCounts{$linePops[$col++]}{$_})++;
    } split(/ +/, $genotypeLine);
    # now work out frequencies
    foreach my $pop1 (@pops){
        # retrieve genotype list for populations individuals
        my $totalGoodCount = $gtCounts{$pop1}{"AA"} +
            $gtCounts{$pop1}{"AC"} + $gtCounts{$pop1}{"CC"};
        my $badCount = $gtCounts{$pop1}{"NN"};
        if($totalGoodCount){
            $pAA = $gtCounts{$pop1}{"AA"} / $totalGoodCount;
            $pAC = $gtCounts{$pop1}{"AC"} / $totalGoodCount;
            $pCC = $gtCounts{$pop1}{"CC"} / $totalGoodCount;
            $pNN = $badCount / ($badCount  + $totalGoodCount);
        } else {
            $pAA = 0; $pAC = 0; $pCC = 0; $pNN = 1;
        }
        if($showFreqs){
            if($showCounts){
                if($alleleStats){
                    # e.g. 10/60/30 -> 80/120
                    printf($nd,$gtCounts{$pop1}{"AA"} * 2 + $gtCounts{$pop1}{"AC"});
                    printf($nd,$gtCounts{$pop1}{"CC"} * 2 + $gtCounts{$pop1}{"AC"});
                } else {
                    printf($nd,$gtCounts{$pop1}{"AA"});
                    printf($nd,$gtCounts{$pop1}{"AC"});
                    printf($nd,$gtCounts{$pop1}{"CC"});
                }
            } else {
                if($alleleStats){
                    # e.g. 0.1/0.6/0.3 -> 0.4/0.6
                    printf($nw, $pAA + ($pAC / 2));
                    printf($nw, $pCC + ($pAC / 2));
                } else {
                    printf($nw, $pAA);
                    printf($nw, $pAC);
                    printf($nw, $pCC);
                }
            }
        }
        if($showNull){
            if($showCounts){
                printf($nd, $badCount);
            } else {
                printf($nw, $pNN);
            }
        }
    }
    if($showDelta){
        my $delta = -1;
        if(scalar(@pops) == 2){
            my $c0AA = $gtCounts{$pops[0]}{"AA"};
            my $c0AC = $gtCounts{$pops[0]}{"AC"};
            my $c0CC = $gtCounts{$pops[0]}{"CC"};
            my $c1AA = $gtCounts{$pops[1]}{"AA"};
            my $c1AC = $gtCounts{$pops[1]}{"AC"};
            my $c1CC = $gtCounts{$pops[1]}{"CC"};
            if((($c0AA + $c0AC + $c0CC) * ($c1AA + $c1AC + $c1CC)) != 0){
                # just allele delta for the moment
                my $p0A = ($c0AA * 2 + $c0AC) / (2 * ($c0AA + $c0AC + $c0CC));
                my $p1A = ($c1AA * 2 + $c1AC) / (2 * ($c1AA + $c1AC + $c1CC));
                $delta = abs($p0A - $p1A);
            }
        }
        printf($nw, $delta);
    }
    print("\n");
}
