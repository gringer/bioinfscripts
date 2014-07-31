#!/usr/bin/perl

# indiv_minimal.pl -- Make sure all individuals are different using
# a minimal set of information from the genotyped individuals.
# Practically, this means extracting additional lines from an input
# file until every individual appears different (or the total input
# file, whichever comes first)

# Author: David Eccles (gringer), 2007 <programming@gringer.org>

use strict;
use warnings;

sub usage {
  print("usage: ./indiv_minimal.pl <simplegt / tped file>\n");
  print("\nselects markers until every individual looks different\n");
  print("\nOther Options:\n");
  print("-help        : Only display this help message\n");
  print("\n");
}

my @files = ();

# extract command line arguments
while(@ARGV){
    my $argument = shift @ARGV;
    if(-f $argument){ # file existence check
	push(@files, $argument);
    } else {
        if($argument eq "-help"){
            usage();
            exit(0);
        }
    }
}

@ARGV = @files;


my $nw = '%-0.5f '; #number width string formatting
my $sw = '%-7s '; #number width string formatting

my $nmarkers = 0;
my $informativeMarkers = 0;
my @diffcounts = ();
my @indLabels = ();

my $filetype = "";
my $oldIDP = -1;


while(<>){
    chomp;
    my $line = $_;
    if(!$filetype){
        if($line =~ /IDs:\s+(.*?)\s*>/){
            print(STDERR "simplegt format detected\n");
            $filetype = "simplegt";
            @indLabels = ();
            @diffcounts = ();
            my $tmpLine = $line;
            while($tmpLine =~ /IDs:\s+(.*?)\s*>/){
                my @lineData = split(/\s+/, $1);
                push(@indLabels, @lineData);
                # Remove first defined group, onto the next one
                $tmpLine =~ s/^.*?>//;
            }
            for(my $c1=0;$c1<(@indLabels-1);$c1++) {
                for(my $c2=$c1+1;$c2<(@indLabels+0);$c2++) {
                    $diffcounts[$c1*(@indLabels+0)+$c2] = 0;
                }
            }
        } elsif($line =~ /^(chr)?([0-9]{1,2}|[XY])\s/){
            print(STDERR "tped format detected\n");
            $filetype = "tped";
        } else {
            print(STDERR "simplegt format assumed\n");
            $filetype = "simplegt";
        }
    }
    if($line =~ /^#/){
        next;
    }
    my $marker = "";
    my @genotypes = ();
    if($filetype eq "tped"){
        my @fields = split(/\s/,$line);
        shift(@fields); # chromosome
        $marker = shift(@fields);
        shift(@fields); # centimorgan location
        shift(@fields); # base-pair location
        while(@fields){
            push(@genotypes, shift(@fields).shift(@fields));
        }
    } else {
        ($marker, @genotypes) = split(/\s/,$line);
    }
    $nmarkers++;

    if(!@indLabels){
        if($ARGV =~ /.tped$/){
            print(STDERR "Reading tfam data\n");
            my $tfamFileName = $ARGV;
            $tfamFileName =~ s/.tped$/.tfam/;
            open(my $tfamFile, "<", $tfamFileName) or die("Cannot open $tfamFileName");
            while(<$tfamFile>){
                my @fields = split(/\s/,$_);
                push(@indLabels, $fields[0]."_".$fields[1]);
            }
            close($tfamFile);
        } else {
            print(STDERR "Generating dummy individual IDs\n");
            @indLabels = (1 .. scalar(@genotypes));
        }
        @diffcounts = ();
        for(my $c1=0;$c1<(@indLabels-1);$c1++) {
            for(my $c2=$c1+1;$c2<(@indLabels+0);$c2++) {
                $diffcounts[$c1*(@indLabels+0)+$c2] = 0;
            }
        }
    }

    # substitute complementary alleles to simplify comparisons later
    for(my $c1=0;$c1<(@indLabels+0);$c1++) {
	$genotypes[$c1] =~ s/(a|t)/t/ig;
	$genotypes[$c1] =~ s/(g|c)/c/ig;
    }
    for(my $c1=0;$c1<(@indLabels-1);$c1++) {
 	for(my $c2=$c1+1;$c2<(@indLabels+0);$c2++) {
 	    my $gt1 = $genotypes[$c1];
 	    my $gt2 = $genotypes[$c2];
	    #same genotype (assuming a=t,c=g)
 	    if ($gt1 eq $gt2){
		$diffcounts[$c1*(@indLabels+0)+$c2]++;
 	    }
	    #same genotype, but switched around
 	    elsif (($gt1 eq "ac" && $gt2 eq "ca") ||
		   ($gt1 eq "ca" && $gt2 eq "ac")){
		$diffcounts[$c1*(@indLabels+0)+$c2]++;
 	    }
 	    # if genotypes are not similar, then don't consider similar
 	    else {
		# @diffcounts[$c1*(@indLabels+0)+$c2] += 0;
 	    }
 	}
    }
    my $different = 1; # true
    my $identicalPairs = 0;
    for(my $c1=0;$c1<(@indLabels-1);$c1++) {
 	for(my $c2=$c1+1;$c2<(@indLabels+0);$c2++) {
	    if ($diffcounts[$c1*(@indLabels+0)+$c2] == $nmarkers){
		$different = 0; # false
                $identicalPairs++;
	    }
 	}
    }
    if($identicalPairs != $oldIDP){
        $informativeMarkers++;
        printf(STDERR "%s: %d identical pairings\n", $marker, $identicalPairs);
        printf("%s\n", $marker);
        $oldIDP = $identicalPairs;
    }
    if ($different){
	last;
    }
}

# Header line 1, shows number of individuals
printf("Number of individuals: %d\n", (@indLabels+0));

# Header line 2, shows number of markers
printf("Number of markers: %d (%d informative)\n", $nmarkers, $informativeMarkers);

print("Pairwise comparisons:\n");
for(my $c1=0;$c1<(@indLabels-1);$c1++) {
    for(my $c2=$c1+1;$c2<(@indLabels+0);$c2++) {
	printf("%10s %10s %s\n", $indLabels[$c1], $indLabels[$c2], ($diffcounts[$c1*(@indLabels+0)+$c2] / $nmarkers));
    }
}
