#!/usr/bin/perl

# rsfilter.pl -- hunts (quickly) for a set of markers in one pass of a
# file.

# Note: marker names *must not* be the same as an existing file

# Author: David Eccles (gringer), 2008 <programming@gringer.org>

# This program tries to be fast and memory-efficient, so decisions are
# pushed to outer loops. Using the ordered option, '-o', retains
# genotype lines in memory for the specified markers until all input
# has been read, which can consume large amounts of memory when the
# requested marker set is large.

use strict;
use warnings;

use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

sub usage {
  print("usage: ./rsfilter.pl <marker1> <marker2> ... [options] < <file name>\n");
  print("\nOther Options:\n");
  print("-r : invert filter (i.e. select markers to exclude)\n");
  print("-o : order by selection\n");
  print("-- : no more files on command line (treat file names as markers)\n");
  print("\n");
}

my %markers = ();
my $complement = 0; # false # determines if complement, rather
                            # than intersect, should be chosen
my $order = 0; # false # order markers by given ordering
my @markerorder = ();

my @inFiles = ();

my $filesFinished = 0; # false

while(@ARGV){
    my $arg = shift(@ARGV);
    if(!$filesFinished && (-e $arg)){
        if(!keys(%markers)){
            print(STDERR "Retrieving marker names from $arg...");
            my $markFile = 0;
            my $fileName = $arg;
            $markFile = new IO::Uncompress::Gunzip "$fileName" or
                die "Unable to open $fileName\n";
            while(<$markFile>){
                if(/^\"?(.*?)\"?[\s,]+/){
                    my $marker = $1;
                    # print(STDERR "adding marker $marker\n");
                    $markers{$marker} = 1;
                    push(@markerorder, $marker);
                }
            }
            print(STDERR keys(%markers)." marker names extracted\n");
        } else {
            push(@inFiles, $arg);
        }
    } elsif ($arg =~ /^-r/){
        $complement = 1; # true
        print(STDERR "Markers will be excluded, rather than filtered\n");
    } elsif ($arg =~ /^-o/){
        $complement = 0; # false
        $order = 1; # true
        print(STDERR "Markers will output in specified order\n");
    } elsif ($arg =~ /^--/){
        $filesFinished = 1; # true
        print(STDERR "No more files will be read\n");
    } elsif ($arg =~ /^-help/){
        usage();
        exit(0);
    } else {
        $markers{$arg} = 1;
        push(@markerorder, $arg);
    }
}

@ARGV = @inFiles;

my %markerlines = ();

#print keys(%markers)." marker names extracted\n";

if($complement){ # tests are in outer loop to slightly reduce processor effort
    while (<>){
        my $line = $_;
        if($line =~ /^(\"?.*?\"?)[\s,]+/){
            if (!$markers{$1}){
                print($line);
            }
        }
    }
} else {
    if(!$order){
        while (<>){
            my $line = $_;
	    if($line =~ /^(\"?.*?\"?)[\s,]+/){
                if ($markers{$1}){
                    print($line);
                }
            }
        }
    } else {
        while (<>){
            my $line = $_;
	    if($line =~ /^(\"?.*?\"?)[\s,]+/){
                if ($markers{$1}){
                    my $marker = $1;
                    $markerlines{$1} = $line;
                }
            }
        }
        foreach(@markerorder){
            my $marker = $_;
            if($markerlines{$marker}){
                print($markerlines{$marker});
            }
        }
    }
}
