#!/usr/bin/perl

# snpchip2linkage.pl -- converts from a simplegt formatted file into a
# linkage formatted (PED) file.

# Author: David Eccles (gringer), 2007 <programming@gringer.org>

# updated 2008 to accept individual label header
# updated 2010 to accept MAP file, affected/unaffected change and more
#    standard argument parsing

use warnings;
use strict;
use FileHandle;

sub usage {
    print(STDERR "usage: ./snpchip2linkage.pl <marker location file> ".
        "<genotype file>\n\n");
  print(STDERR "\nOther Options:\n");
  print(STDERR "-help          : display this information\n");
  print(STDERR "-uaf <integer> : Define last individual with unaffected status\n");
  print(STDERR "-aff <integer> : Define last individual with affected status\n");
}

my %genogroup = ();

my $numgts = 0;
my %markers = ();
my @markerOrder = ();
my $markersFound = 0; # false

my @inputFiles = ();

my %includedMarkers = ();

my $uafChange = 0; # false
my $affChange = 0; # false


while(@ARGV){
    my $arg = shift(@ARGV);
    if(-f $arg){
        if(!$markersFound){
            print(STDERR "Reading in marker file... ");
            my $markerColumn = -1;
            open MARKFILE, "< $arg"
                or die("cannot open $arg for reading");
            while(<MARKFILE>){
                my @data = split(/\s+/, $_);
                if($data[0] !~ /^#/){
                    if($markerColumn == -1){
                        for(my $i = 0; $i < scalar(@data); $i++){
                            if($data[$i] =~ /(snp|rs)/i){
                                $markerColumn = $i;
                            }
                        }
                    }
                    $markers{$data[$markerColumn]} = 1;
                    push(@markerOrder,$data[$markerColumn]);
                }
            }
            $markersFound = 1; # true;
            print(STDERR "done!\n");
        } else {
            push(@inputFiles, $arg);
        }
    } else {
        if($arg eq "-help"){
            usage();
            exit(0);
        } elsif($arg eq "-aff"){
            $affChange = shift(@ARGV);
            print(STDERR "Status will change from affected to unaffected after ".
                  "individual $affChange\n");
        } elsif($arg eq "-uaf"){
            $uafChange = shift(@ARGV);
            print(STDERR "Status will change from unaffected to affected after ".
                  "individual $uafChange\n");
        } else {
            print(STDERR "Unknown argument '$arg'");
            usage();
            exit(1);
        }
    }
}

if(!$markersFound){
    print STDERR "Error: no marker files specified on command line\n";
    usage();
    exit(2);
}

if($affChange && $uafChange){
    print STDERR "Error: Please select only '-aff' or '-uaf'\n";
    usage();
    exit(3);
}

@ARGV = @inputFiles;

my @indLabels = ();
my @phenoVals = ();
my $phenoNum = 1;

print(STDERR "Reading in input... ");
while (<>){
    my $line = $_;
    if($line =~ /^##/){
        ## Determine individual labels
        ## This works even if more than one <ID> region is present in the
        ## header line, as might be the case in a 'join'ed file
        if($line =~ /IDs:\s+(.*?)\s*>/){
            @indLabels = ();
            @phenoVals = ();
        }
        while($line =~ /IDs:\s+(.*?)\s*>/){
            my @lineData = split(/\s+/, $1);
            push(@indLabels, @lineData);
            ## generate an array of (@lineData) copies of $phenoNum
            push(@phenoVals, (($phenoNum) x scalar(@lineData)) );
            $phenoNum++;
            $line =~ s/^.*?>//;
        }
    } else{
        my ($marker, $rest) = split(/\s+/, $line, 2);
        $rest = lc($rest);
        $rest =~ s/\s+/ /g; #replace whitespace with spaces
        $rest =~ tr/acgt/1234/;
        $rest =~ s/[^0-9 ]/0/g;
        $rest =~ s/ ([34])([12])/ $2$1/ig; # order heterozygotes so
                                           # smallest number is first
        my @genotypes = split(/\s+/,$rest);
        if(defined($markers{$marker})){
            $includedMarkers{$marker} = 1;
            @{$genogroup{$marker}} = @genotypes;
            if($numgts && ($numgts != (@genotypes))){
                die("Number of genotypes does not match");
            }
            else{
                $numgts = (@genotypes);
            }
        }
    }
}
print(STDERR "done!\n");

if(scalar(keys(%markers)) > scalar(keys(%includedMarkers))){
    printf(STDERR "Warning: map file contains more markers (%d) ".
           "than have been retrieved from simplegt file (%d).\n",
        scalar(keys(%markers)), scalar(keys(%includedMarkers)));
}

print(STDERR "Creating individual labels... ");
if(@indLabels){
    if(scalar(@indLabels) != $numgts){
        die("Number of genotypes does not match labels in header line");
    }
} else {
    @indLabels = (1 .. $numgts);
    my $rest = scalar(@indLabels) - ($affChange + $uafChange);
    if($affChange){
        @phenoVals = ((2) x $affChange);
        push(@phenoVals, (1) x $rest);
    } elsif($uafChange){
        @phenoVals = ((1) x $uafChange);
        push(@phenoVals, (2) x $rest);
    } else {
        @phenoVals = ((0) x $rest);
    }
}
print(STDERR "done!\n");

print(STDERR "Writing PED file..");
for(my $i=0; $i < $numgts; $i++){
    if($i % 100 == 0){
        print(STDERR ".");
    }
    my @genotypes = ();
    foreach my $marker (@markerOrder){
        if($includedMarkers{$marker}){
            my $writeval = @{$genogroup{$marker}}[$i];
            $writeval =~ s/(.)(.)/$1 $2/;
            push(@genotypes, $writeval);
        }
    }
    printf("%s 1 0 0 0 %d  ", $indLabels[$i], $phenoVals[$i]);
    print join("  ",@genotypes)."\n";
}
print(STDERR " done!\n");
