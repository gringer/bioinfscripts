#!/usr/bin/perl

## quantile_subset.pl -- filters out items from a sorted list (with
##   repeated values in the first field) that are the closest to the
##   desired quantiles

## For example, if there are 100 values, then a request for the 5th
## and 95th quantiles will produce the 5th and 95th value respectively.

## usage: cat <csv_file> | quantile_subset.pl -quantiles 0.05,0.95

use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);

my @values = ();
my $marker = "";

my @quantiles = (0.05,0.95);
my $quantileStr = "";

GetOptions("quantiles=s" => \$quantileStr) or
  die("Error in command line arguments");

if($quantileStr){
  @quantiles = split(/,/, $quantileStr);
}

while(<>){
  chomp;
  my @F = split(/,/, $_);
  my $lastLine = $_;
  if($marker ne $F[0]){
    if(@values){
      foreach my $quantile (@quantiles){
        if($quantile == 1){
          printf("%s\n", $values[$#values]);
        } else {
          my $rankPos = $quantile * $#values;
          my $rankInt = int($rankPos);
          my $rankFrac = $rankPos - $rankInt;
          ## standard quantile calculation doesn't make sense for
          ## fractional locations when actual position is in between
          ## two arbitrary text fields, so choose 0.5 as a threshold
          ## [e.g. SNP at chromosome 4, location 30 Mb vs
          ##       SNP at chromosome 7, location 10 Mb]
          if($rankFrac < 0.5){
            printf("%s\n", $values[$rankInt]);
          } else {
            printf("%s\n", $values[$rankInt+1]);
          }
        }
      }
      @values = ();
    }
  }
  $marker = $F[0];
  if(!/[0-9]$/){
    ## Write header line(s) [end with non-numeric character] out directly
    printf("%s\n", $lastLine);
  } else {
    push(@values, $lastLine);
  }
}

if(@values){
  foreach my $quantile (@quantiles){
    if($quantile == 1){
      printf("%s\n", $values[$#values]);
    } else {
      my $rankPos = $quantile * $#values;
      my $rankInt = int($rankPos);
      my $rankFrac = $rankPos - $rankInt;
      if($rankFrac < 0.5){
        printf("%s\n", $values[$rankInt]);
      } else {
        printf("%s\n", $values[$rankInt+1]);
      }
    }
  }
}
