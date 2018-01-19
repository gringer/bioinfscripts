#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);

my $startName = "";
my $startPos = -1;
my $lastPos = -1;
my $lastDesc = "";

my $covAdjust = 5;
my $printedHeader = 0; # false
my $accuracyDP = 0;
my $DCOnly = 0; # false
my $threshold = 0; # Log2 threshold to report DC (otherwise DC is set to 0)
my $header = 1; # true

my $compare = "";

GetOptions("adjustment=i" => $covAdjust, "header!" => $header,
	   "compare=s" => \$compare, "dpaccuracy=i" => \$accuracyDP,
	   "threshold=f" => \$threshold, "onlydc!" => \$DCOnly ) or
    die("Error in command line arguments");

my @comp1 = ();
my @comp2 = ();
my @comps = ();

while(<>){
  chomp;
  if(!$_){
    next;
  }
  my ($refName, $pos, $refAllele, $cov, $bases, $qual, $rest) =
      split(/\t/, $_, 7);
  my $skip = ($bases =~ tr/<>//);
  my @adjCovs = ($cov - $skip);
  while($rest){
    ($cov, $bases, $qual, $rest) = split(/\t/, $rest, 4);
    $skip = ($bases =~ tr/<>//);
    push(@adjCovs, $cov - $skip);
  }
  if(!$printedHeader){
    ## set up comparisons
    if(!$compare){
      for(my $x = 0; $x < $#adjCovs; $x++){
	for(my $y = $x+1; $y <= $#adjCovs; $y++){
	  push(@comp1, $x);
	  push(@comp2, $y);
	  push(@comps, "$x,$y");
	}
      }
    } else {
      while($compare =~ s/([0-9]+?)[-,]([0-9]+?)[;,\s]?//){
	push(@comp1, $1);
	push(@comp2, $2);
	push(@comps, "${1}vs${2}");
      }
    }
    if($header){
      if($DCOnly){
        print(join("\t", "##Ref", "Start", "End", @comps)."\n");
      } else {
        print(join("\t", "##Ref", "Start", "End",
                   (1..($#adjCovs+1)), @comps)."\n");
      }
    }
    $printedHeader = 1; # true
  }
  my @DCov = ();
  for(my $i = 0; $i <= $#comp1; $i++){
    my $x = $comp1[$i]-1;
    my $y = $comp2[$i]-1;
    my $val = (log($adjCovs[$x]+$covAdjust) -
	       log($adjCovs[$y]+$covAdjust)) / log(2);
    push(@DCov, sprintf("%0.${accuracyDP}f",
			abs($val) < $threshold ? 0 : $val));
  }
  my $descLine = ($DCOnly) ?
      join("\t", @DCov) :
      join("\t", @adjCovs, @DCov);
  $descLine =~ s/-0/0/g;
  if(($refName ne $startName) || ($descLine ne $lastDesc)){
    # print sequence (if any)
    if($startName){
      print(join("\t", $startName, $startPos, $lastPos, $lastDesc)."\n");
    }
    $startName = $refName;
    $startPos = $pos;
    $lastDesc = $descLine;
  }
  $lastPos = $pos;
}

if($startName){
  print(join("\t", $startName, $startPos, $lastPos, $lastDesc)."\n");
}
