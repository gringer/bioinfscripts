#!/usr/bin/perl

use warnings;
use strict;

use List::Util;

sub processCounts{
  my ($marker, $countRef) = @_;
  if($marker eq ""){
    return;
  }
  #print(STDERR "Processing $marker\n");
  foreach my $sex (keys(%{$countRef})){
    my %tCounts = %{$countRef->{$sex}};
    my $total = 0;
    my $minA = "0";
    my $minC = undef;
    while(my ($a, $b) = each %tCounts) {
      if(!defined($minC) || ($b < $minC)){
        $minC = $b;
        $minA = $a;
      }
      $total += $b;
    }
    my $minF = ($total > 0) ? $minC / $total : 0;
    printf("%s,%s,%s,%0.4f\n", $marker, $sex, $minA, $minF);
  }
}

my $lastMarker = "";
my $counts = {};

while(<>){
  chomp;
  #print $_."\n";
  my ($marker, $sex, $type, $bases, $count) = split(/,/);
  if($marker ne $lastMarker){
    processCounts($lastMarker, $counts);
    $counts = {"m" => {}, "f" => {}, "a" => {}};
    $lastMarker = $marker;
  }
  #print(join(";",%{$counts})."\n");
  if(($type ne "a") || ($bases eq "0")){
    next;
  }
  $counts->{$sex}->{$bases} = $count;
}
processCounts($lastMarker, $counts);
