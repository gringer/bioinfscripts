#!/usr/bin/perl

use warnings;
use strict;

my $read = "";
my @tigs = ();
my %links = ();
while(<>){
  my @F = split(/\s/);
  if($F[0] eq $read){
    push(@tigs, $F[5]);
  } else {
    foreach my $t1 (@tigs){
      foreach my $t2 (@tigs){
        $links{$t1}{$t2}++;
      }
    }
    $read = $F[0];
    @tigs = ();
  }
}

foreach my $t1 (@tigs){
  foreach my $t2 (@tigs){
    $links{$t1}{$t2}++;
  }
}

print("tig1,tig2,count\n");
foreach my $t1 (sort(keys(%links))){
  foreach my $t2 (sort(keys(%{$links{$t1}}))){
    my $line = sprintf("%s,%s,%d\n", $t1, $t2, $links{$t1}{$t2});
    $line =~ s/Consensus_//g;
    print $line;
  }
}
