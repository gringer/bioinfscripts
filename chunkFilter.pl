#!/usr/bin/perl

use warnings;
use strict;

my $mw = 10; # match window

my $lastContig="";
my $start=-($mw+1);
my $length=-($mw+1);

my @lineBuffer = ();

while(<>){
  chomp;
  my @F = split(/ /);
  if(($F[0] eq $lastContig) &&
     ((abs($start-$F[1]) <= $mw) || (abs($length-$F[2]) <= $mw))){
    $F[1] = $start if $start < $F[1];
    $F[2] = $length if $length > $F[2];
  } else {
    push(@lineBuffer,"");
    my $sig = join(".", ($lastContig, $start, $length));
    print(join(" | $sig \n", @lineBuffer));
    @lineBuffer = ();
    $lastContig="";
    $start=-($mw+1);
    $length=-($mw+1);
  }
  push(@lineBuffer, $_);
  ($lastContig, $start, $length)=@F;
}

my $sig = join(".", ($lastContig, $start, $length));
print(join(" | $sig \n", @lineBuffer));
