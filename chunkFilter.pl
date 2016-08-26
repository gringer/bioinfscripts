#!/usr/bin/perl

use warnings;
use strict;

my $mw = 10; # match window
my $maxChunkSize = 5000; # limit chunk size

my $lastContig="";
my $start=-($mw+1);
my $length=-($mw+1);

my @lineBuffer = ();

while(<>){
  chomp;
  my @F = split(/ /);
  my $newStart = $F[1];
  my $newEnd = $F[1] + $F[2];
  my $oldEnd = $start + $length;
  if(($F[0] eq $lastContig) && ($length < ($maxChunkSize-$mw)) &&
     (($start - $mw) < $newStart && ($oldEnd + $mw) > $newStart) ||
     (($oldEnd + $mw) > $newEnd && ($start - $mw) < $newEnd)){
    $start = $newStart if ($start > $newStart);
    $newEnd = $oldEnd if ($oldEnd > $newEnd);
    $F[1] = $start;
    $F[2] = $newEnd - $start;
  } else {
    push(@lineBuffer,"");
    my $sig = join(" ", ($lastContig, $start, $length));
    print(join(" $sig\n", @lineBuffer));
    @lineBuffer = ();
    $lastContig = "";
    $start = -($mw+1);
    $length = -($mw+1);
  }
  push(@lineBuffer, $_);
  ($lastContig, $start, $length) = @F;
}

my $sig = join(" ", ($lastContig, $start, $length));
push(@lineBuffer,"");
print(join(" $sig\n", @lineBuffer));
