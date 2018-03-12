#!/usr/bin/perl

# stl2svg.pl -- create an SVG file from a top-down view of an STL object
# usage: ~/scripts/stl2svg.pl input.stl > output.svg

# Note: There is no orientation / camera setting. This script just takes
#       a view from the top. Use OpenSCAD (or similar) to re-orient STL
#       files so that the top-down view is the desired view.

use warnings;
use strict;

use Math::Trig qw(:pi);

my @pts = ();
my $lastX = undef;
my $lastY = undef;
my $cull = 0; # false
my $zAng = undef;
my $colourVal = undef;
my $ns = undef;

my %lines = ();

my %polys = ();
my %normColours = ();
my %colours = ();
my $minX = undef;
my $minY = undef;
my $minZ = undef;
my $maxX = undef;
my $maxY = undef;
my $maxZ = undef;


while (<>) {
  chomp;
  s/^\s+//;
  if (/^facet normal (([0-9.\-e]+) ([0-9.\-e]+) ([0-9.\-e]+))/) {
    $ns = $1;
    my ($nx, $ny, $nz) = ($2, $3, $4);
    $ny = -$ny;
    $zAng = atan2($ny,$nx);
    if ((abs($nx) < 0.01) && (abs($ny) < 0.01)) {
      $normColours{$ns} = sprintf("%02x", 240);
    } else {
      $normColours{$ns} = sprintf("%02x",(40 + (abs($zAng) / pi) * 200));
    }
    @pts = ();
    $lastX = undef;
    $lastY = undef;
    $minZ = undef;
    $maxZ = undef;
    if (($nz) < 0.01) {
      $cull = 1;                # true
    } else {
      $cull = 0;
    }
  }
  if (/^endfacet/) {
    if (!$cull) {
      for(my $i = 0; $i < scalar(@pts); $i++){
        my $pt1 = $pts[$i];
        my $pt2 = $pts[($i+1)%(scalar(@pts))];
        if($pt2 lt $pt1){
          ($pt2, $pt1) = ($pt1, $pt2);
        }
        $lines{$ns}{sprintf("%s %s", $pt1, $pt2)}++;
      }
    }
  }
  if (/^vertex (.*$)/) {
    my ($x, $y, $z) = split(/ /, $1);
    $y = -$y;
    if (!defined($minX) || ($minX > $x)) {
      $minX = $x;
    }
    if (!defined($minY) || ($minY > $y)) {
      $minY = $y;
    }
    if (!defined($minZ) || ($minZ > $z)) {
      $minZ = $z;
    }
    if (!defined($maxX) || ($maxX < $x)) {
      $maxX = $x;
    }
    if (!defined($maxY) || ($maxY < $y)) {
      $maxY = $y;
    }
    if (!defined($maxZ) || ($maxZ < $z)) {
      $maxZ = $z;
    }
    push(@pts, sprintf("%s,%s,%s", $x, $y, $z));
  }
}

foreach my $normal (keys(%lines)) {
  my %newLines = ();
  foreach my $lineStr (keys(%{$lines{$normal}})){
    my $count = $lines{$normal}{$lineStr};
    if($count == 1){
      my ($pt1, $pt2) = split(/ /, $lineStr);
      if(!exists($newLines{$pt1})){
        $newLines{$pt1} = ();
      }
      if(!exists($newLines{$pt2})){
        $newLines{$pt2} = ();
      }
      push(@{$newLines{$pt1}}, $pt2);
      push(@{$newLines{$pt2}}, $pt1);
      #printf("  %s\n", $lineStr);
    }
  }
  my @startPts = keys(%newLines);
  while(@startPts){
    my $lastPt = shift(@startPts);
    my @polyPts = ();
    my $pathMinZ = undef;
    my $pathMaxZ = undef;
    my $lastX = 0;
    my $lastY = 0;
    while(@{$newLines{$lastPt}}){
      my ($x, $y, $z) = split(/,/, $lastPt);
      if(!defined($pathMinZ) || $pathMinZ > $z){
        $pathMinZ = $z;
      }
      if(!defined($pathMaxZ) || $pathMaxZ < $z){
        $pathMaxZ = $z;
      }
      push(@polyPts, sprintf("%s,%s", $x-$lastX, $y-$lastY));
      ($lastX, $lastY) = ($x, $y);
      my $newPt = shift(@{$newLines{$lastPt}});
      # have removed $lastPt -> $newPt, but also need to remove $newPt -> $lastPt
      my @pts2 = grep {$_ ne $lastPt} @{$newLines{$newPt}};
      $newLines{$newPt} = \@pts2;
      $lastPt = $newPt;
    }
    if(@polyPts){
      my $midZ = $pathMaxZ + $pathMinZ;
      if(!exists($polys{$midZ})){
        $polys{$midZ} = ();
        $colours{$midZ} = ();
      }
      push(@{$polys{$midZ}}, join(" ", @polyPts));
      push(@{$colours{$midZ}}, $normColours{$normal});
    }
  }
}

my $figWidth = ($maxX - $minX);
my $figHeight = ($maxY - $minY);
my $svgWidth = sprintf("%0.0f", $figWidth * 1.05);
my $svgHeight = sprintf("%0.0f", $figHeight * 1.05);
my $ox = $figWidth * 0.025 - $minX;
my $oy = $figHeight * 0.025 - $minY;

printf("<svg xmlns:svg=\"http://www.w3.org/2000/svg\"
   xmlns=\"http://www.w3.org/2000/svg\"
   preserveAspectRatio=\"xMidYMid meet\"
   viewBox=\"0 0 $svgWidth $svgHeight\"
   version=\"1.1\">\n");
printf(" <g fill=\"none\" ".
       "stroke-width=\"0.03\" stroke-linejoin=\"round\">\n");
foreach my $zPos (sort {$a <=> $b} (keys(%polys))) {
  foreach my $line (@{$polys{$zPos}}) {
    my ($fx, $fy, $rest) = split(/[, ]/, $line, 3);
    my $colour = shift(@{$colours{$zPos}});
    $colour = "#$colour$colour$colour";
    $fx += $ox;
    $fy += $oy;
    printf("  <path d=\"m%0.4f,%0.4f %s z\" fill=\"%s\" stroke=\"%s\"/>\n",
           $fx, $fy, $rest, $colour, $colour);
  }
}
printf(" </g>\n");
printf("</svg>\n");
