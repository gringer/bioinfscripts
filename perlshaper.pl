#!/usr/bin/perl

# perlshaper.pl -- convert shapefile data to SVG format, for use in
# wikimedia maps.

# This code is designed with wikimedia in mind, so defaults should
# produce SVG images that follow the general map conventions (see
# http://en.wikipedia.org/wiki/Wikipedia:WikiProject_Maps/Conventions). In
# particular, the expected input format is that of Natural Earth
# files (http://www.naturalearthdata.com).

# Author: David Eccles (gringer) 2010-2013 <programming@gringer.org>

#  -- Begin GPL license blurb --

# Copyright 2010-2013 David Eccles (gringer)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#  -- End GPL license blurb --

# The most recent version of this code can be found at
# http://en.wikipedia.org/wiki/User:Gringer/perlshaper

# Just as a head's up, latitude lines are horizontal in the typical
# projection of the world, while longitude lines are
# vertical. Latitude is usually represented by a number (in degrees)
# in the range -90..90, while longitude is usually represented by a
# number (in degrees) in the range -180..180.
#
#  /---------\  --- latitude (Y)
# | | | | | | |
# |-|-|-|-|-|-|
# | | | | | | |
#  \---------/
#     |
#     \-- longitude (X)

use warnings;
use strict;

use Math::Trig;
# http://search.cpan.org/~jasonk/Geo-ShapeFile-2.52/lib/Geo/ShapeFile.pm
use Geo::ShapeFile;
use SVG;

my $progVersion = 1.83;

sub usage {
  print(STDERR qq{usage: ./perlshaper.pl <shapefile(s)> [options] > output.svg
\n---- Basic Options ----
-help                 : Only display this help message
-type <string>        : Type of map (location|locator|area|world|orthographic)
-res <float>          : Resolution (in adjusted degrees)
-round <int>          : Round values to given number of decimal places
-centre <long>,<lat>  : Identify centre of map (by longitude/latitude)
-centre <ISO 3A-code> : Identify centre of map (by target country)
-data <file>          : Colour countries based on numerical data in file
\n---- Advanced Options ----
-psize <float>        : Radius of small points
-colour <string>      : Change colour theme to a different map type
-proj <string>        : Change target projection
-landcol <string>     : Land colour
-seacol <string>      : Sea colour
-bordcol <string>     : Border colour
-sub <ISO 3A-code>    : Identify subject region
-pol <ISO 3A-code>    : Identify related political region
-only <ISO 3A-code>   : Only display specified shape(s)
-zoom <ISO 3A-code>   : Zoom to projected extents of shape(s)
-nokey                : Don't display heatmap key (if heatmap is used)
-[no]lines            : [don't] print lines of latitude and longitude
-v                    : Verbose output (also '-vv')
});
}

# Transform a point from {[-180,180],[-90,90]} to
# {[-0.5,0.5],[-0.5,0.5]} relative to the desired projection. This
# replaces the proj4 transform function with something that will
# always return a valid point. Orthographic transformations that
# appear on the opposite side of the globe will be converted to a
# point with the same angle, but sitting at the circle
# edge. Adjustments need to be made post-transform (e.g. multiply by
# SVG width, add to fit in 0..max range) to convert to SVG
# coordinates.

sub transform {
  my ($options, $inPoint) = @_;
  my $oldLong = $inPoint->[0];
  my $oldLat = $inPoint->[1];
  my $projection = $options->{"projection"};
  my $lambda = ($oldLong - $options->{"centreLn"});
  while($lambda <= -180){
    $lambda += 360;
  }
  while($lambda > 180){
    $lambda -= 360;
  }
  my $include = 1;
  my $phi = $oldLat;
  my $phi1 = $options->{"centreLt"};
  my ($x, $y);
  if ($projection eq "equirectangular") {
    # see http://en.wikipedia.org/wiki/Equirectangular_projection
    $x = ($lambda * cos($phi1 * pi / 180)) / 360;
    $y = ($phi / 180);
  } elsif ($projection eq "orthographic") {
    # see http://en.wikipedia.org/wiki/Orthographic_projection_(cartography)
    $phi = $phi * pi / 180;
    $phi1 = $phi1 * pi / 180;
    $lambda = $lambda * pi / 180;
    $x = cos($phi) * sin($lambda);
    $y = cos($phi1) * sin($phi) - sin($phi1) * cos($phi) * cos($lambda);
    my $theta = atan2($y,$x);
    my $cosc = sin($phi1) * sin($phi) +
      cos($phi1) * cos($phi) * cos($lambda);
    if (($cosc) < 0) {
      $include = 0;
      $x = cos($theta);         # put on edge of map (i.e. r = 1)
      $y = sin($theta);
    }
    if ($options->{"rotation"}) {
      # there may be a quicker way to do this above...
      my $r = sqrt($x*$x + $y*$y);
      $theta += $options->{"rotation"} * pi / 180;
      $x = $r * cos($theta);
      $y = $r * sin($theta);
    }
    # adjust by a factor of 0.5 to fit in -0.5..0.5 that other projections use
    $x *= 0.5;
    $y *= 0.5;
  } elsif ($projection eq "wintri") {
    $phi = $phi * pi / 180;
    my $cphi1 = (2 / pi);
    $lambda = $lambda * pi / 180;
    # see http://en.wikipedia.org/wiki/Winkel_Tripel
    my $alpha = acos(cos($phi) * cos($lambda / 2));
    my $sincalpha = ($alpha == 0) ? 1 : (sin($alpha) / $alpha);
    $x = (($lambda * $cphi1) +
             (2 * cos($phi) * sin($lambda / 2) / $sincalpha)) / (4 + 2 * pi);
    $y = ($phi + sin($phi) / $sincalpha) / (2 * pi);
  } else {
    die("Error: projection must be one of the following: ".
        "equirectangular, wintri, orthographic");
  }
  return([$x, $y, $include]);
}

# Convert from standard equirectangular projection ("+proj=latlong
# +datum=WGS84" in proj4-speak) to another projection. Natural Earth
# Data uses this equirectangular format, where coordinates are
# specified based on their latitude and longitude locations in
# degrees.

# This method attempts to fit the map into a rectangle (or square)
# that fits in a 1100x550 box.

sub project {
  my ($options, $pointRef) = @_;
  my @input = @{$pointRef};
  my $newProj = $options->{"projection"};
  my $projWidth = $options->{"svgWidth"};
  my $projHeight = $options->{"svgHeight"};
  my $padding = $options->{"padding"};
  my $xAdj = $options->{"xAdj"};
  my $yAdj = $options->{"yAdj"};
  my ($minX, $minY) = (1.5 - $xAdj, 1.5 - $yAdj);
  my ($maxX, $maxY) = ($minX + $projWidth + $padding * 2,
                       $minY + $projHeight + $padding * 2);
  if(defined($options->{"minX"})){
    ($minX, $maxX) = ($options->{"minX"}, $options->{"maxX"});
    ($minY, $maxY) = ($options->{"minY"}, $options->{"maxY"});
  }
  my $xScale = $options->{"xScale"};
  my $yScale = $options->{"yScale"};
  my @output = ();
  my ($oldX, $oldY);
  my $oldLat;
  foreach my $inPoint (@input) {
    my $newLong = 0;
    my $newLat = 0;
    if (UNIVERSAL::can($inPoint,'isa')) {
      $newLong = $inPoint->X;
      $newLat = $inPoint->Y;
    } else {
      $newLong = $inPoint->[0];
      $newLat = $inPoint->[1];
    }
    my $pt = transform($options, [$newLong, $newLat]);
    my $px = ($pt -> [0]) * $xScale * $projWidth;
    # Y inverted because SVG file is inverted
    my $py = ($pt -> [1]) * $yScale * -$projHeight;
    # transformed points should fit in the box (0,0)-(width,height) after adjustment
    if(($px + $xAdj < $minX) || ($px + $xAdj > $maxX) ||
       ($py + $yAdj < $minY) || ($py + $yAdj > $maxY)){
      # $pt->[2] = 0;
    }
    $oldX = $px if !defined($oldX);
    $oldY = $py if !defined($oldY);
    $oldLat = $newLat if !defined($oldLat);
    my ($xd, $yd) = ($px - $oldX, $py - $oldY);
    if (sqrt($xd * $xd + $yd * $yd) > $projWidth / 2) {
      # don't connect lines that have wrapped around the map (over half the image width)
      if(($options->{"zoomed"}) ||
         ($options->{"projection"} eq "orthographic")){
        # zoomed and orthographic projections don't get additional points added
        push(@output, 0);
      } else {
        # add additional points on the border edge
        my ($minLong, $maxLong) = ($options->{"minLong"}, $options->{"maxLong"});
        my $oldEdgePoint = transform($options, [($px > $oldX) ? $minLong : $maxLong, $oldLat]);
        my $newEdgePoint = transform($options, [($px > $oldX) ? $maxLong : $minLong, $newLat]);
        $oldEdgePoint->[0] = $oldEdgePoint->[0] * $projWidth;
        $newEdgePoint->[0] = $newEdgePoint->[0] * $projWidth;
        $oldEdgePoint->[1] = $oldEdgePoint->[1] * -$projHeight;
        $newEdgePoint->[1] = $newEdgePoint->[1] * -$projHeight;
        push(@output, $oldEdgePoint);
        push(@output, 0);
        push(@output, $newEdgePoint);
      }
    }
    $oldX = $px;
    $oldLat = $newLat;
    $px = $px;
    $py = $py;
    push(@output, [$px, $py, $pt->[2]]);
  }
  if (scalar(grep($_ && $_->[2],@output)) > 0) {
    return(@output);
  } else { # if everything is clipped, output 'undefined' for everything
    return((0) x scalar(@input));
  }
}

sub orthDist {
  # calculates minimum distance between point (xP,yP) and line [(x1,y1) - (x2,y2)]
  my ($point1, $pointP, $point2) = @_;
  # distance formula from http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
  # with modifications for division by zero, and points outside the line
  my $x1 = $point1 -> [0];
  my $y1 = $point1 -> [1];
  my $xP = $pointP -> [0];
  my $yP = $pointP -> [1];
  my $x2 = $point2 -> [0];
  my $y2 = $point2 -> [1];
  if (($x2 == $x1) && ($y2 == $y1)) {
    # exclude case where denominator is zero
    my $p1Dist = sqrt(($x1-$xP)**2 + ($y1 - $yP)**2);
    return($p1Dist);
  }
  my $dist = abs(($x2 - $x1)*($y1 - $yP) - ($x1 - $xP)*($y2 - $y1)) /
    sqrt(($x2 - $x1)**2 + ($y2 - $y1)**2);
  if ($dist == 0) {
    # on the line, but need to consider points outside the line
    # this fixes a problem where the equator line is clipped
    my $p1Dist = sqrt(($x1-$xP)**2 + ($y1 - $yP)**2);
    my $p2Dist = sqrt(($x2-$xP)**2 + ($y2 - $yP)**2);
    my $p12Dist = sqrt(($x1-$x2)**2 + ($y1 - $y2)**2);
    my $sigma = 0.0001;
    if (($p1Dist + $p2Dist) > ($p12Dist + $sigma)) {
      # point is outside the line, use smallest distance from line end points
      $dist = ($p1Dist < $p2Dist) ? $p1Dist : $p2Dist;
    }
  }
  return($dist);
}

# Simplify curves to reduce number of points in SVG file.  This method
# should keep points if they have already been included as part of the
# simplifcation of another curve, so that shapes that overlap won't
# have any gaps. This requires a global hash of already added points.

# The linear simplifcation used here is the following algorithm:
# http://en.wikipedia.org/wiki/Ramer%E2%80%93Douglas%E2%80%93Peucker_algorithm

# Note that this algorithm is somewhat slow O(n^2), and faster (but
# less accurate) simplifcation algorithms exist. An 'accumulated
# error' algorithm could be implemented if more speed is desired.

sub simplify {
  my ($options, $pointHash, $pointRef) = @_;
  my @input = @{$pointRef};
  my $res = $options->{"resolution"};
  my @validPoints = ();
  my %pointListHash = ();
  if($res){
    # round points to specified resolution
    map {
      if($_){
        # note: ($x <=> 0) is the sign function
        $_->[0] = int($_->[0] / $res + 0.5 * ($_->[0] <=> 0)) * $res;
        $_->[1] = int($_->[1] / $res + 0.5 * ($_->[1] <=> 0)) * $res;
      }
    } @input;
  }
  foreach my $i (0..$#input) {
    if ($input[$i]) {
      push(@validPoints, $i);
    }
  }
  if (!(@validPoints) || (scalar(@validPoints) <= 2) || ($res == 0)) {
    return(@input);
  }
  $pointListHash{$validPoints[0]} = 1;
  $pointListHash{$validPoints[-1]} = 1;
  my @pointList = sort {$a <=> $b} keys(%pointListHash);
  foreach my $i (($pointList[0] + 1) .. ($pointList[-1] - 1)) {
    # add in points already added from other curves
    my $checkPoint = $input[$i];
    if ($checkPoint) {
      if ($pointHash->{$checkPoint}) {
        $pointListHash{$i} = 1;
      }
    }
  }
  @pointList = sort {$a <=> $b} keys(%pointListHash);
  my $added;
  do {
    $added = 0;
    my @firstPoints = @pointList[0..($#pointList - 1)];
    my @nextPoints = @pointList[1..$#pointList];
    foreach my $ip (0..$#firstPoints) {
      my $maxDist = 0;
      my $maxPos = 0;
      my $startPoint = $input[$firstPoints[$ip]];
      my $endPoint = $input[$nextPoints[$ip]];
      foreach my $ic (($firstPoints[$ip] + 1) .. ($nextPoints[$ip] - 1)) {
        my $checkPoint = $input[$ic];
        if ($checkPoint) {
          my $pointDist = orthDist($startPoint, $checkPoint, $endPoint);
          if ($pointDist > $maxDist) {
            $maxDist = $pointDist;
            $maxPos = $ic;
          }
        }
      }
      if (($maxDist) && ($maxDist >= $res)) {
        $pointListHash{$maxPos} = 1;
        $added++;
      }
    }
    @pointList = sort {$a <=> $b} keys(%pointListHash);
  } while ($added > 0);
  return(@input[@pointList]);
}

sub intersect {
}

sub clip {
  # clip the shape to the limits of the zoom box (if any). The
  # algorithm used here is the Sutherland-Hodgman algorithm:
  # http://en.wikipedia.org/wiki/Sutherland-Hodgman

  # Note that this algorithm can leave some zero-width polygon
  # sections in the resultant clip region.
  my ($options, $pointHash) = @_;
  my $xAdj = $options->{"xAdj"};
  my $yAdj = $options->{"yAdj"};
  my @oldPoints = @{$pointHash};
  # if all points have been previously excluded, clip them out
  if(!grep($_ && $_->[2], @oldPoints)){
    return((0) x scalar(@oldPoints));
  }
  if(!$options->{"zoomed"}){
    return(@oldPoints);
  }
  # map edge name to the array location it refers to, and the
  # direction of checking
  my %edges = ("minX" => 0, "maxX" => 1, "minY" => 2, "maxY" => 3);
  for my $edge (keys(%edges)){
    my $edgeLoc = int($edges{$edge} / 2);
    my $checkLt = $edges{$edge} % 2;
    my @processedPoints = ();
    my $lastPt = 0;
    while(@oldPoints && !$oldPoints[$#oldPoints]){
      push(@processedPoints, pop(@oldPoints));
    }
    if(@oldPoints){
      $lastPt = $oldPoints[$#oldPoints];
    }
    while(@oldPoints){
      my $pt = shift(@oldPoints);
      if (!$pt) {
        push(@processedPoints, $pt);
      } else {
        if ($checkLt xor ($pt->[$edgeLoc] >= $options->{$edge})) {
          if (!$checkLt xor ($lastPt->[$edgeLoc] >= $options->{$edge})) {
            # calculates the intersection of the line between pt and
            # lastPt and the boundary line
            my $ptA = $pt->[$edgeLoc];
            my $lptA = $lastPt->[$edgeLoc];
            my $ptB = $pt->[1 - $edgeLoc];
            my $lptB = $lastPt->[1 - $edgeLoc];
            my $edgeCoord = $options->{$edge};
            my $locProp = ($ptA - $edgeCoord) / ($ptA - $lptA);
            my $locCoord = $locProp * ($ptB - $lptB) + $ptB;
            push(@processedPoints,
                 [($edgeLoc == 0) ? $edgeCoord : $locCoord,
                  ($edgeLoc == 1) ? $edgeCoord : $locCoord, 1]);
          }
          push(@processedPoints, $pt);
        } elsif ($checkLt xor ($lastPt->[$edgeLoc] >= $options->{$edge})) {
          # same as above, not wrapped into a function
          my $ptA = $pt->[$edgeLoc];
          my $lptA = $lastPt->[$edgeLoc];
          my $ptB = $pt->[1 - $edgeLoc];
          my $lptB = $lastPt->[1 - $edgeLoc];
          my $edgeCoord = $options->{$edge};
          my $locProp = ($ptA - $edgeCoord) / ($ptA - $lptA);
          my $locCoord = $locProp * ($ptB - $lptB) + $ptB;
          push(@processedPoints,
               [($edgeLoc == 0) ? $edgeCoord : $locCoord,
                ($edgeLoc == 1) ? $edgeCoord : $locCoord, 1]);
        }
        $lastPt = $pt;
      }
    } # ends while
    push(@oldPoints, @processedPoints);
  } # ends for
  return(@oldPoints);
}

sub boundBox {
  # Determines post-projection rectangular bounding box for a polygon
  # (e.g. for determining if borders are useful, or zooming in to a
  # particular region)
  my @input = @_;
  @input = grep($_, @input);
  my $minX = 0;
  my $minY = 0;
  my $maxX = 0;
  my $maxY = 0;
  if (@input) {
    my $minPoint = $input[0];
    $minX = $minPoint -> [0];
    $minY = $minPoint -> [1];
    $maxX = $minPoint -> [0];
    $maxY = $minPoint -> [1];
    foreach my $point (@input) {
      my $px = $point -> [0];
      my $py = $point -> [1];
      $minX = $px if ($px < $minX);
      $minY = $py if ($py < $minY);
      $maxX = $px if ($px > $maxX);
      $maxY = $py if ($py > $maxY);
    }
  }
  return(($minX, $minY, $maxX, $maxY));
}

sub shapeNames {
  # retrieve shape name data from a shape database
  # for countries, three-letter ISO codes are preferred
  my ($shapeData) = @_;
  my $outCountry = "";
  my $outSov = ""; # 'sovereignty' is easy to mis-spell, and 'sov' is shorter
  my $outRegion = "";
  if (exists($shapeData->{"Name2"})) {
    # Name2 seems to be 3-letter code for Point data
    $outSov = $shapeData->{"Name2"};
    $outCountry = $shapeData->{"Name2"};
  } elsif (exists($shapeData->{"sov_a3"}) && exists($shapeData->{"adm0_a3"})) {
    # 3-letter code for polygon data
    $outSov = $shapeData->{"sov_a3"};
    $outCountry = $shapeData->{"adm0_a3"};
  } elsif (exists($shapeData->{"SOV_A3"}) && exists($shapeData->{"ADM0_A3"})) {
    # 3-letter code for polygon data
    $outSov = $shapeData->{"SOV_A3"};
    $outCountry = $shapeData->{"ADM0_A3"};
  } elsif (exists($shapeData->{"SOV"}) && exists($shapeData->{"COUNTRY"})) {
    # 10m data uses full name, rather than 3-letter code
    # should a lookup table be used for conversion?
    $outSov = $shapeData->{"SOV"};
    $outCountry = $shapeData->{"COUNTRY"};
  } elsif (exists($shapeData->{"ISO"})) {
    $outSov = $shapeData->{"ISO"}; # no 10m Sov data in state/province file
    $outCountry = $shapeData->{"ISO"};
  } elsif (exists($shapeData->{"NAME_0"})) {
    $outSov = $shapeData->{"NAME_0"}; # this is used in GADM data
    $outCountry = $shapeData->{"NAME_0"};
  } elsif (exists($shapeData->{"featurecla"})) {
    if (($shapeData->{"featurecla"} eq "River") && exists($shapeData->{"name"})) {
      $outRegion = $shapeData->{"name"};
    }
  } elsif (exists($shapeData->{"FEATURECLA"})) {
    if (($shapeData->{"FEATURECLA"} eq "River") && exists($shapeData->{"NAME"})) {
      $outRegion = $shapeData->{"NAME"};
    }
  }
  if (exists($shapeData->{"name_1"})) {
    $outRegion = $shapeData->{"name_1"};
  } elsif (exists($shapeData->{"NAME_1"})) {
    $outRegion = $shapeData->{"NAME_1"};
  }
  utf8::encode($outRegion);
  utf8::encode($outSov);
  utf8::encode($outCountry);
  return($outSov, $outCountry, $outRegion);
}

sub chopPoly {
  # closes up lines in a polygon that contain discontinuous points
  # this usually happens when part of a polygon goes outside the
  # projection area
  # input is a sequence of Geo::Points, output is a string that can be
  # used as a path description in an SVG file
  my ($options, $pointHash, $joinEnds) = @_;
  my @oldPoints = @{$pointHash};
  if (!@oldPoints) {
    return (""); # nothing to process
  }
  my $dp = $options->{"roundDP"};
  my $xAdj = $options->{"xAdj"};
  my $yAdj = $options->{"yAdj"};
  # trim off start and end undefined / missing points
  while(@oldPoints && !$oldPoints[0]){
    shift(@oldPoints);
  }
  while(@oldPoints && !$oldPoints[$#oldPoints]){
    pop(@oldPoints);
  }
  # create list of contiguous sequences
  my @contigLists = ();
  my $currentList = [];
  push(@contigLists, $currentList);
  foreach(@oldPoints){
    if($_){
      push(@{$currentList}, $_);
    } else {
      if(@{$currentList}){
        $currentList = [];
        push(@contigLists, $currentList);
      }
    }
  }
  # If a discontinuity exists, check if start/end point are the
  # same. If so, rotate the array to start/end on the first
  # discontinuity. This fixes an issue where polygons that cross the
  # boundaries of a full-sized location / world projection are
  # inappropriately joined
  if(scalar(@contigLists) > 1){
    my $distX = $oldPoints[0]->[0] - $oldPoints[$#oldPoints]->[0];
    my $distY = $oldPoints[0]->[1] - $oldPoints[$#oldPoints]->[1];
    my $distFL = sqrt($distX*$distX + $distY*$distY);
    my $minDist = 0.01;
    if($distFL < $minDist){
      # tack start list on the end of the end list, remove start list
      push(@{$contigLists[$#contigLists]}, @{shift @contigLists});
    }
  }
  if (!(@oldPoints) || !(@{$contigLists[0]})) {
    return ("");
  }
  my $printfString = ($dp > -1) ? ("%.".$dp."g,%.".$dp."g") : ("%g,%g");

  my @subPaths = ();
  foreach (@contigLists) {
    my @currentList = @{$_};
    if(scalar(@currentList) > 1){ # ignore single-point subpaths
      my @pointStrs = map {
        sprintf($printfString, ($_->[0] + $xAdj), ($_->[1] + $yAdj));
      } @currentList;
      push(@subPaths, "M ".join(" ",@pointStrs).($joinEnds?" Z":""));
    }
  }
  my $pointPath = join(" ", @subPaths);

  return($pointPath);
}

my $projOpts = {
    "resolution" => 0,
    "roundDP" => -1,
    "projection" => "", # was sourceProj, targetProj
    "centreLn" => "",  # projection centre X / longitude
    "centreLt" => "",  # projection centre Y / latitude
    "svgWidth" => "", # width of SVG file, excluding padding
    "svgHeight" => "", # height of SVG file, excluding padding
    "xAdj" => 0, # X adjustment to fit in SVG bounds
    "yAdj" => 0,  # Y adjustment to fit in SVG bounds
    "xScale" => "",  # Scale factor of points (for zoom box)
    "yScale" => "",  # Scale factor of points (for zoom box)
    "padding" => "" # amount of padding to add (for zoom box)
};

my $mapType = "";

my $colourTheme = "";

my $centreCountry = "";

my $landColour = ""; # outside area (or default land colour)
my $seaColour = ""; # ocean / sea / lake
my $borderColour = "";

# Extra colours (only used if land/sea/border not specified)

my $subLandColour = ""; # subject colour
my $othLandColour = ""; # other lands of the same political unity
my $coastColour = ""; # lake's coast, rivers, sea coasts
my $polBordColour = ""; # Other minor political borders
my $intBordColour = ""; # border colour for areas of interest

my $printLines = "";

my $printKey = 1; # true

my $latSep = 30; # separation (degrees) for lines of latitude
my $longSep = 30; # separation (degrees) for lines of longitude

my $latLimit = ""; # maximum latitude to display

my @dataFiles = (); # files containing numerical data

my %subjectNames = ();
my %politicalNames = ();
my %onlyNames = ();
my %zoomNames = ();

my %zoomExtents = ();

my $debugLevel = 0;

my @shpFileBases = ();

my @commandLine = @ARGV;

# Return codes
# 0: success
# 1: unknown command-line option
# 2: file error
# 3: projection error

# extract command line arguments
while (@ARGV) {
  my $argument = shift @ARGV;
  if (((-f $argument) && ($argument =~ /\.shp$/)) ||
      (-f ($argument.".shp"))) { # file existence check
    if ($argument =~ /\.shp$/) {
      $argument =~ s/\.shp$//;
    }
    printf(STDERR "Adding '%s.shp' to list of input files\n", $argument)
      if ($debugLevel > 0);
    push (@shpFileBases, $argument);
  } elsif (-f $argument) {
    print(STDERR "Error: Invalid file extension for '$argument'. ".
          "Please use files with '.shp' extension\n");
    usage();
    exit(2);
  } else {
    if ($argument eq "-help") {
      # flag options [no arguments]
      usage();
      exit(0);
    } elsif ($argument eq "-v") {
      $debugLevel = 1;
      print(STDERR "Enabling verbose output\n");
    } elsif ($argument eq "-vv") {
      $debugLevel = 2;
      print(STDERR "Enabling very verbose output\n");
    } elsif ($argument eq "-lines") {
      $printLines = 3;          # true
      print(STDERR "Printing lines of latitude and longitude\n");
    } elsif ($argument eq "-nolines") {
      $printLines = 0;          # false
      print(STDERR "Printing lines of latitude and longitude\n");
    } elsif ($argument eq "-nokey") {
      $printKey = 0;            # false
      print(STDERR "Will not display a key for the heatmap\n");
    } elsif ($argument eq "-psize") {
      # floating point options
      $projOpts->{"pointSize"} = shift @ARGV;
      printf(STDERR "Point size changed to %s\n",
             $projOpts->{"pointSize"});
    } elsif ($argument eq "-round") {
      $projOpts->{"roundDP"} = shift @ARGV;
      printf(STDERR "will round to %s decimal places\n",
             $projOpts->{"roundDP"});
    } elsif ($argument eq "-res") {
      $projOpts->{"resolution"} = shift @ARGV;
      printf(STDERR "Resolution changed to %s\n",
             $projOpts->{"resolution"});
    } elsif ($argument eq "-type") {
      # string options
      my $newArg = shift(@ARGV);
      if ($newArg =~ /(location|locator|area|world|orthographic)/) {
        $mapType = $newArg;
      } else {
        print(STDERR "Error: '$newArg' is not a valid map type. ".
              "Only accepts '(location|locator|area|world|orthographic)'\n");
        usage();
        exit(1);
      }
      print(STDERR "Map type changed to '$mapType'\n");
    } elsif ($argument =~ /^-colo(u)?r$/) {
      $colourTheme = shift @ARGV;
      print(STDERR "Colour theme changed to '$colourTheme'\n");
    } elsif (($argument eq "-centre") || ($argument eq "-center")) {
      my $newArg = shift(@ARGV);
      if ($newArg =~ /,/) {
        my @centre = split(/,/,$newArg,2);
        $projOpts->{"centreLn"} = $centre[0];
        $projOpts->{"centreLt"} = $centre[1];
        printf(STDERR "Map centre changed to '%s,%s'\n",
               $projOpts->{"centreLn"},
               $projOpts->{"centreLt"});
      } else {
        $centreCountry = $newArg;
        print(STDERR "Map centre changed to centre of country ".
              "'$newArg' (if found)\n");
      }
    } elsif ($argument =~ /^-sub(ject)?/) {
      my $newArg = shift @ARGV;
      $subjectNames{$newArg} = 1;
      print(STDERR "Adding '$newArg' as a map subject\n");
    } elsif ($argument eq "-pol") {
      my $newArg = shift @ARGV;
      $politicalNames{$newArg} = 1;
      print(STDERR "Adding '$newArg' as another politically-related area\n");
    } elsif ($argument eq "-data") {
      my $newArg = shift @ARGV;
      push(@dataFiles, $newArg);
      print(STDERR "Will extract numerical data from ".
            "'$newArg' (assuming <ISO 3A-code>,<float>)\n");
    } elsif ($argument eq "-only") {
      my $newArg = shift @ARGV;
      $onlyNames{$newArg} = 1;
    } elsif ($argument eq "-zoom") {
      my $newArg = shift @ARGV;
      if ($newArg =~ /([0-9\.]+),([0-9\.]+),([0-9\.]+),([0-9\.]+)/) {
        $projOpts->{"zoomed"} = 1;
        $projOpts->{"manualZoom"} = 1;
        $projOpts->{"minX"} = $1;
        $projOpts->{"minY"} = $2;
        $projOpts->{"maxX"} = $3;
        $projOpts->{"maxY"} = $4;
        printf(STDERR "Setting initial zoom box to (%s,%s)-(%s,%s)\n",
               $projOpts->{"minX"},$projOpts->{"minY"},
               $projOpts->{"maxX"},$projOpts->{"maxY"});
      } else {
        $zoomNames{$newArg} = 1;
        $projOpts->{"zoomed"} = 1;
        print(STDERR "Zooming map to include '$newArg'\n");
      }
    } elsif ($argument eq "-landcol") {
      $landColour = shift @ARGV;
      print(STDERR "Land colour changed to '$landColour'\n");
    } elsif ($argument eq "-seacol") {
      $seaColour = shift @ARGV;
      print(STDERR "Sea colour changed to '$seaColour'\n");
    } elsif ($argument eq "-bordcol") {
      $borderColour = shift @ARGV;
      print(STDERR "Border colour changed to '$borderColour'\n");
    } elsif ($argument eq "-proj") {
      $projOpts->{"projection"} = shift @ARGV;
      printf(STDERR "Projection changed to '%s'\n",
             $projOpts->{"projection"});
    } else {
      print(STDERR "Error: Unknown command-line option or non-existent file, ".
            "'$argument'\n");
      usage();
      exit(1);
    }
  }
}

if(scalar(@shpFileBases) == 0){
  print(STDERR "Error: No input files specified\n");
  usage();
  exit(2);
}

if(keys(%onlyNames)){
  print(STDERR "Will only display the following features on the map:\n");
  print(STDERR "  ".join("\n  ",keys(%onlyNames))."\n");
}


$mapType = "location" unless $mapType;

if ($mapType =~ /^(location|locator|area)/) {
  $projOpts->{"projection"} = "equirectangular"
    unless $projOpts->{"projection"};
  $projOpts->{"svgWidth"} = 1100 unless $projOpts->{"svgWidth"};
  $projOpts->{"svgHeight"} = 550 unless $projOpts->{"svgHeight"};
  $printLines = 0 if ($printLines eq "");
} elsif ($mapType =~ /^world/) {
  # Winkel Tripel projection for world maps
  $projOpts->{"projection"} = "wintri"
    unless $projOpts->{"projection"};
  $projOpts->{"svgWidth"} =
    ((4 + 2 * pi) / (2 * pi) * 550) unless $projOpts->{"svgWidth"};
  $projOpts->{"svgHeight"} = 550 unless $projOpts->{"svgHeight"};
  $printLines = 0 if ($printLines eq "");
} elsif ($mapType =~ /^ortho(graphic)?/) {
  # Orthographic projection for locator maps
  $projOpts->{"projection"} = "orthographic"
    unless $projOpts->{"projection"};
  $projOpts->{"xScale"} = 1 unless $projOpts->{"xScale"};
  $projOpts->{"yScale"} = 1 unless $projOpts->{"yScale"};
  $projOpts->{"svgWidth"} = 550 unless $projOpts->{"svgWidth"};
  $projOpts->{"svgHeight"} = 550 unless $projOpts->{"svgHeight"};
  $printLines = 3 if ($printLines eq "");
}

# specify width / height if not modified by type definitions
$projOpts->{"pointSize"} = 0.5 unless $projOpts->{"pointSize"};
$projOpts->{"svgWidth"} = 1100 unless $projOpts->{"svgWidth"};
$projOpts->{"svgHeight"} = 1100 unless $projOpts->{"svgHeight"};
$projOpts->{"xScale"} = 1 unless $projOpts->{"xScale"};
$projOpts->{"yScale"} = 1 unless $projOpts->{"yScale"};
# adjustment allows for 1.5px border around SVG image
$projOpts->{"xAdj"} = $projOpts->{"svgWidth"} / 2 + 1.5;
$projOpts->{"yAdj"} = $projOpts->{"svgHeight"} / 2 + 1.5;
$projOpts->{"padding"} = 0 unless $projOpts->{"padding"};

$colourTheme = $mapType unless $colourTheme;

# protect user-specified land/sea/border colours

my $landColTemp = $landColour;
my $seaColTemp = $seaColour;
my $bordColTemp = $borderColour;

# default colours (a consensus of consensus colours, if you like)

$landColour = "#FEFEE9";
$borderColour = "#646464";
$seaColour = "#C6ECFF";

# set up colours for particular map types.
# location: used as backgrounds for automatic geo-localisation
# locator: display an article's subject area of occupation
# area: focus on group their respective area of control
# world: grey world map
# orthographic: grey/green orthographic map\n");

if ($colourTheme eq "location") {
  $subLandColour = "#FEFEE9";
  $othLandColour = "#F6E1B9";
  $landColour = "#E0E0E0";
  $borderColour = "#646464";
  $seaColour = "#C6ECFF";
  $coastColour = "#0978AB";
} elsif (($colourTheme eq "locator") || ($colourTheme eq "area")) {
  $subLandColour = "#F07568";
  $othLandColour = "#FFFFD0";
  $landColour = "#F7D3AA";
  $polBordColour = "#D0C0A0";
  $intBordColour = "#E0584E";
  $borderColour = "#A08070";
  # assumed from locator map definition
  $seaColour = "#9EC7F3";
  $coastColour = "#1821DE";
} elsif ($colourTheme eq "world") {
  $landColour = "#E0E0E0";
  $coastColour = "#646464";
  $seaColour = "#FFFFFF";
  $borderColour = 'none';
} elsif ($colourTheme eq "orthographic") {
  $subLandColour = "#336733";
  $othLandColour = "#C6DEBD";
  $landColour = "#E0E0E0";
  $seaColour = "#FFFFFF";
  $intBordColour = "#335033";
  $borderColour = "#646464";
  $coastColour = "#646464";
  #$borderColour = "#FFFFFF";
}

# replace subtype colours with colours of super if not explicitly specified

$subLandColour = $landColour if !$subLandColour;
$othLandColour = $landColour if !$othLandColour;

$coastColour = $seaColour if !$coastColour;

$intBordColour = $borderColour if !$intBordColour;
$polBordColour = $borderColour if !$polBordColour;

# overrides when colours are specified as arguments

$landColour = $landColTemp if $landColTemp;
$subLandColour = $landColTemp if $landColTemp;
$othLandColour = $landColTemp if $landColTemp;

$seaColour = $seaColTemp if $seaColTemp;
$coastColour = $seaColTemp if $seaColTemp;

$borderColour = $bordColTemp if $bordColTemp;
$intBordColour = $bordColTemp if $bordColTemp;
$polBordColour = $bordColTemp if $bordColTemp;

if(!keys(%subjectNames) && !keys(%politicalNames)){
  # no subjects or other political areas, so treat entire world as subject
  $landColour = $subLandColour;
  $borderColour = $intBordColour;
}


my %addedPoints = ();

# centre of map is a feature, so hunt for it
if ($centreCountry) {
  my ($minX, $minY, $maxX, $maxY);
   foreach my $shpFileBase (@shpFileBases) {
    print(STDERR "Loading $shpFileBase to hunt for features... ")
      if ($debugLevel> 0);
    my $shp = new Geo::ShapeFile($shpFileBase)
      or die("Unable to load $shpFileBase.shp");
    my $shapeCount = $shp -> shapes();
    # iterate through shapes in shapefile
    for my $shapeNum (1 .. ($shapeCount)) { # 1-based indexing
      my %data = $shp->get_dbf_record($shapeNum);
      my ($sovName, $countryName, $regionName) = shapeNames(\%data);
      if (($countryName eq $centreCountry) || ($sovName eq $centreCountry)) {
        # found it, so find polygon extents in degrees
        my $shape = $shp->get_shp_record($shapeNum);
        my @shapePoints = $shape->points();
        foreach my $point (@shapePoints) {
          my $px = $point->X;
          my $py = $point->Y;
          if(!defined($minX)){
            ($minX, $minY, $maxX, $maxY) = ($px, $py, $px, $py);
          }
          if(($px - $minX) > 270){ # deal with pesky -180/180 wrapping issues
            $px -= 360;
          }
          if(($maxX - $px) > 270){
            $px += 360;
          }
          $minX = $px if($px < $minX);
          $minY = $py if($py < $minY);
          $maxX = $px if($px > $maxX);
          $maxY = $py if($py > $maxY);
        }
        my $midLong = ($minX + $maxX) / 2;
        if($midLong > 180){
          $midLong -= 360;
        }
        if($midLong < -180){
          $midLong += 360;
        }
        my $midLat = ($minY + $maxY) / 2;
        $projOpts->{"centreLn"} = $midLong;
        $projOpts->{"centreLt"} = $midLat;
        printf(STDERR "Found centre feature for '%s' at point (%s,%s)... ",
               $countryName, $midLong, $midLat) if ($debugLevel> 0);
      }
    }
    print(STDERR "done.\n") if ($debugLevel> 0);
  }
}

if ((!$projOpts->{"projection"}) ||
    ($projOpts->{"projection"} !~
     /^(equirectangular|wintri|orthographic)$/)) {
  printf(STDERR "Error: Unknown projection '%s'\n".
         "Should be one of: equirectangular, wintri, orthographic\n",
         $projOpts->{"projection"});
  usage();
  exit(3);
}

if (($projOpts->{"centreLn"} eq "") ||
    ($projOpts->{"centreLt"} eq "")) {
  $projOpts->{"centreLn"} = 0;
  $projOpts->{"centreLt"} = 0;
}

if ($projOpts->{"manualZoom"}) {
  # user specified zoom box, so convert to non-adjusted values
  # TODO: make sure this is the right thing to do
  $projOpts->{"minX"} -= $projOpts->{"xAdj"};
  $projOpts->{"minY"} -= $projOpts->{"yAdj"};
  $projOpts->{"maxX"} -= $projOpts->{"xAdj"};
  $projOpts->{"maxY"} -= $projOpts->{"yAdj"};
}

# Determine projection extents, if necessary
if (keys(%zoomNames)) {
  foreach my $shpFileBase (@shpFileBases) {
    if (-f ($shpFileBase.".prj")) {
      open(PROJFILE, "< $shpFileBase".".prj")
        or die("Unable to load $shpFileBase.prj");
      while (<PROJFILE>) {
        chomp;
        $projOpts->{"sourceProj"} = $_;
      }
    }
    print(STDERR "Loading $shpFileBase to hunt for features... ")
      if ($debugLevel> 0);
    my $shp = new Geo::ShapeFile($shpFileBase)
      or die("Unable to load $shpFileBase.shp");
    my $shapeCount = $shp -> shapes();
    # iterate through shapes in shapefile
    for my $shapeNum (1 .. ($shapeCount)) { # 1-based indexing
      my %data = $shp->get_dbf_record($shapeNum);
      my ($sovName, $countryName, $regionName) = shapeNames(\%data);
      if (($zoomNames{$countryName}) || ($zoomNames{$sovName})) {
        my $shape = $shp->get_shp_record($shapeNum);
        my @tmpPoints = $shape->points();
        my $pointCount = scalar(@tmpPoints);
        @tmpPoints = project($projOpts, \@tmpPoints);
        my ($tminX, $tminY, $tmaxX, $tmaxY) = boundBox(@tmpPoints);
        if($tminX || $tminY || $tmaxX || $tmaxY){ # make sure there are points
          $projOpts->{"minX"} = $tminX
            if ((!exists($projOpts->{"minX"})) ||
                ($projOpts->{"minX"} > $tminX));
          $projOpts->{"minY"} = $tminY
            if ((!exists($projOpts->{"minY"})) ||
                ($projOpts->{"minY"} > $tminY));
          $projOpts->{"maxX"} = $tmaxX
            if ((!exists($projOpts->{"maxX"}))
                || ($projOpts->{"maxX"} < $tmaxX));
          $projOpts->{"maxY"} = $tmaxY
            if ((!exists($projOpts->{"maxY"}))
                || ($projOpts->{"maxY"} < $tmaxY));
          printf(STDERR "Found extent for '$countryName ($regionName)', ".
                 "adjusting zoom box to (%s,%s)-(%s,%s)...",
                 $projOpts->{"minX"} + $projOpts->{"xAdj"},
                 $projOpts->{"minY"} + $projOpts->{"yAdj"},
                 $projOpts->{"maxX"} + $projOpts->{"xAdj"},
                 $projOpts->{"maxY"} + $projOpts->{"yAdj"}) if ($debugLevel> 0);
        }
      }
    }
    print(STDERR "done.\n") if ($debugLevel> 0);
  }
}

$projOpts->{"pointSize"} = 1 unless $projOpts->{"pointSize"};
$printLines = 0 if ($printLines eq "");

# adjust SVG dimensions to fit in with zoom extents
if ($projOpts->{"zoomed"}) {
  printf(STDERR "old SVG width: %0.2f\n", $projOpts->{"svgWidth"});
  printf(STDERR "old SVG height: %0.2f\n", $projOpts->{"svgHeight"});
  printf(STDERR "Zooming to projected region (%0.2f,%0.2f)-(%0.2f,%0.2f)\n",
         $projOpts->{"minX"} + $projOpts->{"xAdj"},
         $projOpts->{"minY"} + $projOpts->{"yAdj"},
         $projOpts->{"maxX"} + $projOpts->{"xAdj"},
         $projOpts->{"maxY"} + $projOpts->{"yAdj"})
     if ($debugLevel > 0);
  my $width = $projOpts->{"maxX"} - $projOpts->{"minX"};
  my $height = $projOpts->{"maxY"} - $projOpts->{"minY"};
  my $ratio = $width / $height;
  # requre width to be over 2x height before changing reference dimension
  my $refDim = ($width > ($height)) ? $width : $height;
  my $oldW = $projOpts->{"svgWidth"};
  my $oldH = $projOpts->{"svgHeight"};
  $projOpts->{"svgWidth"} = ($refDim == $width) ? 1100 : 550 * $ratio;
  $projOpts->{"svgHeight"} = ($refDim == $height) ? 550 : 1100 / $ratio;
  $projOpts->{"xScale"} = $oldW / $width;
  $projOpts->{"yScale"} = $oldH / $height;
  my $relMagX = ($projOpts->{"svgWidth"} * $projOpts->{"xScale"}) / $oldW;
  my $relMagY = ($projOpts->{"svgHeight"} * $projOpts->{"yScale"}) / $oldH;
  printf(STDERR "[Relative SVG magnification is (%0.3g,%0.3g)x]\n",
         $relMagX, $relMagY) if ($debugLevel > 0);
  printf(STDERR "[SVG scale is (%0.3g,%0.3g)x]\n",
         $projOpts->{"xScale"}, $projOpts->{"yScale"}) if ($debugLevel > 0);
  # modify x/y adjust for new map boundaries plus a bit of padding
  $projOpts->{"padding"} = $projOpts->{"pointSize"} * 11;
  $projOpts->{"xAdj"} = $projOpts->{"minX"} * -$relMagX
    + 1.5 + $projOpts->{"padding"};
  $projOpts->{"yAdj"} = $projOpts->{"minY"} * -$relMagY
    + 1.5 + $projOpts->{"padding"};
}

$projOpts->{"minX"} = 1.5 - $projOpts->{"xAdj"};
$projOpts->{"minY"} = 1.5 - $projOpts->{"yAdj"};
$projOpts->{"maxX"} = $projOpts->{"minX"} + $projOpts->{"svgWidth"}
  + $projOpts->{"padding"} * 2;
$projOpts->{"maxY"} = $projOpts->{"minY"} + $projOpts->{"svgHeight"}
  + $projOpts->{"padding"} * 2;

my $svg = SVG->new('height' => $projOpts->{"svgHeight"} + $projOpts->{"padding"} * 2 + 3,
                   'width' => $projOpts->{"svgWidth"} + $projOpts->{"padding"} * 2 + 3);
$svg->comment("Created using David Eccles' (gringer) perlshaper.pl script, ".
              "version $progVersion, ".
              "http://en.wikipedia.org/wiki/User:Gringer/perlshaper");
$svg->comment("Command line: ".join(" ",@commandLine));

# extract / colour based on data from files
my $countryColours = "";
my %countryValues = ();
my $minCValue = 100;
my $maxCValue = 0;
foreach my $dataFile (@dataFiles) {
  $countryColours .= sprintf("\n/* Using data from '%s' */",
                             $dataFile);
  if (!(-f $dataFile)) {
    printf(STDERR "Warning: data file, '%s', was not found",
           $dataFile);
    next;
  }
  open(INFILE, "< $dataFile") or
    die("Unable to open $dataFile");
  printf(STDERR "Extracting data from '%s'...", $dataFile) if ($debugLevel> 0);
  my $addedCountries = 0;
  while (<INFILE>) {
    if (!(/^#/)) {
      if (/^([a-zA-Z]{3}),([0-9\.\-]*)$/) {
        my $country = uc($1);
        my $value = $2;
        if ($value < $minCValue) {
          $minCValue = $value;
        }
        if ($value > $maxCValue) {
          $maxCValue = $value;
        }
        $countryValues{$country} = $value;
        $addedCountries++;
      }
    }
  }                             # closes while(<INFILE>)
  printf(STDERR " done (data for %d countries found)!\n",
         $addedCountries) if ($debugLevel> 0);
}
# can now determine proper scaling as the range is known
# [$minCValue .. $maxCValue]
foreach my $country (keys(%countryValues)) {
  # gradient is cyan->yellow to avoid colour confusion lines
  # (see http://www.colblindor.com/2009/01/19/colorblind-colors-of-confusion/)
  # (i.e. from #00FFFF to #FFFF00)
  my $value = $countryValues{$country};
  # convert to 0.0-255.0 and round to nearest integer
  my $redAmt =
    sprintf("%02X", (($value - $minCValue) * 255) /
            ($maxCValue - $minCValue) + 0.5);
  my $blueAmt =
    sprintf("%02X", (($maxCValue - $value) * 255) /
            ($maxCValue - $minCValue) + 0.5);
  my $colour = sprintf("#%sFF%s", $redAmt, $blueAmt);
  $countryColours .= sprintf("\n.%s{fill: %s;} /* %s */",
                             $country, $colour, $value);
}
if ($countryColours) {
  $countryColours .= "\n".(qq|
rect.heatmap{
  fill: url(#grHeatmap);
  stroke: black;
  stroke-width: 1px;
}

text.heatmap{
  fill: black;
  stroke: none;
  font-size: 20px;
}
|);
}

# style sheet
my $styleSheet = $svg->style('id' => 'mapStyle', 'type' => 'text/css');
my $highlightWidth = $projOpts->{"pointSize"} * 2.5;
my $riverWidth = $projOpts->{"pointSize"} / 2;
$styleSheet->CDATA(qq|
/* Cascading Style Sheet (CSS) definitions for region colours */

/* land: $landColour -- outside area (or default land colour)
 * sea: $seaColour -- ocean / sea / lake
 * border: $borderColour -- land borders
 * subject land: $subLandColour -- subject colour
 * other land: $othLandColour -- other lands of the same political unity
 * coast: $coastColour -- lake's coast, rivers, sea coasts
 * political border: $polBordColour -- Other minor political borders
 * interest border: $intBordColour -- border colour for areas of interest
 */

/* basic styles */
path{
   stroke-linejoin: bevel;
}

path,circle{
   stroke-width: $projOpts->{"pointSize"};
}

/* Main region definition */

.region{
   fill: $landColour;
   stroke: $borderColour;
   stroke-width: $projOpts->{"pointSize"};
}

/* Political groupings */

.political{
   fill: $othLandColour;
   stroke: $polBordColour;
}

/* Subject / area of interest */

.subject{
   fill: $subLandColour;
   stroke: $intBordColour;
}

circle.highlight{
   stroke-width: $highlightWidth;
}

/* Rivers */

.river{
   fill: none;
   stroke: $coastColour;
   stroke-width: $riverWidth;
}

/* Sea and decoration */
.seabase{
   fill: $seaColour;
   stroke: none;
}

.mapborder{
   fill: none;
   stroke: $coastColour;
   stroke-width: 1.25;
}

.mapshadow{
  fill: url(#grSea);
  stroke: none;
}

/* Grid lines */
.grid{
  fill: none;
  stroke: $coastColour;
}

/* Colouring a particular country, use '*' with group ID */
/*
 * #gNZL *{
 *   fill: #ff0000;
 *   stroke: black;
 * }
 */
/* ...or '.' with country name */
/*
 * .NZL *{
 *   fill: #ff0000;
 *   stroke: black;
 * }
 */
|.$countryColours);

# place sea layer beneath image

my $seaPath = "";

if ($mapType =~ /^ortho/ && !$projOpts->{"zoomed"}) {
  my $svgDefs = $svg->defs('id' => 'defSea');
  my $seaGradient = $svgDefs->gradient('id' => "grSea", '-type' => "radial",
                                       'fx' => 0.625, 'fy' => 0.375,
                                       'class' => 'seabase gradient');
  $seaGradient->stop('id' => 'sEdge', offset => 0,
                     'style' => 'stop-color:#FFFFFF;stop-opacity:0;');
  $seaGradient->stop('id' => 'sMid', offset => (1/sqrt(2)),
                     'style' => 'stop-color:#808080;stop-opacity:0.0625;');
  $seaGradient->stop('id' => 'sCentre', offset => 1,
                     'style' => 'stop-color:#000000;stop-opacity:0.25;');
  my $seaGroup = $svg->group('id' => "gSeaBase");
  $seaGroup->circle(id => "cSeaBase", 'cx' => 0 + $projOpts->{"xAdj"},
                    'cy' => 0 + $projOpts->{"yAdj"},
                    'r' => ($projOpts->{"svgWidth"} / 2),
                    'class' => 'seabase');
}

# place sea -- a line around the edge of the world for unzoomed maps

if (($mapType =~ /^(location|locator|area|world)/) &&
    !$projOpts->{"zoomed"}) {
  print(STDERR "Finding world edge...") if ($debugLevel> 0);
  my $seaGroup = $svg->group('id' => "gSeaBase");
  my $minLong = 0;
  my $minLongVal;
  my $maxLong = 0;
  my $maxLongVal;
  my $maxOffset = 0;
  for (my $longPoint = 180; $longPoint <= 540; $longPoint += 0.1) {
    my $point = transform($projOpts, [$longPoint - 360, 0])
      or die("cannot transform");
    if ($point->[2]) {
      if (!defined($minLongVal) || (($point -> [0]) < $minLongVal)) {
        $minLongVal = $point->[0];
        $minLong = sprintf("%.4f", $longPoint)+0;
      }
      if (!defined($maxLongVal) || (($point -> [0]) > $maxLongVal)) {
        $maxLongVal = $point->[0];
        $maxLong = sprintf("%.4f", $longPoint)+0;
      }
    }
  }
  if($minLong > $maxLong){
    $maxLong += 360;
  }
  my @linePoints = ();
  for (my $latPoint = -90; $latPoint <= 90; $latPoint += 0.1) {
    push(@linePoints, [$minLong - 360, $latPoint]);
  }
  for (my $longPoint = $minLong; $longPoint < $maxLong; $longPoint += 0.1) {
    push(@linePoints, [$longPoint - 360, 90]);
  }
  for (my $latPoint = 90; $latPoint >= -90; $latPoint -= 0.1) {
    push(@linePoints, [$maxLong - 360, $latPoint]);
  }
  for (my $longPoint = $maxLong; $longPoint > $minLong; $longPoint -= 0.1) {
    push(@linePoints, [$longPoint - 360, -90]);
  }
  while($minLong > 180){
    $minLong -= 360;
  }
  while($maxLong > 180){
    $maxLong -= 360;
  }
  $projOpts->{"minLong"} = $minLong;
  $projOpts->{"maxLong"} = $maxLong;
  @linePoints = project($projOpts, \@linePoints);
  @linePoints = simplify($projOpts, \%addedPoints, \@linePoints);
  $seaPath = chopPoly($projOpts, \@linePoints, 1);
  $seaGroup->path('id' => 'pSeaBase', 'd' => $seaPath, 'class' => 'seabase');
  printf(STDERR " done (%s,%s)!\n", $minLong, $maxLong) if ($debugLevel> 0);
}

# place sea -- a rectangle around the image for zoomed maps
if ($projOpts->{"zoomed"}) {
  my $seaGroup = $svg->group('id' => "gSeaBase");
  my @linePoints = ([$projOpts->{"minX"},$projOpts->{"minY"},1],
                    [$projOpts->{"maxX"},$projOpts->{"minY"},1],
                    [$projOpts->{"maxX"},$projOpts->{"maxY"},1],
                    [$projOpts->{"minX"},$projOpts->{"maxY"},1]);
  $seaPath = chopPoly($projOpts, \@linePoints, 1);
  $seaGroup->path('id' => 'pSeaBase', 'd' => $seaPath, 'class' => 'seabase');
}


my $worldGroup = $svg->group('id' => "gCountries", 'class' => 'countries');

my $fileNum = 0;
foreach my $shpFileBase (@shpFileBases) {
  $worldGroup->comment("Using data from ${shpFileBase}.shp");
  if (-f ($shpFileBase.".prj")) {
    open(PROJFILE, "< $shpFileBase".".prj")
      or die("Unable to load $shpFileBase.prj");
    while (<PROJFILE>) {
      chomp;
      $projOpts->{"sourceProj"} = $_;
    }
  }
  $fileNum++;
  # The ShapeFile 'new' procedure doesn't seem to have good error
  # checking (i.e. it will happily load a random file of a
  # sufficiently large size), so 'die' may not be useful here.
  print(STDERR "Loading $shpFileBase...") if ($debugLevel> 0);
  my $shp = new Geo::ShapeFile($shpFileBase)
    or die("Unable to load $shpFileBase.shp");
  my $shapeCount = $shp -> shapes();
  print(STDERR " done ($shapeCount shapes loaded).\n") if ($debugLevel> 0);
  # iterate through shapes in shapefile
  print(STDERR "Reprojecting and simplifying shapes...") if ($debugLevel> 0);
  for my $shapeNum (1 .. ($shapeCount)) { # 1-based indexing
    my $shape = $shp->get_shp_record($shapeNum);
    my %data = $shp->get_dbf_record($shapeNum);
    my $type = $shp->type($shape->shape_type());
    my $partCount = $shape->num_parts();
    my @tmpPoints = $shape->points();
    my $pointCount = scalar(@tmpPoints);
    printf(STDERR "Object #%s (which is a %s) has %d part",
           $shapeNum, $type, $partCount) if ($debugLevel > 1);
    if (($partCount > 1) || ($partCount == 0)) {
      printf(STDERR "s") if ($debugLevel > 1);
    }
    printf(STDERR " and %d point", $pointCount) if ($debugLevel > 1);
    if (($pointCount > 1) || ($pointCount == 0)) {
      printf(STDERR "s") if ($debugLevel > 1);
    }
    my ($sovName, $countryName, $regionName) = shapeNames(\%data);
    if (keys(%onlyNames) &&
        !$onlyNames{$countryName} && !$onlyNames{$sovName}) {
      # if this country shouldn't be displayed, then don't proceed further
      next;
    }
    if ($shapeNum > 3) {
      print(STDERR ".") if ($debugLevel> 0);
    }
    # replace spaces if necessary so group IDs aren't messed up too much
    $sovName =~ s/ /_/g;
    $countryName =~ s/ /_/g;
    $regionName =~ s/ /_/g;
    # sort out sovereign territory / country code mismatch.  If the
    # country and sovereign territory differ, add country to the
    # group of the sovereign territory
    my $countryGroup = 0;       # false
    my $shapeGroup = 0;         # false
    my $sovGroup = 0;           # false
    my $sovString = sprintf("g%s",$sovName);
    my $countryString = sprintf("g%s",$countryName);
    my $regionString = sprintf("g%s",$regionName);
    if ($sovName eq "") {
      # can't find a country name, so improvise with a unique ID
      $sovString = sprintf("gF%dS%d", $fileNum, $shapeNum);
      $countryString = $sovString;
    }
    $sovGroup = $worldGroup->getElementByID($sovString);
    if (!$sovGroup) {
      $sovGroup = $worldGroup->group('id' => $sovString);
    }
    if ($sovName eq $countryName) {
      $countryGroup = $sovGroup;
    } else {
      $countryGroup = $sovGroup->getElementByID($countryString);
      if (!$countryGroup) {
        $countryGroup = $sovGroup->group('id' => $countryString);
      }
    }
    if ($regionName) {          # a region inside a country
      $shapeGroup = $countryGroup->getElementByID($regionString);
      if (!$shapeGroup) {
        $shapeGroup = $countryGroup->group('id' => $regionString);
      }
    } else {
      $shapeGroup = $countryGroup;
    }
    # iterate through component parts of shape
    for my $partNum (1 .. $partCount) { # 1-based indexing
      my @newPart = $shape->get_part($partNum);
      @newPart = project($projOpts, \@newPart);
      @newPart = clip($projOpts, \@newPart);
      @newPart = simplify($projOpts, \%addedPoints, \@newPart);
      my @partBoundBox = ();
      if(scalar(@newPart) < 20){
        @partBoundBox = boundBox(@newPart);
      }
      my $polyID = sprintf("polyF%dS%dP%d", $fileNum, $shapeNum, $partNum);
      my $partClass = " $countryName";
      if($countryName ne $sovName){
        $partClass = " $sovName $countryName";
      }
      my $fillCol = $landColour;
      my $strkCol = $borderColour;
      if (($politicalNames{$sovName}) || ($politicalNames{$countryName}) ||
          ($politicalNames{$regionName})) {
        # order makes sure if subject is political as well, it
        # will be coloured as a subject
        $fillCol = $othLandColour;
        $strkCol = $polBordColour;
        $partClass .= ' political';
      }
      if (($subjectNames{$sovName}) || ($subjectNames{$countryName}) ||
          ($subjectNames{$regionName})) {
        $fillCol = $subLandColour;
        $strkCol = $intBordColour;
        $partClass .= ' subject';
      }
      $partClass =~ s/ +$//;
      if ($type eq "Polygon") {
        if (@newPart) {
          my $pointString = chopPoly($projOpts, \@newPart, 1);
          if (@partBoundBox &&
              ($partBoundBox[2] - $partBoundBox[0] < 10) &&
              ($partBoundBox[3] - $partBoundBox[1] < 10)) {
            my $point = [($partBoundBox[0] + $partBoundBox[2]) / 2,
                         ($partBoundBox[1] + $partBoundBox[3]) / 2,1];
            my $pointID = sprintf("pointF%dS%dP%d", $fileNum,
                                  $shapeNum, $partNum);
            if($partClass =~ /(subject|political)/) {
              # subject/political points get a highlight ring placed around them
              my $highlightID = sprintf("pthglF%dS%dP%d", $fileNum,
                                        $shapeNum, $partNum);
              $shapeGroup->circle('cx' => (($point->[0]) + $projOpts->{"xAdj"}),
                                  'cy' => (($point->[1]) + $projOpts->{"yAdj"}),
                                  'r' => $projOpts->{"pointSize"} * 10,
                                  'id' => $highlightID,
                                  'opacity' => 0.25,
                                  'class' => 'highlight'.$partClass);
            }
            if(scalar(@newPart) == 1){
              $shapeGroup->circle('cx' => (($point->[0]) + $projOpts->{"xAdj"}),
                                  'cy' => (($point->[1]) + $projOpts->{"yAdj"}),
                                  'r' => $projOpts->{"pointSize"},
                                  'id' => $pointID,
                                  'class' => 'region'.$partClass);
            }
          }
          if ($pointString && ($pointString ne "Z")) {
            $shapeGroup->path('d' => $pointString, 'id' => $polyID,
                              'class' => 'region'.$partClass);
          }
        }
      } elsif ($type eq "PolyLine") {
        # assume this is a river
        if (@newPart) {
          # polylines don't get closed
          my $pointString = chopPoly($projOpts, \@newPart, 0);
          if ($pointString) {
            $shapeGroup->path('d' => $pointString, 'id' => $polyID,
                              'class' => 'river'.$partClass);
          }
        }
      } else {
        die("Cannot deal with shape type '$type'");
      }
      printf(STDERR " (%s points)", scalar(@newPart)) if ($debugLevel> 1);
    }
    if ($partCount == 0) {
      if ($type eq "Point") {
        # used for small countries
        my @shapePoints = $shape->points();
        my $pointNum = 0;
        foreach my $oldPoint (@shapePoints) {
          $pointNum++;
          my @pointArr = ();
          push(@pointArr, $oldPoint);
          @pointArr = project($projOpts, \@pointArr);
          @pointArr = clip($projOpts, \@pointArr);
          my $point = $pointArr[0];
          if ($point) {
            # it is assumed that points are so small that
            # displaying a separate border would overwhelm
            # the object
            my $pointID = sprintf("pointF%dS%dP%d", $fileNum,
                                  $shapeNum, $pointNum);
            my $partClass = " $countryName";
            if($countryName ne $sovName){
              $partClass = " $sovName $countryName";
            }
            my $fillCol = $borderColour; # was $landColour
            my $strkCol = 'none';
            if (($politicalNames{$sovName}) ||
                ($politicalNames{$countryName}) ||
                ($politicalNames{$regionName})) {
              # order makes sure if subject is political
              # as well, it will be coloured as a subject
              $fillCol = $othLandColour;
              $strkCol = $polBordColour;
              $partClass .= ' political';
            }
            if (($subjectNames{$sovName}) ||
                ($subjectNames{$countryName}) ||
                ($subjectNames{$regionName})) {
              $fillCol = $subLandColour;
              $strkCol = $intBordColour;
              $partClass .= ' subject';
            }
            if ($partClass =~ /(subject|political)/) {
              # subject points get a highlight ring placed around them
              my $highlightID = sprintf("pthglF%dS%dP%d", $fileNum,
                                        $shapeNum, $pointNum);
              $shapeGroup->circle('cx' => (($point->[0]) + $projOpts->{"xAdj"}),
                                  'cy' => (($point->[1]) + $projOpts->{"yAdj"}),
                                  'r' => $projOpts->{"pointSize"} * 10,
                                  'id' => $highlightID,
                                  'opacity' => 0.25,
                                  'class' => 'highlight'.$partClass);
            }
            $shapeGroup->circle('cx' => (($point->[0]) + $projOpts->{"xAdj"}),
                                'cy' => (($point->[1]) + $projOpts->{"yAdj"}),
                                'r' => $projOpts->{"pointSize"},
                                'id' => $pointID,
                                'class' => 'region'.$partClass);
          }
        }
      }
    }
    print(STDERR "\n") if ($debugLevel> 1);
  }                             # ends for(shapenum)
  print(STDERR " done!\n") if ($debugLevel> 0);
}

$latLimit = 80 unless $latLimit;

if ($printLines && !$projOpts->{"zoomed"}) {
  print(STDERR "Printing lines of latitude and longitude...")
    if ($debugLevel> 0);
  # printLines is a bit-packed integer
  # 1: print Latitude [only]
  # 2: print Longitude [only]
  # 2: print both Latitude and Longitude
  my $lineGroup = $svg->group('id' => "gLatLongLines");
  # min and max lat make sure there are no sticks extending below the
  # lines of latitude
  my $minLat = -$latLimit + (($latSep - (-$latLimit % $latSep)) % $latSep);
  my $maxLat = $latLimit - (($latLimit % $latSep) % $latSep);
  if (($printLines & 1) != 0) {
    # % < 1 allows for non-integer separations, but makes sure always
    # on a degree line
    for my $latPoint (grep (($_ % $latSep) < 1, $minLat..$maxLat)) {
      print(STDERR ".") if ($debugLevel> 0);
      my @linePoints = ();
      for (my $longPoint = -180; $longPoint <= 180; $longPoint += 0.1) {
        push(@linePoints, [$longPoint, $latPoint]);
      }
      my $polyID = sprintf("lat%d", $latPoint);
      @linePoints = project($projOpts, \@linePoints);
      @linePoints = simplify($projOpts, \%addedPoints, \@linePoints);
      my $pathString = chopPoly($projOpts, \@linePoints, 0);
      if ($pathString && ($pathString ne "Z")) {
        $lineGroup->path('d' => $pathString, 'id' => $polyID, 'opacity' => 0.5,
                         'class' => 'grid');
      }
    }
  }
  if (($printLines & 2) != 0) {
    for my $longPoint (grep (($_ % $longSep) < 1, -180..179)) {
      print(STDERR ".") if ($debugLevel> 0);
      my @linePoints = ();
      for (my $latPoint = $minLat; $latPoint <= $maxLat; $latPoint += 0.1) {
        push(@linePoints, [$longPoint, $latPoint]);
      }
      my $polyID = sprintf("long%d", $longPoint);
      @linePoints = project($projOpts, \@linePoints);
      @linePoints = simplify($projOpts, \%addedPoints, \@linePoints);
      my $pathString = chopPoly($projOpts, \@linePoints, 0);
      if ($pathString && ($pathString ne "Z")) {
        $lineGroup->path('id' => $polyID, 'd' => $pathString, 'opacity' => 0.5,
                         'class' => 'grid');
      }
    }
  }
  print(STDERR " done!\n") if ($debugLevel> 0);
}

if ($projOpts->{"zoomed"}) {
  my $seaGroup = $svg->group('id' => "gGlobeDecoration");
  my @linePoints = ([$projOpts->{"minX"},$projOpts->{"minY"},1],
                    [$projOpts->{"maxX"},$projOpts->{"minY"},1],
                    [$projOpts->{"maxX"},$projOpts->{"maxY"},1],
                    [$projOpts->{"minX"},$projOpts->{"maxY"},1]);
  $seaPath = chopPoly($projOpts, \@linePoints, 1);
  $seaGroup->path('id' => 'pGlobeBorder', 'd' => $seaPath,
                  'class' => 'mapborder');
  printf(STDERR "Zoomed to ((%0.2f,%0.2f)-(%0.2f,%0.2f))!\n",
         $projOpts->{"minX"}, $projOpts->{"minY"},
         $projOpts->{"maxX"}, $projOpts->{"maxY"}) if ($debugLevel> 0);
}

if (($mapType =~ /^ortho/) && !$projOpts->{"zoomed"}) {
  print(STDERR "Printing border...") if ($debugLevel> 0);
  my $seaGroup = $svg->group('id' => "gGlobeDecoration");
  # shading
  $seaGroup->circle('id' => "cGlobeShade", 'cx' => $projOpts->{"xAdj"},
                    'cy' => $projOpts->{"yAdj"}, 'r' => ($projOpts->{"svgWidth"} / 2),
                    'class' => 'mapshadow');
  $seaGroup->circle('id' => "cGlobeBorder", 'cx' => $projOpts->{"xAdj"},
                    'cy' => $projOpts->{"yAdj"}, 'r' => ($projOpts->{"svgWidth"} / 2),
                    'class' => 'mapborder');
  print(STDERR " done!\n") if ($debugLevel> 0);
}

if (($mapType =~ /^(location|locator|area|world)/) &&
    !$projOpts->{"zoomed"}) {
  print(STDERR "Printing border...") if ($debugLevel> 0);
  my $seaGroup = $svg->group('id' => "gGlobeDecoration");
  $seaGroup->path('id' => 'pGlobeBorder', 'd' => $seaPath,
                  'class' => 'mapborder');
  print(STDERR " done!\n") if ($debugLevel> 0);
}

# show a heatmap key if heatmap colours are used
if (keys(%countryValues) && $printKey) {
  print(STDERR "Printing heatmap key...") if ($debugLevel> 0);
  my $svgDefs = $svg->defs('id' => 'defHeatmap');
  my $hmGradient = $svgDefs->gradient('id' => "grHeatmap", '-type' => "linear",
                                      'fx' => 0.625, 'fy' => 0.375,
                                      'class' => 'seabase gradient');
  $hmGradient->stop('id' => 'sMin', offset => 0,
                    'style' => 'stop-color:#00FFFF;');
  $hmGradient->stop('id' => 'sMax', offset => 1,
                    'style' => 'stop-color:#FFFF00;');
  my $gradientGroup = $svg->group('id' => 'gHeatmapScale');
  $gradientGroup->rect('id' => 'rHeatmapScale', 'class' => 'heatmap gradient',
                       'x' => ($projOpts->{"svgWidth"} / 2) - 150, 'y' => $projOpts->{"svgHeight"}-40,
                       'width' => 300, 'height' => 35);
  $gradientGroup->text('id' => 'tScaleMin', 'class' => 'heatmap text min',
                       'x' => ($projOpts->{"svgWidth"} / 2) - 95, 'y' => $projOpts->{"svgHeight"}-15,
                       -cdata => $minCValue, 'style' => 'text-anchor: start');
  $gradientGroup->text('id' => 'tScaleMax', 'class' => 'heatmap text max',
                       'x' => ($projOpts->{"svgWidth"} / 2) + 95, 'y' => $projOpts->{"svgHeight"}-15,
                       -cdata => $maxCValue, 'style' => 'text-anchor: end');
  print(STDERR " done!\n") if ($debugLevel> 0);
}

print $svg->xmlify();
