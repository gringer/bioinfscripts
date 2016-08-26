#!/usr/bin/perl

use warnings;
use strict;

use Pod::Usage; ## uses pod documentation in usage code
use Getopt::Long qw(:config auto_version auto_help pass_through); # for option parsing

use Math::Trig qw(asin_real acos_real pi);
# http://search.cpan.org/~jasonk/Geo-ShapeFile-2.52/lib/Geo/ShapeFile.pm
use Geo::ShapeFile;
use SVG;

use POSIX qw(fmod);

## more verbose traceback on errors
use Carp 'verbose'; 
$SIG{ __DIE__ } = \&Carp::confess;

our $VERSION = "1.96";

=head1 NAME

perlshaper.pl - convert shapefile data to SVG format, for use in
wikimedia maps.

=head1 SYNOPSIS

./perlshaper.pl <shapefile(s)> -type <mapType> [options] > output.svg


=head2 Basic Options

=over 2

=item B<-help>

Only display this help message

=item B<-type> I<string>

Type of map (location|locator|area|world|orthographic)

=item B<-round> I<int>

Round values (in SVG file) to given number of decimal places

=item B<-centre> I<long>,I<lat>

Identify centre of map (by longitude/latitude)

=item B<-centre> I<ISO 3A-code>

Identify centre of map (by target country). The target can be specified as 'I<sovCode>/I<countryCode>' for a stricter match.

=item B<-data> I<file>

Colour countries based on numerical data in file

=back


=head1 ARGUMENTS

Most of the map-based variables used in the code can be customised to
fit particular needs using command-line arguments.

=head2 Advanced Options

=over 2

=item B<-psize> I<float>

Radius of small points

=item B<-colour> I<string>

Change colour theme to a different map type

=item B<-proj> I<string>

Change target projection

=item B<-landcol> I<string>

Land colour

=item B<-seacol> I<string>

Sea colour

=item B<-bordcol> I<string>

Border colour

=item B<-sub> I<ISO 3A-code>

Identify subject region

=item B<-pol> I<ISO 3A-code>

Identify related political region

=item B<-only> I<ISO 3A-code>

Only display specified shape(s)

=item B<-zoom> I<ISO 3A-code>

Zoom to projected extents of shape(s)

=item B<-zoom> I<minX,Y,maxX,Y>

Zoom to projected extents, typically with I<X> = longitude, I<Y> = latitude

=item B<-nokey>

Don't display heatmap key (if heatmap is used)

=item B<-[no]lines>

[don't] print lines of latitude and longitude

=item B<-v>

Verbose output (also '-v -v')

=back

=head1 DESCRIPTION

This code is designed with wikimedia in mind, so the default options
should produce SVG images that follow the general map conventions (see
http://en.wikipedia.org/wiki/Wikipedia:WikiProject_Maps/Conventions). In
particular, the expected input format is that of Natural Earth files
(http://www.naturalearthdata.com).

Just as a head's up, latitude lines are horizontal in the typical
projection of the world, while longitude lines are
vertical. Latitude is usually represented by a number (in degrees)
in the range -90..90, while longitude is usually represented by a
number (in degrees) in the range -180..180.

  /---------\  --- latitude (Y)
 | | | | | | |
 |-|-|-|-|-|-|
 | | | | | | |
  \---------/
     |
     \-- longitude (X)

=head1 METHODS

=cut

=head2 transform(optionRef, point)

Transform a point from {[-180,180],[-90,90]} to
{[-0.5,0.5],[-0.5,0.5]} relative to the desired projection. This
replaces the proj4 transform function with something that will
always return a valid point. Orthographic transformations that
appear on the opposite side of the globe will be converted to a
point with the same angle, but sitting at the circle
edge. Adjustments need to be made post-transform (e.g. multiply by
SVG width, add to fit in 0..max range) to convert to SVG
coordinates.

=cut

sub transform {
  my ($options, $inPoint) = @_;
  my $oldLong = $inPoint->[0];
  my $oldLat = $inPoint->[1];
  my $projection = $options->{"projection"};
  my $lambda = fmod(($oldLong - $options->{"centreLn"} + 540), 360) - 180;
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
    my $alpha = acos_real(cos($phi) * cos($lambda / 2));
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

=head2 project(optionRef, pointRef)

Convert a set of points from standard equirectangular projection
("+proj=latlong +datum=WGS84" in proj4-speak) to another
projection. Natural Earth Data uses this equirectangular format, where
coordinates are specified based on their latitude and longitude
locations in degrees.

This method attempts to fit the map into a rectangle (or square)
that fits in a 1100x550 box.

=cut

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
  my $maxDist2 = ($projWidth / 4) ** 2; # maximum distance between points
  foreach my $inPoint (@input) {
    my $newLong = 0;
    my $newLat = 0;
    if ((ref($inPoint) eq "Geo::ShapeFile::Point") || 
	(ref($inPoint) eq "HASH")) {
      $newLong = $inPoint->X;
      $newLat = $inPoint->Y;
    } elsif (ref($inPoint) eq "ARRAY") {
      $newLong = $inPoint->[0];
      $newLat = $inPoint->[1];
    } else {
	my $ptRefType = ref($inPoint);
	die("Error: expecting a reference to an array or a hash [actual reference type: $ptRefType]");
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
    my $xd = $px - $oldX;
    if (($xd * $xd) > $maxDist2) {
      # don't connect lines that have wrapped around the map (over half the image width)
      if(($options->{"zoomed"}) ||
         ($options->{"projection"} eq "orthographic")){
        # zoomed and orthographic projections don't get additional points added
        ### removed to fix latitude/longitide lines not printing
        # push(@output, 0);
      } else {
        # add additional points on the border edge
        my ($minLong, $maxLong) = ($options->{"minLong"}, $options->{"maxLong"});
        my $oldEdgePoint =
          transform($options, [($px > $oldX) ? $minLong : $maxLong, $oldLat]);
        my $newEdgePoint =
          transform($options, [($px > $oldX) ? $maxLong : $minLong, $newLat]);
        $oldEdgePoint->[0] = $oldEdgePoint->[0] * $xScale * $projWidth;
        $newEdgePoint->[0] = $newEdgePoint->[0] * $xScale * $projWidth;
        $oldEdgePoint->[1] = $oldEdgePoint->[1] * $yScale * -$projHeight;
        $newEdgePoint->[1] = $newEdgePoint->[1] * $yScale * -$projHeight;
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

=head2 orthDist2(point1, pointP, point2)

calculates the square of the minimum distance between the point (xP,yP) and the line [(x1,y1) - (x2,y2)]

The distance formula is from
L<wolfram|http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html>
with modifications for division by zero, and for points outside the line

=cut

sub orthDist2 {
  my ($point1, $pointP, $point2) = @_;
  my $x1 = $point1 -> [0];
  my $y1 = $point1 -> [1];
  my $xP = $pointP -> [0];
  my $yP = $pointP -> [1];
  my $x2 = $point2 -> [0];
  my $y2 = $point2 -> [1];
  if (($x2 == $x1) && ($y2 == $y1)) {
    # exclude case where denominator is zero
    my $dist2 = ($x1-$xP)**2 + ($y1 - $yP)**2;
    return($dist2);
  }
  my $dist2 = (($x2 - $x1)*($y1 - $yP) - ($x1 - $xP)*($y2 - $y1))**2 /
    (($x2 - $x1)**2 + ($y2 - $y1)**2);
  if ($dist2 == 0) {
    # on the line, but need to consider points outside the line
    # this fixes a problem where the equator line is clipped
    my $p1Dist = (($x1-$xP)**2 + ($y1 - $yP)**2);
    my $p2Dist = (($x2-$xP)**2 + ($y2 - $yP)**2);
    my $p12Dist = (($x1-$x2)**2 + ($y1 - $y2)**2);
    my $sigma = 0.0001;
    if (($p1Dist + $p2Dist) > ($p12Dist + $sigma)) {
      # point is outside the line, use smallest distance from line end points
      $dist2 = ($p1Dist < $p2Dist) ? $p1Dist : $p2Dist;
    }
  }
  return($dist2);
}

=head2 distGr2(point1, point2, res2)

Determines if the distance between the point (x1,y1) and the point
(x2,y2) is greater than the square root of res2.

=cut

sub distGr2 {
  my ($point1, $point2, $res2) = @_;
  my $dx = $point2 -> [0] - $point1 -> [0];
  my $dy = $point2 -> [1] - $point1 -> [1];
  return(($dx * $dx + $dy * $dy) > ($res2));
}


=head2 simplify(optionRef, pointHash, pointRef)

Simplify curves to reduce number of points in SVG file.  This method
should keep points if they have already been included as part of the
simplifcation of another curve, so that shapes that overlap won't
have any gaps. This requires a global hash of already added points.

This is a quick process based around resolving distance -- if points
are closer than the resolution distance to another included point,
then they are removed.

=cut

sub simplify {
  my ($options, $pointHash, $pointRef) = @_;
  if($options->{"roundDP"} == -1){
    return @{$pointRef};
  }
  my $formatString = "%0.".$options->{"roundDP"}."f";
  my $dps = 10**$options->{"roundDP"};
  my $last = undef;
  my ($lastXRounded, $lastYRounded, $xRounded, $yRounded);
  my @input = map {
    if($_){
      $xRounded = int($_->[0] * $dps);
      $yRounded = int($_->[1] * $dps);
      if(!defined($last) || $pointHash->{($xRounded, $yRounded)}){
        $pointHash->{($xRounded,$yRounded)} = 1;
        $last = $_;
        $lastXRounded = $xRounded;
        $lastYRounded = $yRounded;
      } elsif(($xRounded != $lastXRounded) || ($yRounded != $lastYRounded)){
        $pointHash->{($xRounded,$yRounded)} = 1;
        $last = $_;
        $lastXRounded = $xRounded;
        $lastYRounded = $yRounded;
      }
    } else {
      $last = undef;
    }
    $_;
  } @{$pointRef};
  return(@input);
}

=head2 clip(optionRef, pointRef)

clip the shape to the limits of the zoom box (if any). The
algorithm used here is the Sutherland-Hodgman algorithm:
http://en.wikipedia.org/wiki/Sutherland-Hodgman

Note that this algorithm can leave some zero-width polygon
sections in the resultant clip region.

=cut

sub clip {
  my ($options, $pointRef) = @_;
  my $xAdj = $options->{"xAdj"};
  my $yAdj = $options->{"yAdj"};
  my @oldPoints = @{$pointRef};
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

=head2 boundBox(pointList)

Determines post-projection rectangular bounding box for a polygon
(e.g. for determining if borders are useful, or zooming in to a
particular region)

=cut

sub boundBox {
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

=head2 shapeNames(shapeData)

retrieve shape name data from a shape database. For countries, three-letter ISO codes are preferred

=cut

sub shapeNames {
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

=head2 makeRelative(optionRef, pointRef)

Creates relative point definitions (excluding the first point). This
is used to reduce the SVG output file size.

=cut

sub makeRelative {
  my ($options, $pointRef) = @_;
  my $round = ($options->{"roundDP"} != -1);
  my $dps = $round ? 10**$options->{"roundDP"} : 1;
  my $xAdj = $round ? int($options->{"xAdj"} * $dps) : $options->{"xAdj"};
  my $yAdj = $round ? int($options->{"yAdj"} * $dps) : $options->{"yAdj"};
  my @points = map{$_->[0] *= $dps; $_->[1] *= $dps; $_} @{$pointRef};
  my $lastX = 0;
  my $lastY = 0;
  foreach my $point (@points) {
    my $nextX = $round ? int($point->[0] + $xAdj) : ($point->[0] + $xAdj);
    my $nextY = $round ? int($point->[1] + $yAdj) : ($point->[1] + $yAdj);
    $point->[0] = $nextX - $lastX;
    $point->[1] = $nextY - $lastY;
    $lastX = $nextX;
    $lastY = $nextY;
  }
  map{$_->[0] /= $dps; $_->[1] /= $dps} @points;
}

=head2 chopPoly(optionRef, pointRef, joinEnds)

Closes up lines in a polygon that contain discontinuous points -- this
usually happens when part of a polygon goes outside the projection
area. If I<joinEnds> is true, then the ends of the polygon are joined
together.

The input is a sequence of Geo::Points, and the output is a string
that can be used as a path description in an SVG file

=cut

sub chopPoly {
  my ($options, $pointRef, $joinEnds) = @_;
  my @oldPoints = @{$pointRef};
  if (!@oldPoints) {
    return (""); # nothing to process
  }
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
  # note: %f rounds to decimal places, %g rounds to significant figures
  #my $printfString = ($dp > -1) ? ("%.".$dp."f,%.".$dp."f") : ("%f,%f");
  my $printfString = ("%s,%s");

  my @subPaths = ();
  foreach (@contigLists) {
    my @currentList = @{$_};
    if(scalar(@currentList) > 1){ # ignore single-point subpaths
      makeRelative($options, \@currentList);
      my @pointStrs = map {
        sprintf($printfString, ($_->[0]), ($_->[1]));
      } @currentList;
      my $subPath = "M".join("l",@pointStrs).($joinEnds?"Z":"");
      $subPath =~ s/l0,0(?=[^\.])//g; # remove 'no change' relative movements
      push(@subPaths, $subPath);
    }
  }
  my $pointPath = join(" ", @subPaths);

  return($pointPath);
}

my $mapType = "";

my $colourTheme = "";

my $centreCountry = "";

my $landColour = ""; # outside area (or default land colour)
my $seaColour = ""; # ocean / sea / lake
my $borderColour = "";

# Extra colours (only used if land/sea/border not specified)

my $subLandColour = ""; # subject colour
my $othLandColour = ""; # other lands of the same political unity
my $coastColour = ""; # coast along lake, rivers, sea coasts
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

# set default options
my $projOpts =
  {
   "roundDP" => -1,
   "projection" => "", # map projection (equirectangular / orthographic)
   "centreLn" => "",  # projection centre X / longitude
   "centreLt" => "",  # projection centre Y / latitude
   "svgWidth" => "", # width of SVG file, excluding padding
   "svgHeight" => "", # height of SVG file, excluding padding
   "xAdj" => 0, # X adjustment to fit in SVG bounds
   "yAdj" => 0,  # Y adjustment to fit in SVG bounds
   "xScale" => "",  # Scale factor of points (for zoom box)
   "yScale" => "",  # Scale factor of points (for zoom box)
    "padding" => "", # amount of padding to add (for zoom box)
   "lines" => 1, # draw lines of latitude / longitude
   "key" => 1, # print a key for the heatmap
  };

GetOptions($projOpts, 'lines!', 'key!', 'pointSize|psize=f', 'roundDP|round=i',
           'centre|center=s',
           'projection=s',
           'zoomed|zoom=s@',
           'colourTheme|color|colour=s', => \$colourTheme,
           'v|verbose+', => \$debugLevel,
           'mapType|type=s' => \$mapType,
           'landcol=s' => \$landColour,
           'seacol=s' => \$seaColour,
           'bordcol=s' => \$borderColour,
           'subjects|sub=s' => sub { my ($opt,$val) = @_;
				     $subjectNames{$val} = 1},
           'politicals|pol=s' => sub { my ($opt,$val) = @_;
				       $politicalNames{$val} = 1},
           'only=s' => sub { my ($opt,$val) = @_;
			     $onlyNames{$val} = 1},
           'data=s@' => \@dataFiles,
           'man' => sub { pod2usage({-verbose => 3}) }
          );

# process remaining command line arguments (hopefully only shapefiles)
while (@ARGV) {
  my $argument = shift @ARGV;
  if (((-f $argument) && ($argument =~ /\.shp$/)) ||
      (-f ($argument.".shp"))) { # file existence check
    $argument =~ s/\.shp$//;
    printf(STDERR "Adding '%s.shp' to list of input files\n", $argument)
      if ($debugLevel > 0);
    push (@shpFileBases, $argument);
  } elsif (-f $argument) {
    pod2usage({-exitVal => 2, -message => "Error: Invalid file extension for '$argument'. ".
               "Please use files with '.shp' extension\n"});
  } else {
    pod2usage({-exitVal => 3, -message => "Error: Unknown command-line option or non-existent file, ".
            "'$argument'\n", -verbose => 0});
  }
}

if($debugLevel > 0){
  print(STDERR "Command-line options:\n");
  foreach my $key (keys(%{$projOpts})){
    if(ref($projOpts->{$key}) eq 'ARRAY'){
      printf(STDERR "  %s: %s\n", $key, join("; ", @{$projOpts->{$key}}));
    } elsif($projOpts->{$key} ne "") {
      printf(STDERR "  %s: %s\n", $key, $projOpts->{$key});
    }
  }
}

if($projOpts->{"lines"}){
  $projOpts->{"lines"} = 3;
}

if(scalar(@shpFileBases) == 0){
  print(STDERR "Error: No input files specified\n");
  pod2usage(2);
}

if(keys(%onlyNames)){
  print(STDERR "Will only display the following features on the map:\n");
  print(STDERR "  ".join("\n  ",keys(%onlyNames))."\n");
}

if(@dataFiles){
  warn("Will extract numerical data (assuming <ISO 3A-code>,<float>) from the following files:\n");
  foreach my $fName (@dataFiles){
    warn("  $fName\n");
  }
}

if($projOpts->{"zoomed"}){
  my @zoomArgs = @{$projOpts->{"zoomed"}};
  foreach my $newArg (@zoomArgs){
    if ($newArg =~ /([0-9\.\-]+),([0-9\.\-]+),([0-9\.\-]+),([0-9\.\-]+)/) {
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
      print(STDERR "Zooming map to include '$newArg'\n");
    }
  }
}


if($mapType !~ /^(location|locator|area|world|orthographic)$/){
  pod2usage({ -exitval => 1, -message => "Error: '".$projOpts->{"mapType"}.
              "' is not a valid map type. ".
              "Only accepts '(location|locator|area|world|orthographic)'",
              -verbose => 0});
}

if($projOpts->{"centre"}){
  my $newArg = $projOpts->{"centre"};
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
}

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
      if (($countryName eq $centreCountry) || ($sovName eq $centreCountry) ||
	  ("$sovName/$countryName" eq $centreCountry)) {
        # found it, so find polygon extents in degrees
        my $shape = $shp->get_shp_record($shapeNum);
        my @shapePoints = $shape->points();
        foreach my $point (@shapePoints) {
          my $px = $point->X;
          my $py = $point->Y;
          if (!defined($minX)) {
            ($minX, $minY, $maxX, $maxY) = ($px, $py, $px, $py);
          }
          if (($px - $minX) > 270) { # deal with pesky -180/180 wrapping issues
            $px -= 360;
          }
          if (($maxX - $px) > 270) {
            $px += 360;
          }
          $minX = $px if($px < $minX);
          $minY = $py if($py < $minY);
          $maxX = $px if($px > $maxX);
          $maxY = $py if($py > $maxY);
        }
        my $midLong = ($minX + $maxX) / 2;
        if ($midLong > 180) {
          $midLong -= 360;
        }
        if ($midLong < -180) {
          $midLong += 360;
        }
        my $midLat = ($minY + $maxY) / 2;
        $projOpts->{"centreLn"} = $midLong;
        $projOpts->{"centreLt"} = $midLat;
	if($countryName eq $sovName){
	    printf(STDERR "Found centre feature for ".
		   "'%s' at point (%s,%s)... ",
		   $countryName,
		   $midLong, $midLat) if ($debugLevel> 0);
	} else {
	    printf(STDERR "Found centre feature for ".
		   "'%s/%s' at point (%s,%s)... ",
		   $sovName, $countryName,
		   $midLong, $midLat) if ($debugLevel> 0);
	}
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
  pod2usage(3);
}

if (($projOpts->{"centreLn"} eq "") ||
    ($projOpts->{"centreLt"} eq "")) {
  $projOpts->{"centreLn"} = 0;
  $projOpts->{"centreLt"} = 0;
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
      if (($zoomNames{$countryName}) || ($zoomNames{$sovName}) ||
	  ($zoomNames{"$sovName/$countryName"})) {
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
	  my $findString = ($sovName eq $countryName) ? $sovName
	      : "$sovName/$countryName";
	  $findString .= " ($regionName)" if $regionName;
          printf(STDERR "Found extent for '$findString', ".
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
  if($projOpts->{"manualZoom"}){
    warn("Detected a manual zoom\n") if $debugLevel > 0;
    my @tmpPoints = ([$projOpts->{"minX"},$projOpts->{"minY"}],
                     [$projOpts->{"maxX"},$projOpts->{"maxY"}]);
    my @projectedPoints = project($projOpts, \@tmpPoints);
    printf(STDERR "  min: (%0.2f,%0.2f) -> (%0.2f,%0.2f)\n", 
	   $tmpPoints[0][0], $tmpPoints[0][1],
	   $projectedPoints[0][0], $projectedPoints[0][1]);
    printf(STDERR "  max: (%0.2f,%0.2f) -> (%0.2f,%0.2f)\n", 
	   $tmpPoints[1][0], $tmpPoints[1][1],
	   $projectedPoints[1][0], $projectedPoints[1][1]);
  }
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
  printf(STDERR "[Relative SVG magnification is (%0.3f,%0.3f)x]\n",
         $relMagX, $relMagY) if ($debugLevel > 0);
  printf(STDERR "[SVG scale is (%0.3f,%0.3f)x]\n",
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

## pre-processed variables have been set up, so can now start writing to the SVG

my $svg = SVG->new('height' => $projOpts->{"svgHeight"} + $projOpts->{"padding"} * 2 + 3,
                   'width' => $projOpts->{"svgWidth"} + $projOpts->{"padding"} * 2 + 3,
                   -indent => "  ",
                  -nocredits => 1); # multi-line credits make Inkscape editing harder
$svg->comment("Created using David Eccles' (gringer) perlshaper.pl script, ".
              "version $VERSION, ".
              "http://en.wikipedia.org/wiki/User:Gringer/perlshaper");
$svg->comment("Generated using the Perl SVG Module V". $SVG::VERSION .
              "by Ronan Oger [see http://www.roitsystems.com]"); # add in credit line
$svg->comment("Command line: ".join(" ",@commandLine));

# extract / colour based on data from files
my $countryColours = "";
my %countryValues = ();
my %groupNames = ();
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
$styleSheet->CDATA( <<_EOCSS_
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
_EOCSS_
                    .$countryColours);

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
    my $point = transform($projOpts, [$longPoint, 0])
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
  for (my $latPoint = -90; $latPoint <= 90; $latPoint += 1) { # left edge
    push(@linePoints, [$minLong - 360, $latPoint]);
  }
  for (my $longPoint = $minLong; $longPoint < $maxLong; $longPoint += 1) { # line across top
    push(@linePoints, [$longPoint - 360, 90]);
  }
  for (my $latPoint = 90; $latPoint >= -90; $latPoint -= 1) { # right edge
    push(@linePoints, [$maxLong - 360, $latPoint]);
  }
  for (my $longPoint = $maxLong; $longPoint > $minLong; $longPoint -= 1) { # line across bottom
    push(@linePoints, [$longPoint - 360, -90]);
  }
  while($minLong > 180){
    $minLong -= 360;
  }
  while($maxLong > 180){
    $maxLong -= 360;
  }
  printf(STDERR " done (%s,%s)!\n", $minLong, $maxLong) if ($debugLevel> 0);
  $projOpts->{"minLong"} = $minLong;
  $projOpts->{"maxLong"} = $maxLong;
  printf(STDERR "Drawing border points... ");
  printf(STDERR "Projecting border points... ");
  @linePoints = project($projOpts, \@linePoints);
  printf(STDERR "simplifying %d border points... ", scalar(@linePoints));
  @linePoints = simplify($projOpts, \%addedPoints, \@linePoints);
  printf(STDERR "chopping border points... ");
  $seaPath = chopPoly($projOpts, \@linePoints, 1);
  $seaGroup->path('id' => 'pSeaBase', 'd' => $seaPath, 'class' => 'seabase');
  printf(STDERR "done!\n");
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
    if (keys(%onlyNames) && !$onlyNames{"$sovName/$countryName"} &&
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
      $groupNames{$sovGroup->{"id"}} = 1;
    }
    if ($sovName eq $countryName) {
      $countryGroup = $sovGroup;
    } else {
      $countryGroup = $sovGroup->getElementByID($countryString);
      if (!$countryGroup) {
        $countryGroup = $sovGroup->group('id' => $countryString);
        $groupNames{$countryGroup->{"id"}} = 1;
      }
    }
    if ($regionName) {          # a region inside a country
      $shapeGroup = $countryGroup->getElementByID($regionString);
      if (!$shapeGroup) {
        $shapeGroup = $countryGroup->group('id' => $regionString);
        $groupNames{$shapeGroup->{"id"}} = 1;
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
          ($politicalNames{"$sovName/$countryName"}) ||
          ($politicalNames{$regionName})) {
        # order makes sure if subject is political as well, it
        # will be coloured as a subject
        $fillCol = $othLandColour;
        $strkCol = $polBordColour;
        $partClass .= ' political';
      }
      if (($subjectNames{$sovName}) || ($subjectNames{$countryName}) ||
          ($subjectNames{"$sovName/$countryName"}) ||
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

__END__

=head1 ERRORS

=head2 Return codes

=over 2

=item B<0>: success

=item B<1>: unknown command-line option

=item B<2>: file error

=item B<3>: projection error

=back

=head1 BUGS

=over 4

=item * The zoom box needs a defined centre to work properly

=back

=head1 AUTHOR

David Eccles (gringer) 2010-2014 <programming@gringer.org>

=head1 LICENSE

Copyright (c) 2010-2014 David Eccles (gringer) <programming@gringer.org>

Permission to use, copy, modify, and/or distribute this software for
any purpose with or without fee is hereby granted, provided that the
above copyright notice and this permission notice appear in all
copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL
WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE
AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL
DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER
TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
PERFORMANCE OF THIS SOFTWARE.

=head1 AVAILABILITY

The most recent version of this code can be found at

https://github.com/gringer/bioinfscripts/blob/master/perlshaper.pl
