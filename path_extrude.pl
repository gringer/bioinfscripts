#!/usr/bin/perl
use warnings;
use strict;

use Math::Trig qw(acos :radial :pi);

$\ = $/;

our $VERSION = "2013.11.20.0";

=pod

=head1 DESCRIPTION

Extrude a polygon along a path and produce a polyhedron definition
suitable for compiling in OpenSCAD. The polygon should be closed,
defined by the points (in clockwise order) along the first line of
input. The path is defined by the points along the second line of
input.

[points are represented as a reference to an array]

At each point along the path, the polygon will be placed and rotated
so that a line from the previous point to the next point will be
parallel to the normal vector of the polygon. The extrusion is carried
out by generating two triangles that connect polygon edges between
different adjacent points on the path. End faces will be placed at the
start and end points and rotated to face the next or previous point
respectively.

If the start and end points of the path are exactly the same (by text
matching), then the created polyhedron is joined at the start / end
point as if it were a mid point, and the end faces are not generated.

=head1 SYNOPSIS

cat arguments.txt | ./path_extrude.pl > out.scad

The input file is a line-delimited file of the form B<option> =
B<value>. See L<ARGUMENTS> for details.

=head1 ARGUMENTS

Arguments are specified in the input file as B<option> = B<value>
(e.g. C<polygon = [[-1,-1],[-1,1],[1,1],[1,-1]]>).

=over 2

B<polygon> A list of 2D polygon points in clockwise order

B<path> A list of 3D path points in the order of extrusion

B<pathfn> A function returning a 3-element array that depends on the
parameter t = 0 .. 2*pi

B<$fn> The number of path points to generate (used with I<pathfn>)

=back

Points are represented as an OpenSCAD array of arrays (i.e.
[[x1,y1,z1],[x2,y2,z2],...] for 3D points).

The output of this script is an OpenSCAD polyhedron definition that
extrudes the polygon over the specified path. The polygon normal at
each point in the path is the same as the vector between the path
point and the next path point.

=head1 CAVEATS

When the shift to the next point is too extreme, the polygon may be
rotated incorrectly and triangles joining edges may intersect with one
another. This may make certain programs sad.

The ends of the extruded polygon are closed off by creating triangles
from adjacent polygon points to the bounding box centre. It is
therefore assumed that all rays travelling from the origin to infinity
will either intersect with a single point, or cross a single line of
the polygon. This assumption can be safely violated if the start and
end points are the same, because in that case no ends are filled in.

=head1 METHODS

=cut

=head2 pathDiff(I<3DPoint>, I<3DPoint>)

Calculate the difference vector between two rectangular points.

=cut

sub pathDiff{
  my ($x0, $y0, $z0) = @{$_[0]};
  my ($x1, $y1, $z1) = @{$_[1]};
  my @out = ($x0 - $x1, $y0 - $y1, $z0 - $z1);
  return(\@out);
}

=head2 rToS(I<pointArray>)

Convert all points in the array from rectangular (cartesian) to
spherical coordinates. If the array describes 2D points, the I<Y>
component is set to 0.

=cut

sub rToS{
  my @point = @{$_[0]};
  # convert 2D to 3D by setting Z to 0
  # this makes the plane normal become 1,0,0[sph]
  if(@point < 3){
    push(@point, 0);
  }
  my ($x, $y, $z) = @point;
  my @out = cartesian_to_spherical($x,$y,$z);
  #my $r = sqrt($x*$x+$y*$y+$z*$z);
  #my @out = ($r,atan2($y,$x),acos($z/$r));
  printf(STDERR "%0.3f,%0.3f,%0.3f -> %0.3f,%0.3f,%0.3f\n",
         @point, @out);
  return(\@out);
}

=head2 rToS(I<pointArray>)

Convert all points in the array from spherical to rectangular
coordinates.

=cut

sub sToR{
  my @point = @{$_[0]};
  my ($r, $theta, $phi) = @point;
  my @out = spherical_to_cartesian($r, $theta, $phi);
  grep{$_ = sprintf("%0.3f", $_) } @out; # round to 3 d.p.
  return(\@out);
}

=head2 rotate(I<3DPoint>,I<pointArray>)

Rotate all XY-planar rectangular points in the array about the origin
by the angles defined in the difference vector.

See L<http://www.youtube.com/watch?v=cUv-5unQxtE>

=cut

sub rotate{
  my ($rDiff, @rectPoints) = @_;
  my (undef, $theta, $phi) = cartesian_to_spherical(@{$rDiff});
  my @newPoints = map {
    my @n = @{$_};
    die("Can only rotate 2D points") if(@n != 2);
    # set Z to 0 to point normal up (same direction as sph(1,0,0)
    push(@n, 0);
    # rotate around y axis to point normal to x axis
    @n = ($n[0] * cos($phi) + $n[2] * sin($phi),
          $n[1],
         -$n[0] * sin($phi) + $n[2] * cos($phi));
    # rotate around z axis to point normal to difference vector on XY plane
    @n = ($n[0] * cos($theta) - $n[1] * sin($theta),
          $n[0] * sin($theta) + $n[1] * cos($theta),
          $n[2]);
    grep{$_ = sprintf("%0.3f", $_) } @n; # round to 3 d.p.
    \@n;
  } @rectPoints;
  return(@newPoints);
}

=head2 translate(I<3DPoint>,I<pointArray>)

Translate all rectangular points in the array by the I<X>, I<Y> and
I<Z> components of the rectangular point / vector.

=cut

sub translate{
  my ($tPos, @rectPoints) = @_;
  my @newPoints = map {
    my @point = @{$_};
    $point[0] += $tPos->[0];
    $point[1] += $tPos->[1];
    $point[2] += $tPos->[2];
    \@point;
  } @rectPoints;
  return(@newPoints);
}

=head2 vecMean(I<3DPoint>,I<3DPoint>)

Calculate the mean of two vectors defined by rectangular coordinates.

=cut

sub vecMean{
  my ($v1, $v2) = @_;
  my @res = @{$v1};
  $res[$_] += $v2->[$_] for 0 .. $#res;
  return(\@res);
}

=head2 equal(I<3DPoint>,I<3DPoint>)

Determine if two points are equal.

=cut

sub equal{
  my ($p1, $p2) = @_;
  for(my $i = 0; $i < @{$p1}; $i++){
    if($p1->[$i] ne $p2->[$i]){
      return 0; # false
    }
  }
  return 1; # true
}


=head2 pointToText(I<point>)

Convert a point into an OpenSCAD array representation (comma
separated, with square brackets)

=cut

sub pointToText{
  my @out = @{$_[0]};
  return sprintf("[%s".(",%s" x (@out-1))."]", @out);
}

my ($polyLine, $pathLine) = ("","");
my @polyPoints = ();
my @pathPoints = ();
my $numPoints;
my $pathFunction = ""; # function depending on t = 0 .. 2*pi

while(<>){
  if(/^poly(gon)?\s*=\s*(.*)$/i){
    $polyLine = $2;
  }
  if(/^path\s*=\s*(.*)$/i){
    $pathLine = $1;
  }
  if(/^pathf(unctio)?n\s*=\s*(.*)$/i){
    $pathFunction = $2;
  }
  if(/^\$fn\s*=\s*(.*)$/){
    $numPoints = $1;
  }
}

# convert polygon into 2D point notation
while($polyLine =~ s/\[?\[([^\]]+)\],?//){
    my @point = split(/,/, $1);
    if(@point != 2){
        die("Polygon points should be 2D points: [[x1,y1],[x2,y2],...]");
    }
    push(@polyPoints, \@point);
}

# convert path into internal 3D point notation
while($pathLine =~ s/\[?\[([^\]]+)\],?//){
    my @point = split(/,/, $1);
    if(@point != 3){
        die("Path points should be 3D points: [[x1,y1,z1],[x2,y2,z2],...]");
    }
    push(@pathPoints, \@point);
}

if($pathFunction && !$numPoints){
  $numPoints = 20; # default number of points when unspecified
}


if($pathFunction){
    printf(STDERR "Processing function: %s\n", $pathFunction);
  my $tInc = (2 * pi) / ($numPoints - 1);
  # end point is trickier because floating point equality is tricky
  for(my $t = 0; $t < (2 * pi + $tInc/2); $t += $tInc){
    my @point = eval($pathFunction);
    if(@point != 3){
        print(STDERR @point);
        die("Path point function should return 3D points, f(t) = (x,y,z)");
    }
    grep{$_ = sprintf("%0.3f", $_); $_ = 0 if ($_ == 0) } @point; # round to 3 d.p.
    push(@pathPoints, \@point);
  }
}

if($numPoints && (@pathPoints != $numPoints)){
  die("Path points specified, but count doesn't match \$fn");
}

if(@pathPoints < 2){
  die("Need at least two points along a path for extrusion");
}

my ($rDiffNext, $rDiffLast) = (0, 0);

my @scadPoints = ();

for(my $p = 0; $p < @pathPoints; $p++){
  # add a polygon at each path point, rotated to face the next point
  if(($p+1) < @pathPoints){
    $rDiffNext = pathDiff($pathPoints[$p+1],$pathPoints[$p]);
  } elsif(equal($pathPoints[0], $pathPoints[$#pathPoints])) {
    $rDiffNext = pathDiff($pathPoints[1],$pathPoints[$p]);
  }
  if(($p) > 0){
    $rDiffLast = pathDiff($pathPoints[$p],$pathPoints[$p-1]);
  } elsif(equal($pathPoints[0], $pathPoints[$#pathPoints])) {
    $rDiffLast = pathDiff($pathPoints[0],$pathPoints[$#pathPoints-1]);
  } else {
    $rDiffLast = $rDiffNext;
  }
  # the actual rotation is the mean of the rotation from the previous
  # point and the rotation to the next point
  my $rDiff = vecMean($rDiffLast, $rDiffNext);
  push(@scadPoints,
       join(",",map {pointToText($_)}
            translate($pathPoints[$p],rotate($rDiff, @polyPoints))));
}
# add on centre points for ends
push(@scadPoints, pointToText($pathPoints[0]));
push(@scadPoints, pointToText($pathPoints[$#pathPoints]));
grep {s/\],\[/\],\n            \[/g} @scadPoints;

print("polyhedron(");
print("  points = [".join(",\n            ", @scadPoints)."],");

my @trianglePoints = ();

# start surface
if(!equal($pathPoints[0], $pathPoints[$#pathPoints])){
  for(my $i = 0; $i < @polyPoints; $i++){
    push(@trianglePoints,
         sprintf("[%d,%d,%d]", ($i+1)%(@polyPoints),
                 $i, @polyPoints*@pathPoints));
  }
}

for(my $pi = 0; $pi < @pathPoints - 1; $pi++){
  my $p = $pi * @polyPoints;
  for(my $i = 0; $i < @polyPoints; $i++){
    push(@trianglePoints,
         # connecting adjacent edges with two triangles
         sprintf("[%d,%d,%d],", $i+$p, ($i+1)%(@polyPoints)+$p,
                                  ($i+1)%(@polyPoints)+@polyPoints+$p).
         sprintf("[%d,%d,%d]",($i+@polyPoints)+$p,
                                  $i+$p, ($i+1)%(@polyPoints)+@polyPoints+$p));
  }
}

# end surface
if(!equal($pathPoints[0], $pathPoints[$#pathPoints])){
  my $p = (@pathPoints - 1) * @polyPoints;
  for(my $i = 0; $i < @polyPoints; $i++){
    push(@trianglePoints,
         sprintf("[%d,%d,%d]", $i+$p,
                 ($i+1)%(@polyPoints)+$p,
                 @polyPoints*@pathPoints+1));
  }
}

print("  triangles = [".join(",\n               ", @trianglePoints)."]);");

=head1 LICENSE

Copyright (c) 2013, David Eccles (gringer) <bioinformatics@gringene.org>

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

=cut
