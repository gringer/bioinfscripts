#!/usr/bin/perl

use warnings;
use strict;

use Pod::Usage; ## uses pod documentation in usage code
use Getopt::Long qw(:config auto_version auto_help pass_through);
use Math::Trig qw(:pi);

use Graphics::ColorObject;

our $VERSION = "0.5";
our $DEBUG = 0;

=head1 NAME

stl2svg.pl -- converts an STL file into an SVG vector image

=head1 SYNOPSIS

./stl2svg.pl <input.stl>[:#colour]

=head2 Options

=over 2

=item B<-help>

Only display this help message

=back

=head1 DESCRIPTION

This program currently creates an SVG file from a top-down view of an
STL object. There is no orientation / camera setting. This script just
takes a view from the top. Use OpenSCAD (or similar) to re-orient STL
files so that the top-down view is the desired view.

=cut

# https://stackoverflow.com/questions/12472036/how-can-i-convert-rgb-to-lab

sub HL2RGB {
  my ($H, $L) = @_;
  $L = ($L * (2.55-0.78) * 3.1) - 1.5; # normalise for range 0..1
  my $a = cos($H * pi/180) * 2.5; # adjust to make values more saturated
  my $b = sin($H * pi/180) * 2.5;
  return(Lab2RGB($L, $a, $b));
}

# http://www.brucelindbloom.com/index.html?Eqn_Lab_to_XYZ.html
sub Lab2RGB {
  my ($L, $a, $b) = @_;
  my ($Xr, $Yr, $Zr) = _RGB2XYZitu(0.9504,1.0000,1.0888); # D65

  my $epsilon = 216/24389;
  my $kappa = 24389/27;

  my $fy = ($L + 16) / 116;
  my $fx = $a/500 + $fy; my $fx3 = $fx*$fx*$fx;
  my $fz = $fy - $b/200; my $fz3 = $fz*$fz*$fz;

  my $xr = ($fx3 > $epsilon) ? $fx3 : ((116 * $fx - 16) / $kappa);
  my $yr = ($L > ($kappa * $epsilon)) ?
    (($L+16)/116)*(($L+16)/116)*(($L+16)/116) : ($L/$kappa);
  my $zr = ($fz3 > $epsilon) ? $fz3 : ((116 * $fz - 16) / $kappa);

  return _XYZitu2RGB($xr * $Xr, $yr * $Yr, $zr * $Zr);
}

sub _RGB2XYZitu {
  my ($r, $g, $b) = @_;
  return (
    0.431*$r + 0.342*$g + 0.178*$b,
    0.222*$r + 0.707*$g + 0.071*$b,
    0.020*$r + 0.130*$g + 0.939*$b
  );
}

sub _XYZitu2RGB {
  my ($x, $y, $z) = @_;
  my ($r, $g, $b) = map { $_ > 1 ? 255 : ($_ * 255) } (
    3.063*$x - 1.393*$y - 0.476*$z,
   -0.969*$x + 1.876*$y + 0.042*$z,
    0.068*$x - 0.229*$y + 1.069*$z
                                                      );
  $r *= 255; $g *= 255; $b *= 255;
  $r = ($r < 0)?0:(($r>255)?255:$r);
  $g = ($g < 0)?0:(($g>255)?255:$g);
  $b = ($b < 0)?0:(($b>255)?255:$b);
  return((int($r), int($g), int($b)));
}

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

####################################################
# Command line parsing and verification starts here
####################################################

my $argLine = join(" ",@ARGV);

my $options =
  {
   "hue" => 0,
  };

my %fileHues = ();

GetOptions($options,
           'hue=s',
           'debug!' => \$DEBUG,
          ) or pod2usage(1);

grep {if(s/^(.*?)(:([0-9]+))$/$1/){$fileHues{$1} = $3}} @ARGV;

#### Approximate hue angles
## Red: 30; Green: 140; Blue: 250; Yellow: 90; Magenta: 330; Cyan: 200

## For debugging purposes: testing out different lighting values
# for(my $i = 0; $i < 360; $i++){
#   printf(STDERR "#%02x%02x%02x\n", HL2RGB($i, 0));
# }
# for(my $i = 0; $i < 360; $i++){
#   printf(STDERR "#%02x%02x%02x\n", HL2RGB($i, 0.25));
# }
# for(my $i = 0; $i < 360; $i++){
#   printf(STDERR "#%02x%02x%02x\n", HL2RGB($i, 0.5));
# }
# for(my $i = 0; $i < 360; $i++){
#   printf(STDERR "#%02x%02x%02x\n", HL2RGB($i, 0.75));
# }
# for(my $i = 0; $i < 360; $i++){
#   printf(STDERR "#%02x%02x%02x\n", HL2RGB($i, 1));
# }

while (<>) {
  my $fileHue = ($fileHues{$ARGV}) if (exists($fileHues{$ARGV}));
  chomp;
  s/^\s+//;
  if (/^facet normal (([0-9.\-e]+) ([0-9.\-e]+) ([0-9.\-e]+))/) {
    $ns = $1;
    my ($nx, $ny, $nz) = ($2, $3, $4);
    $ny = -$ny;
    $zAng = atan2($ny,$nx); # determine angle on Z axis
    my $colour = sprintf("#%02x%02x%02x",
                         HL2RGB($fileHue, (sin($zAng + 7 * pi/4)+1) / 2));
    $normColours{$ns} = $colour;
    #if ((abs($nx) < 0.01) && (abs($ny) < 0.01)) {
    #  $normColours{$ns} = sprintf("%02x", 240);
    #} else {
    #  $normColours{$ns} = sprintf("%02x",(40 + (abs($zAng) / pi) * 200));
    #}
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
       "stroke-width=\"0.04\" stroke-linejoin=\"round\">\n");
foreach my $zPos (sort {$a <=> $b} (keys(%polys))) {
  foreach my $line (@{$polys{$zPos}}) {
    my ($fx, $fy, $rest) = split(/[, ]/, $line, 3);
    my $colour = shift(@{$colours{$zPos}});
    $fx += $ox;
    $fy += $oy;
    printf("  <path d=\"m%0.4f,%0.4f %s z\" fill=\"%s\" stroke=\"%s\"/>\n",
           $fx, $fy, $rest, $colour, $colour);
  }
}
printf(" </g>\n");
printf("</svg>\n");
