#!/usr/bin/perl

use warnings;
use strict;

use Math::Trig qw(:pi);

my @pts = ();
my $lastX = undef;
my $lastY = undef;
my $cull = 0; # false
my $zAng = undef;
my $colourVal = undef;

my %lines = ();
my %colours = ();
my $minX = undef;
my $minY = undef;
my $minZ = undef;
my $maxX = undef;
my $maxY = undef;
my $maxZ = undef;


while(<>){
    chomp;
    s/^\s+//;
    if(/^facet normal ([0-9.\-e]+) ([0-9.\-e]+) ([0-9.\-e]+)/){
        my ($nx, $ny, $nz) = ($1, $2, $3);
        $ny = -$ny;
        $zAng = atan2($ny,$nx);
        if((abs($nx) < 0.01) && (abs($ny) < 0.01)){
            $colourVal = sprintf("%02x", 255);
        } else {
            $colourVal = sprintf("%02x",((0.2+(abs($zAng) / pi)*0.7) * 255));
        }
        @pts = ();
        $lastX = undef;
        $lastY = undef;
        $minZ = undef;
        $maxZ = undef;
        if(($nz) < 0.1){
            $cull = 1; # true
        } else {
            $cull = 0;
        }
    }
    if(/^endfacet/){
        if(!$cull){
            my $zPos = ($maxZ + $minZ) / 2;
            if(!exists($lines{$zPos})){
                $lines{$zPos} = ();
                $colours{$zPos} = ();
            }
            push(@{$lines{$zPos}}, sprintf("%s z", join(" ", @pts)));
            push(@{$colours{$zPos}}, $colourVal);
        }
    }
    if(/^vertex (.*$)/){
        my ($x, $y, $z) = split(/ /, $1);
        $y = -$y;
        if(!defined($minX) || ($minX > $x)){
            $minX = $x;
        }
        if(!defined($minY) || ($minY > $y)){
            $minY = $y;
        }
        if(!defined($minZ) || ($minZ > $z)){
            $minZ = $z;
        }
        if(!defined($maxX) || ($maxX < $x)){
            $maxX = $x;
        }
        if(!defined($maxY) || ($maxY < $y)){
            $maxY = $y;
        }
        if(!defined($maxZ) || ($maxZ < $z)){
            $maxZ = $z;
        }
        my $pushStr = sprintf("%s,%s", $x, $y);
        if(defined($lastX)){
            $pushStr = sprintf("%0.5f,%0.5f", $x - $lastX, $y-$lastY);
        }
        push(@pts, $pushStr);
        $lastX = $x;
        $lastY = $y;
    }
}

my $svgWidth = sprintf("%0.0f", ($maxX - $minX) * 1.05);
my $svgHeight = sprintf("%0.0f", ($maxY - $minY) * 1.05);
my $cx = $svgWidth/2;
my $cy = $svgHeight/2;

printf("<svg xmlns:svg=\"http://www.w3.org/2000/svg\"
   xmlns=\"http://www.w3.org/2000/svg\"
   preserveAspectRatio=\"xMidYMid meet\"
   viewBox=\"0 0 $svgWidth $svgHeight\"
   version=\"1.1\">\n");
printf(" <g fill=\"none\" ".
       "stroke-width=\"0.03\" stroke-linejoin=\"round\">\n");
foreach my $zPos (sort {$a <=> $b} (keys(%lines))){
    foreach my $line (@{$lines{$zPos}}){
        my ($fx, $fy, $rest) = split(/[, ]/, $line, 3);
        my $colour = shift(@{$colours{$zPos}});
        $colour = "#$colour$colour$colour";
        $fx += $cx;
        $fy += $cy;
        printf("  <path d=\"m%g,%g %s\" fill=\"%s\" stroke=\"%s\"/>\n",
               $fx, $fy, $rest, $colour, $colour);
    }
}
printf(" </g>\n");
printf("</svg>\n");
