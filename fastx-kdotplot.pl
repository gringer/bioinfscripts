#!/usr/bin/perl

## kmer dotplot -- creates a self-vs-self dotplot of a sequence using
## kmer sequences, outputs a PNG file to standard out
## example usage: ~/scripts/fastx-kdotplot.pl -s 300 -k 13 out_sseq.fa | \
##                convert png:- -resize 1000x1000 - > repeat_region.png

use warnings;
use strict;

use GD; ## for images
use Getopt::Long qw(:config auto_help pass_through);

sub rc {
  my ($seq) = @_;
  $seq =~ tr/ACGTUYRSWMKDVHBXN-/TGCAARYSWKMHBDVXN-/;
  # work on masked sequences as well
  $seq =~ tr/acgtuyrswmkdvhbxn/tgcaaryswkmhbdvxn/;
  return(scalar(reverse($seq)));
}

sub rev {
  my ($seq) = @_;
  return(scalar(reverse($seq)));
}

my $size = 1024;
my $kmerLength = 17; ## number of bases in hash keys

GetOptions("kmer=i" => \$kmerLength, "size=i" => \$size) or
    die("Error in command line arguments");

my @rlengths = ();
my $inQual = 0; # false
my $seqID = "";
my $qualID = "";
my $seq = "";
my $qual = "";
my $buffer = "";

my $im = new GD::Image($size,$size);
my $white = $im->colorAllocate(255,255,255);
my $black = $im->colorAllocate(0,0,0);
my $red = $im->colorAllocate(255,0,0);
my $green = $im->colorAllocate(0,255,0);
my $blue = $im->colorAllocate(0,0,255);
my $yellow = $im->colorAllocate(255,255,0);
my $magenta = $im->colorAllocate(255,0,255);
my $cyan = $im->colorAllocate(0,255,255);

while(<>){
  chomp;
  chomp;
  if(!$inQual){
    if(/^(>|@)((.+?)( .*?\s*)?)$/){
      my $newSeqID = $2;
      my $newShortID = $3;
      my $len = length($seq);
      my $ppb = ($size > $len) ? 1 : $size / $len;
      my $sseq = "";
      if($seqID && (length($seq) > $kmerLength)){
        my $countTotal = 0;
        my $countMax = 0;
        my $maxKmer = "";
	my @rptCounts = ();
	my %posHash = ();
        my %gapCounts = ();
	for(my $p = 0; ($p + $kmerLength) <= $len; $p++){
          push(@{$posHash{substr($seq, $p, $kmerLength)}}, $p);
        }
	foreach my $kmer (keys(%posHash)){
	  my @posList = @{$posHash{$kmer}};
          foreach my $x (@posList){
            foreach my $y (@posList){
              $im->setPixel($x * $ppb, $y * $ppb, $red);
            }
            foreach my $y (@{rc($kmer)}){
              $im->setPixel($x * $ppb, $y * $ppb, $blue);
            }
            foreach my $y (@{rev($kmer)}){
              $im->setPixel($x * $ppb, $y * $ppb, $green);
            }
          }
	}
      }
      $seq = "";
      $qual = "";
      $buffer = "";
      $seqID = $newSeqID;
    } elsif(/^\+(.*)$/) {
      $inQual = 1; # true
      $qualID = $1;
    } else {
      $seq .= $_;
    }
  } else {
    $qual .= $_;
    if(length($qual) >= length($seq)){
      $inQual = 0; # false
    }
  }
}

my $len = length($seq);
my $ppb = ($size > $len) ? 1 : $size / $len;
if ($seqID && (length($seq) > $kmerLength)) {
  my $countTotal = 0;
  my $countMax = 0;
  my $maxKmer = "";
  my @rptCounts = ();
  my %posHash = ();
  my %gapCounts = ();
  for (my $p = 0; ($p + $kmerLength) <= $len; $p++) {
    push(@{$posHash{substr($seq, $p, $kmerLength)}}, $p);
  }
  foreach my $kmer (keys(%posHash)) {
    my @posList = @{$posHash{$kmer}};
    foreach my $x (@posList) {
      foreach my $y (@posList) {
        $im->setPixel($x * $ppb, $y * $ppb, $red);
      }
      foreach my $y (@{$posHash{rc($kmer)}}) {
        $im->setPixel($x * $ppb, $y * $ppb, $blue);
      }
      foreach my $y (@{$posHash{rev($kmer)}}) {
        $im->setPixel($x * $ppb, $y * $ppb, $green);
      }
    }
  }
}

# make sure we are writing to a binary stream
binmode STDOUT;

# Convert the image to PNG and print it on standard output
print $im->png;
