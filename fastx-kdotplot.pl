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
my ($sizeX, $sizeY) = (0,0);
my $subseq = "";
my @region = ();
my $kmerLength = 17; ## number of bases in hash keys
my $blockPicture = 0; ## false
my $hLines = ""; ## horizontal lines

GetOptions("kmer=i" => \$kmerLength, "size=s" => \$size,
	   "hlines=s" => \$hLines,
           "region=s" => \$subseq, "altview!" => \$blockPicture ) or
  die("Error in command line arguments");

if($size =~ /^([0-9]+)x([0-9]+)$/){
  $sizeX = $1;
  $sizeY = $2;
} else {
  $sizeX = $size;
  $sizeY = $size;
}

if($subseq){
  @region = split(/\-/, $subseq);
  print(STDERR $region[0], " - ", $region[1], "\n");
}

my @rlengths = ();
my $inQual = 0; # false
my $seqID = "";
my $qualID = "";
my $seq = "";
my $qual = "";
my $buffer = "";

my $im = new GD::Image($sizeX,$sizeY);

my $white = $im->colorAllocate(255,255,255);
my $black = $im->colorAllocate(0,0,0);
my $red = $im->colorAllocate(0x8b,0,0);
my $green = $im->colorAllocate(0,0xA0,0);
my $blue = $im->colorAllocate(0,0,255);
my $darkRed = $im->colorAllocate(0x6b,0,0);
my $darkGreen = $im->colorAllocate(0,0x80,0);
my $darkBlue = $im->colorAllocate(0,0,0xa0);
my $yellow = $im->colorAllocate(0xa0,0x90,0);
my $magenta = $im->colorAllocate(0x90,0,0xa0);
my $cyan = $im->colorAllocate(0,0xa0,0x90);
my $grey = $im->colorAllocate(0x70,0x70,0x70);

$im->setThickness(3);

while(<>){
  chomp;
  chomp;
  if(!$inQual){
    if(/^(>|@)((.+?)( .*?\s*)?)$/){
      my $newSeqID = $2;
      my $newShortID = $3;
      my $len = length($seq);
      my $logLen = ($len == 0) ? 1 : log($len);
      my $ppb = ($sizeX > $len) ? 1 : $sizeX / $len; # pixels per base
      my $ppl = $sizeY / $logLen; # pixels per log
      my $sseq = "";
      if($seqID && (length($seq) > $kmerLength)){
        if($subseq){
          $seq = substr($seq, $region[0], ($region[1]-$region[0]));
          $len = length($seq);
	  $logLen = ($len == 0) ? 0 : log($len);
          $ppb = ($sizeX > $len) ? 1 : $sizeX / $len;
          $ppl = $sizeY / $logLen;
        }
	if($hLines){
	  my @poss = split(/[,; ]/, $hLines);
	  foreach my $pos (@poss){
	    if($blockPicture){
	      $im->line(0, $sizeY - log($pos) * $ppl, 
			$sizeX, $sizeY - log($pos) * $ppl, $grey);
	    }
	  }
	}
        my $countTotal = 0;
        my $countMax = 0;
	my $dist = 0;
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
	      if($blockPicture){
		if($x != $y){
		  $dist = log(abs($x - $y));
		  $im->setPixel($x * $ppb, $sizeY - $dist * $ppl, $red);
		  $im->setPixel($y * $ppb, $sizeY - $dist * $ppl, $magenta);
		}
	      } else {
		$im->setPixel($x * $ppb, $y * $ppb, $red);
	      }
	    }
          }
	}
	## make sure reverse and reverse complement explicitly overwrite
	foreach my $kmer (keys(%posHash)){
	  my @posList = @{$posHash{$kmer}};
          foreach my $x (@posList){
	    foreach my $y (grep {$_ < $x} (@{$posHash{rc($kmer)}})){
	      if($blockPicture){
		if($x != $y){
		  $dist = log(abs($x - $y));
		  $im->setPixel($x * $ppb, $sizeY - $dist * $ppl, $blue);
		  $im->setPixel($y * $ppb, $sizeY - $dist * $ppl, $cyan);
		}
	      } else {
		$im->setPixel($x * $ppb, $y * $ppb, $blue);
	      }
	    }
	    foreach my $y (grep {$_ > $x} (@{$posHash{rev($kmer)}})){
	      if($blockPicture){
		if($x != $y){
		  $dist = log(abs($x - $y));
		  $im->setPixel($x * $ppb, $sizeY - $dist * $ppl, $green);
		  $im->setPixel($y * $ppb, $sizeY - $dist * $ppl, $yellow);
		}
	      } else {
		$im->setPixel($x * $ppb, $y * $ppb, $green);
	      }
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
my $logLen = ($len == 0) ? 1 : log($len);
my $ppb = ($sizeX > $len) ? 1 : $sizeX / $len;
my $ppl = $sizeY / $logLen;
my $sseq = "";
if($seqID && (length($seq) > $kmerLength)){
  if($subseq){
    $seq = substr($seq, $region[0], ($region[1]-$region[0]));
    $len = length($seq);
    $logLen = ($len == 0) ? 0 : log($len);
    $ppb = ($sizeX > $len) ? 1 : $sizeX / $len;
    $ppl = $sizeY / $logLen;
  }
  if($hLines){
    my @poss = split(/[,; ]/, $hLines);
    foreach my $pos (@poss){
      if($blockPicture){
	$im->line(0, $sizeY - log($pos) * $ppl, 
		  $sizeX, $sizeY - log($pos) * $ppl, $grey);
      }
    }
  }
  my $countTotal = 0;
  my $countMax = 0;
  my $dist = 0;
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
	if($blockPicture){
	  if($x != $y){
	    $dist = log(abs($x - $y));
	    $im->setPixel($x * $ppb, $sizeY - $dist * $ppl, $red);
	    $im->setPixel($y * $ppb, $sizeY - $dist * $ppl, $magenta);
	  }
	} else {
	  $im->setPixel($x * $ppb, $y * $ppb, $red);
	}
      }
    }
  }
  ## make sure reverse and reverse complement explicitly overwrite
  foreach my $kmer (keys(%posHash)){
    my @posList = @{$posHash{$kmer}};
    foreach my $x (@posList){
      foreach my $y (grep {$_ < $x} (@{$posHash{rc($kmer)}})){
	if($blockPicture){
	  if($x != $y){
	    $dist = log(abs($x - $y));
	    $im->setPixel($x * $ppb, $sizeY - $dist * $ppl, $blue);
	    $im->setPixel($y * $ppb, $sizeY - $dist * $ppl, $cyan);
	  }
	} else {
	  $im->setPixel($x * $ppb, $y * $ppb, $blue);
	}
      }
      foreach my $y (grep {$_ > $x} (@{$posHash{rev($kmer)}})){
	if($blockPicture){
	  if($x != $y){
	    $dist = log(abs($x - $y));
	    $im->setPixel($x * $ppb, $sizeY - $dist * $ppl, $green);
	    $im->setPixel($y * $ppb, $sizeY - $dist * $ppl, $yellow);
	  }
	} else {
	  $im->setPixel($x * $ppb, $y * $ppb, $green);
	}
      }
    }
  }
}

# make sure we are writing to a binary stream
binmode STDOUT;

# Convert the image to PNG and print it on standard output
print $im->png;
