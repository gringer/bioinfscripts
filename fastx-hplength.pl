#!/usr/bin/perl
use warnings;
use strict;

## fastx-hplength.pl -- get statistics on homopolymers in a fastq/fasta file

use Getopt::Long qw(:config auto_help pass_through);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $quiet = 0;

GetOptions("quiet!" => \$quiet) or
  die("Error in command line arguments");

# unknown commands are treated as identifiers
my @files = ();
while(@ARGV){
  my $arg = shift(@ARGV);
  if(-e $arg){
    push(@files, $arg);
  }
}
@ARGV = @files;

# use stdin if no files supplied
if(!@ARGV){
  @ARGV = '-' unless (-t STDIN);
}

my %hpCounts = ();

my $baseCount = 0;
my $inQual = 0; # false
my $seqID = "";
my $qualID = "";
my $seq = "";
my $qual = "";

foreach my $file (@ARGV) {
  # This little gunzip dance makes sure the script can handle both
  # gzip-compressed and uncompressed input, regardless of whether
  # or not it is piped
  my $z = new IO::Uncompress::Gunzip($file, "transparent", 1)
    or die "gunzip failed: $GunzipError\n";
  while(<$z>){
    s/\s+$//; # remove ending whitespace
    if (!$inQual) {
      if (/^(>|@)((.+?)( .*?\s*)?)$/) {
        my $newSeqID = $2;
        my $newShortID = $3;
        $baseCount += length($seq);
        my $cur = "";
        my $cchr = "";
        grep { # collect homopolymers
          if($_ ne $cchr){
            $hpCounts{$cur}++ if($cur);
            $cur = $cchr = $_;
          } else {
            $cur .= $cchr;
          }
        } split(//, $seq);
        $hpCounts{$cur}++ if($cur); # collect remaining homopolymer (if any)
        $seq = "";
        $qual = "";
        $seqID = $newSeqID;
      } elsif (/^\+(.*)$/) {
        $inQual = 1;            # true
        $qualID = $1;
      } else {
        $seq .= uc($_);
      }
    } else {
      $qual .= $_;
      if (length($qual) >= length($seq)) {
        $inQual = 0;            # false
      }
    }
  }
  close($z);
}

$baseCount += length($seq);
my $cur = "";
my $cchr = "";
grep { # collect homopolymers
  if($_ ne $cchr){
    $hpCounts{$cur}++ if($cur);
    $cur = $cchr = $_;
  } else {
    $cur .= $cchr;
  }
} split(//, $seq);
$hpCounts{$cur}++ if($cur); # collect remaining homopolymer (if any)

my $cumCount = 0;
foreach my $hpChar (sort {length($a) <=> length($b) || $a cmp $b} (keys(%hpCounts))){
  my $hpCount = $hpCounts{$hpChar};
  my $hpBaseCount = $hpCount * length($hpChar);
  $cumCount += $hpBaseCount;
  printf("%10d %10d ( %6.2f%% / %6.2f%% ) %10s %s %d\n",
         $hpCount,  $hpBaseCount,
         $hpBaseCount * 100 / $baseCount,
         $cumCount * 100 / $baseCount,
	 (length($hpChar) < 10 ? $hpChar : substr($hpChar,0,1)."........."),
	 (length($hpChar) < 10 ? ":" : "x" ), length($hpChar));
}

printf(STDERR "Total sequence length: %d\n", $baseCount) unless $quiet;
