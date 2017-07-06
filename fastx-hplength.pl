#!/usr/bin/perl
use warnings;
use strict;

## fastx-hplength.pl -- get statistics on homopolymers in a fastq/fasta file

use Getopt::Long qw(:config auto_help pass_through);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $idFileName = "";
my $quiet = 0;
my $base = 33;

GetOptions("idfile=s" => \$idFileName, "quiet!" => \$quiet,
           "base" => \$base ) or
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
  # This little gzip dance makes sure the script can handle both
  # gzip-compressed and uncompressed input, regardless of whether
  # or not it is piped
  my $z = new IO::Uncompress::Gunzip($file, "transparent", 1)
    or die "gunzip failed: $GunzipError\n";
  while(<$z>){
    chomp;
    chomp;
    if (!$inQual) {
      if (/^(>|@)((.+?)( .*?\s*)?)$/) {
        my $newSeqID = $2;
        my $newShortID = $3;
        $baseCount += length($seq);
        while ($seq =~ s/(.)(\1*)//) {
          $hpCounts{"$1$2"}++;
        }
        $seq = "";
        $qual = "";
        $seqID = $newSeqID;
      } elsif (/^\+(.*)$/) {
        $inQual = 1;            # true
        $qualID = $1;
      } else {
        $seq .= $_;
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
while($seq =~ s/(.)(\1*)//){
  $hpCounts{"$1$2"}++;
}

my $cumCount = 0;
foreach my $hpChar (sort {length($a) <=> length($b) || $a cmp $b} (keys(%hpCounts))){
  my $hpCount = $hpCounts{$hpChar};
  my $hpBaseCount = $hpCount * length($hpChar);
  $cumCount += $hpBaseCount;
  printf("%10d %10d ( %6.2f%% / %6.2f%% ) %s\n",
         $hpCount,  $hpBaseCount,
         $hpBaseCount * 100 / $baseCount,
         $cumCount * 100 / $baseCount, $hpChar);
}

printf(STDERR "Total sequence length: %d\n", $baseCount);
