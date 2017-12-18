#!/usr/bin/perl
use warnings;
use strict;

## fastx-hplength.pl -- get statistics on homopolymers in a fastq/fasta file

use Getopt::Long qw(:config auto_help pass_through);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $quiet = 0;
my $mode = "CG";

GetOptions("quiet!" => \$quiet, "mode=s" => \$mode) or
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
        if($mode eq "CG"){
          $seq =~ s/CG/XX/gi;
          $seq =~ tr/TtCcGg/AaAaAa/;
          $seq =~ s/XX/CG/gi;
        }
        if($seqID){
          if($qual){
            printf("@%s\n%s\n+\n%s\n", $seqID, $seq, $qual);
          } else {
            $seq =~ s/(.{100})/$1\n/g;
            $seq =~ s/\n$//;
            printf(">%s\n%s\n", $seqID, $seq);
          }
        }
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

$mode = uc($mode);

$baseCount += length($seq);
my $cur = "";
my $cchr = "";
if($mode eq "CG"){
  $seq =~ s/CG/XX/gi;
  $seq =~ tr/TtCcGg/AaAaAa/;
  $seq =~ s/XX/CG/gi;
}
if($seqID){
  if($qual){
    printf("@%s\n%s\n+\n%s\n", $seqID, $seq, $qual);
  } else {
    $seq =~ s/(.{100})/$1\n/g;
    $seq =~ s/\n$//;
    printf(">%s\n%s\n", $seqID, $seq);
  }
}
