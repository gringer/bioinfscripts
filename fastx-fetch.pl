#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long qw(:config auto_help pass_through);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $idFileName = "";
my $quiet = 0;
my $minLength = 1;
my $maxLength = 10 ** 12; # 1 Tbp
my $count = -1;
my $invert = 0; # invert logic
my $trim = 0;

GetOptions("idfile=s" => \$idFileName, "quiet!" => \$quiet,
           "reverse|v!" => \$invert, "trim=i" => \$trim,
           "minLength=i" => \$minLength, "maxLength=i" => \$maxLength,
           "count=i" => \$count, ) or
  die("Error in command line arguments");

my %idsToGet = ();
if($trim){
  $minLength = $minLength + $trim * 2;
  $maxLength = $maxLength + $trim * 2;
}

# unknown commands are treated as identifiers
my @files = ();
while(@ARGV){
  my $arg = shift(@ARGV);
  if(-e $arg){
    push(@files, $arg);
  } else {
    $idsToGet{$arg} = 1;
  }
}
@ARGV = @files;

# use stdin if no files supplied
if(!@ARGV){
  @ARGV = '-' unless (-t STDIN);
}

if($idFileName){
  # read sequence IDs from input file
  if(!$quiet){
    printf(STDERR "Attempting to read from input file ($idFileName)\n");
  }
  my $idFile = new IO::Uncompress::Gunzip "$idFileName" or
    die "Unable to open $idFileName\n";
  while(<$idFile>){
    chomp;
    s/^[>@]//;
    s/\s.*$//;
    $idsToGet{$_} = 1;
  }
  close($idFile);
}

if(!$quiet){
  printf(STDERR "Read %d identifiers\n", scalar(keys(%idsToGet)));
}

if(!$quiet && $invert){
  printf(STDERR "Excluding IDs, rather than selecting\n");
}

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
    chomp;
    chomp;
    if (!$inQual) {
      if (/^(>|@)((.+?)( .*?\s*)?)$/) {
        my $newSeqID = $2;
        my $newShortID = $3;
        if ($seqID && (length($seq) >= $minLength) && (length($seq) <= $maxLength)) {
          if ($trim > 0) {
            $seq = substr($seq, $trim, length($seq)-($trim * 2));
            if ($qual) {
              $qual = substr($qual, $trim, length($qual)-($trim * 2));
            }
          }
          if ($qual) {
            printf("@%s\n%s\n+\n%s\n", $seqID, $seq, $qual);
          } else {
            $seq =~ s/(.{100})/$1\n/g;
            $seq =~ s/\n$//;
            printf(">%s\n%s\n", $seqID, $seq);
          }
          if (--$count == 0) {
            $seqID = "";
            last;
          }
        }
        $seq = "";
        $qual = "";
        if ((!(keys(%idsToGet)) || exists($idsToGet{$newSeqID}) || exists($idsToGet{$newShortID})) xor $invert) {
          $seqID = $newSeqID;
        } else {
          $seqID = "";
        }
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
}

if($seqID && (length($seq) >= $minLength) && (length($seq) <= $maxLength)){
  if($trim > 0){
    $seq = substr($seq, $trim, length($seq)-($trim * 2));
    if($qual){
      $qual = substr($qual, $trim, length($qual)-($trim * 2));
    }
  }
  if($qual){
    printf("@%s\n%s\n+\n%s\n", $seqID, $seq, $qual);
  } else {
    $seq =~ s/(.{100})/$1\n/g;
    $seq =~ s/\n$//;
    printf(">%s\n%s\n", $seqID, $seq);
  }
}
