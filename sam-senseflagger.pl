#!/usr/bin/perl

## sam-sense_flagger.pl -- re-flag a SAM file using sense and
##   anti-sense mapping information to create a strand-specific
##   output. Query-reverse sequences are set to read #1, while
##   Query-forward are set to read #2

use Getopt::Long qw(:config auto_help pass_through);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my $fwdFileName = "";
my $revFileName = "";
my $prefix = "split";
my $filter = 0; ## false
my $filterFwd = 0; ## false
my $filterRev = 0; ## false

GetOptions("fwdids=s" => \$fwdFileName, "revids=s" => \$revFileName,
           "filter!" => \$filter, "xfwd!" => \$filterFwd,
           "xrev!" => \$filterRev, ) or
  die("Error in command line arguments");

if(!$fwdFileName || !$revFileName){
  die("Error: 'fwdFileName' and 'revFileName' have not been defined");
}

my %fwdIDs;
my %revIDs;

open(my $fwdFile, "<", $fwdFileName) or die("Error loading foward ID file");
while(<$fwdFile>){
  chomp;
  $fwdIDs{$_} = $_;
}
close($fwdFile);

open(my $revFile, "<", $revFileName) or die("Error loading reverse ID file");
while(<$revFile>){
  chomp;
  $revIDs{$_} = $_;
}
close($fwdFile);

my $pos = -1;
my $seqName = "";
my $bestFlags = "";
my $bestID = "";
my $bestSeq = "";
my $bestQual = "";
my $bestLine = "";
my $seenCount = 0;

sub printSeq {
  my ($id, $seq, $qual) = @_;
  if($id){
    printf("@%s\n%s\n+\n%s\n", $id, $seq, $qual);
  }
}

while(<>){
  if(/^@/){
    print;
    next;
  }
  chomp;
  my @F = split(/\t/);
  my $flag = $F[1];
  my $revMap = $flag & 0x10; # mapping orientation to *reference*
  my $sense = 0;
  my $antisense = 0;
  my ($fwdRead, $revRead) = (0, 0); # false
  $flag = $flag &= ~0xC3; # clear 0x80, 0x40, 0x01, and 0x03 flags
  if($revIDs{$F[0]}){ ## reverse-mapped to expected *query*
    $flag |= 0x43;
    $fwdRead = $revMap;
    $revRead = !$revMap;
  } elsif($fwdIDs{$F[0]}){ ## forward-mapped to expected *query*
    $flag |= 0x83;
    $fwdRead = !$revMap;
    $revRead = $revMap;
  } else{
    $flag |= 0x200;
  }
  if(!$filter | !(($flag & 0x200) || ($flag & 0x100) || ($flag & 0x800))){
    $F[1] = $flag;
    if(($filterFwd && $fwdRead) || ($filterRev && $revRead) ||
       !($filterFwd | $filterRev)){
      print(join("\t",@F)."\n");
    }
  }
}
