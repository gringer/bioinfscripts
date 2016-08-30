#!/usr/bin/perl

# gtsubset.pl -- subsets a simplegt file using an individual list

use warnings;
use strict;
use FileHandle;
use Getopt::Long qw(:config auto_version auto_help pass_through);

my $idFileName = "";
my %ids = ();

GetOptions(
           'ids=s' => \$idFileName,
          );

my @fileArgs = ();

foreach my $arg (@ARGV){
    if(-f $arg){
	$idFileName = $arg;
    } else {
	$ids{$arg} = 1;
    }
}

@ARGV = @fileArgs;

if($idFileName){
  print(STDERR "Reading in id file... ");
  open(my $idFile, "< $idFileName")
  or die("cannot open $idFileName for reading");
  while(<$idFile>){
      chomp;
      s/(\s|,).*$//;
      $ids{$_} = 1;
  }
  close($idFile);
  print(STDERR "done!\n");
}

if(!%ids){
    print(STDERR "Error: no IDs specified. Cannot continue.\n");
    exit(1);
}

my %filteredIDs = ();
my @idOrder = ();
my @colOrder = ();

while(<>){
  chomp;
  if(/^#/){
    s/<Individual.*?://g;
    s/##//g;
    s/>//g;
    s/^\s+//;
    s/\s+$//;
    my @origIDs = split(/\s+/);
    for(my $i = 0; $i <= $#origIDs; $i++){
      if($ids{$origIDs[$i]}){
        $filteredIDs{$origIDs[$i]} = $i;
        push(@idOrder, $origIDs[$i]);
        push(@colOrder, $i);
      }
    }
    if(!(@idOrder)){
	print(STDERR "Error: no specified IDs found. Cannot continue.\n");
	exit(1);
    }
    if(scalar(@idOrder) < scalar(keys(%ids))){
      printf(STDERR "Warning: some ID values were not found\n");
    }
    printf("## <Individual/Column IDs: %s > ##\n", join(" ",@idOrder));
    next;
  }
  my ($marker, @gts) = split(/\s+/);
  printf("%-15s %s\n", $marker, join(" ", @gts[@colOrder]));
}

