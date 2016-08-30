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
    push(@fileArgs, $arg);
  } else {
    %ids{$arg} = 1;
  }
}

@ARGV = @fileArgs;

if($idFileName){
  print(STDERR "Reading in id file... ");
  open(my $idFile, "< $idFileName")
  or die("cannot open $idFileName for reading");
  while(<$idFile>){
    chomp;
    my @data = split(/\s+/, $_);
    foreach my $id (@data){
      $ids{$id} = 1;
    }
  }
  close($idFile);
  print(STDERR "done!\n");
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
    if(scalar(@idOrder) < scalar(keys(%ids))){
      printf(STDERR "Warning: some ID values were not found");
    }
    printf("## <Individual/Column IDs: %s > ##\n", join(" ",@idOrder));
    next;
  }
  my ($marker, @gts) = split(/\s+/);
  printf("%-15s %s", @gts[@colOrder]);
}

