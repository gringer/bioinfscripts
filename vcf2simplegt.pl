#!/usr/bin/perl

use warnings;
use strict;

## vcf2simplegt.pl -- convert from VCF file to simplegt format

use Getopt::Long qw(:config auto_help pass_through);
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my %colNums = ();

my $idFileName = "";
my @idOrder = ();
my $idsSpecified = 0; # false

my %excludeCols = (
  "#CHROM" => 1, "ID" => 1, "POS" => 1,  "REF" => 1,
  "ALT" => 1, "QUAL" => 1,  "FILTER" => 1,
  "INFO" => 1, "FORMAT" => 1);

GetOptions("idFile=s" => \$idFileName) or
    die("Error in command line arguments");

if($idFileName){
  print(STDERR "Retrieving id names from $idFileName...");
  $idsSpecified = 1; # true
  my $idFile = 0;
  $idFile = new IO::Uncompress::Gunzip "$idFileName" or
      die "Unable to open $idFileName\n";
  while(<$idFile>){
    if(/^\"?(.*?)\"?[\s,]+/){
      my $id = $1;
      push(@idOrder, $id);
    }
  }
  close($idFile);
  print(STDERR scalar(@idOrder)." id names extracted\n");
}

my $nextAlleleNum = 0;
my %alleleNums = ();
my @alleles = ();

while(<>){
  if(/^##/){
    next;
  }
  chomp;
  my @F = split(/\t/, $_);
  if(/(^|\s)ID(\s|$)/){
    my $colNum = 0;
    foreach my $colName (@F){
	$colNums{$colName} = $colNum;
	if(!$idsSpecified && (!$excludeCols{$colName})){
	    push(@idOrder, $colName);
	}
	$colNum++;
    }
    @idOrder = grep {defined($colNums{$_})} @idOrder;
    printf("## <Individual/Column IDs: %s > ##\n",
	join(" ", @idOrder));
    next;
  } elsif(!defined($colNums{"ID"})){
      die("IDs have not been defined / found");
  }
  ## by this time, @idOrder should be populated with column IDs
  ## identify alleles
  @alleles = ($F[$colNums{"REF"}], split(/,/, $F[$colNums{"ALT"}]), ".");
  grep {$_ =~ s/\./$#alleles/ge} @alleles;
  printf("%-15s %s\n", $F[$colNums{"ID"}],
	 join(" ",map {join("",@alleles[split(/\|/, $_)])} @F[@colNums{@idOrder}]));
}
