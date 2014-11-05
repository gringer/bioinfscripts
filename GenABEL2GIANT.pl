#!/usr/bin/perl

use warnings;
use strict;

use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::Compress::Gzip qw(gzip $GzipError);

my ($snpName, $mafName, @giantNames) = @ARGV;

my %chipSNPs = ();
my %MAFs = ();

print(STDERR "Reading SNP file...");
my $snpFile = new IO::Uncompress::Gunzip("$snpName") or
  die("Unable to open $snpName\n");
my $read = 0;
while(<$snpFile>){
  chomp;
  $chipSNPs{$_} = 1;
  if($read++ > 10000){
    print(STDERR ".");
    $read = 0;
  }
}
close($snpFile);
print(STDERR " done.\n");

print(STDERR "Reading MAF file...");
my $mafFile = new IO::Uncompress::Gunzip("$mafName") or
  die("Unable to open $mafName\n");
$read = 0;
while(<$mafFile>){
  chomp;
  if(/^Marker/){
    next;
  }
  my ($marker, $sex, $base, $MAF) = split(/,/);
  $MAFs{$sex}{$marker}="$base,$MAF";
  if($read++ > 300000){
    print(STDERR ".");
    $read = 0;
  }
}
close($mafFile);
print(STDERR " done.\n");


foreach my $giantName (@giantNames){
  my $sexCode = "";
  if ($giantName =~ /NI_([fmc]).*?_IMPUTE2/) {
    $sexCode = ($1 eq "c") ? "a" : $1;
  }
  my $giantFile = new IO::Uncompress::Gunzip("$giantName") or
    die("Unable to open $giantName\n");
  print(STDERR "Reading GIANT file [$giantName]...");
  $read = 0;
  my $of = new IO::Compress::Gzip "GIANT_${giantName}" or
    die "Unable to open GIANT_${giantName} for writing\n";
  print($of "chrom,markername,strand,n,effect_allele,".
         "other_allele,eaf,imputation_type,chi2.1df,beta,".
         "se,p\n");
  while (<$giantFile>) {
    chomp;
    tr/\"//d;
    if (substr($_,0,1) eq ",") {
      next;
    }
    my ($markername, $chrom, $pos, $strand, $A1, $A2, $n, $beta, $se,
        $chi2_1df, $P1df, $p, $effAB, $effBB, $chi2_2df, $P2df,
        $effect_allele) = split(/,/);
    my ($maa,$maf) = split(/,/, $MAFs{$sexCode}{$markername});
    if (($maf == 1) || ($maf == 0)) {
      next;                     # skip homozygous / missing SNPs
    }
    my $eaf = $maf;
    my $other_allele = $A1;
    if ($maa ne $effect_allele) {
      $eaf = 1 - $maf;
    }
    $eaf = sprintf("%0.4f", $eaf);
    if ($other_allele eq $effect_allele) {
      $other_allele = $A2;
    }
    my $imputation_type = 4;
    if ($chipSNPs{$markername}) {
      $imputation_type = 0;
    }
    if($chi2_1df ne "NA"){
      $chi2_1df = sprintf("%0.4f", $chi2_1df);
    }
    if($beta ne "NA"){
      $beta = sprintf("%0.4f", $beta);
    }
    if($se ne "NA"){
      $se = sprintf("%0.6f", $se);
    }
    if($p ne "NA"){
      $p = sprintf("%0.4f", $p);
    }
    print($of "$chrom,$markername,$strand,$n,$effect_allele,".
           "$other_allele,$eaf,$imputation_type,$chi2_1df,$beta,".
           "$se,$p\n");
    if($read++ > 400000){
      print(STDERR ".");
      $read = 0;
    }
  }
  close($of);
  close($giantFile);
  print(STDERR " done.\n");
}
