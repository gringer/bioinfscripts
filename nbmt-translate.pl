#!/usr/bin/perl
## Translate mitochondrial sequence of Nb mitochondria
## (mostly because UAG is *not* a stop codon)
use warnings;
use strict;

## codon usage table -- consider proportion of mitochondrial codons

## for x in A C G T; do for y in A C G T; do for z in A C G T; do echo -n "${x}${y}${z} => \"" | perl -pe 'tr/T/U/'; transeq <(echo -e ">1\n${x}${y}${z}") /dev/stdout 2>/dev/null | perl -pe 's/$/  \",/;' | grep -v '^>'; done; done; done

my %transTable = ## 3-letter codes are for translations from mtDNA tRNAs
( AAA => "(K)",
  AAC => "(N)",
  AAG => "(K)",
  AAU => "(N)",
  ACA => "(T)",
  ACC => "(T)",
  ACG => "Arg",
  ACU => "(T)",
  AGA => "(S)",
  AGC => "(S)",
  AGG => "(S)",
  AGU => "(S)",
  AUA => "(M)",
  AUC => "(I)",
  AUG => "(M)",
  AUU => "(I)",
  CAA => "(Q)",
  CAC => "(H)",
  CAG => "(Q)",
  CAU => "Met",
  CCA => "(P)",
  CCC => "(P)",
  CCG => "(P)",
  CCU => "(P)",
  CGA => "(R)",
  CGC => "(R)",
  CGG => "(R)",
  CGU => "(R)",
  CUA => "(L)",
  CUC => "(L)",
  CUG => "(L)",
  CUU => "(L)",
  GAA => "Phe",
  GAC => "(D)",
  GAG => "(E)",
  GAU => "Ile",
  GCA => "Cys",
  GCC => "(A)",
  GCG => "(A)",
  GCU => "(A)",
  GGA => "(G)",
  GGC => "(G)",
  GGG => "(G)",
  GGU => "(G)",
  GUA => "Tyr",
  GUC => "Asp",
  GUG => "His",
  GUU => "Asn",
  UAA => "L/*", ## was * (mtDNA/Leu)
  UAC => "Val",
  UAG => "Leu",
  UAU => "(Y)",
  UCA => "Trp",
  UCC => "Gly",
  UCG => "(S)",
  UCU => "Ser",
  UGA => "(W)",
  UGC => "Ala",
  UGG => "Pro",
  UGU => "Thr",
  UUA => "(L)",
  UUC => "Glu",
  UUG => "Gln",
  UUU => "Lys",

);

my $seq = "";
my $seqID = "";
while(<>){
  chomp;
  if(/^>(.*)$/){
    my $newID = $1;
    if($seq){
      printf(">%s\n", $seqID);
      while($seq =~ s/^(.{1,60})//){
	  my $seqHead = $1;
	  for(my $i=0; $i < length($seqHead); $i+=3){
	      my $codon = substr($seqHead, $i, 3);
	      $codon =~ tr/T/U/;
	      if($transTable{$codon}){
		  print($transTable{$codon});
	      } else {
		  print("???");
	      }
	  }
	  print("\n".("|----:----" x 6)."\n");
	  printf("%s\n", $seqHead);
      }
    }
    $seq = "";
    $seqID = $newID;
  } else {
    $seq .= $_;
  }
}
if($seq){
    printf(">%s\n", $seqID);
    my $pos = 0;
    while($seq =~ s/^(.{1,60})//){
	my $seqHead = $1;
	for(my $i=0; $i < length($seqHead); $i+=3){
	    my $codon = substr($seqHead, $i, 3);
	    $codon =~ tr/T/U/;
	    if($transTable{$codon}){
		print($transTable{$codon});
	    } else {
		print("???");
	    }
	}
	print("\n ");
	for(my $i=1; $i < (length($seqHead)-3); $i+=3){
	    my $codon = substr($seqHead, $i, 3);
	    $codon =~ tr/T/U/;
	    if($transTable{$codon}){
		print($transTable{$codon});
	    } else {
		print("???");
	    }
	}
	print("\n  ");
	for(my $i=2; $i < (length($seqHead)-3); $i+=3){
	    my $codon = substr($seqHead, $i, 3);
	    $codon =~ tr/T/U/;
	    if($transTable{$codon}){
		print($transTable{$codon});
	    } else {
		print("???");
	    }
	}
	print("\n".("|----:----" x 6)."\n");
	printf("%s\n", $seqHead);
	printf("%-10s",$pos);
	for(my $posInc = 10; $posInc < 60; $posInc += 10){
	    printf("%-10s",substr(sprintf("%10s",$pos+$posInc),8,2));
	}
	$pos += 60;
	print("\n\n");
    }
}

## TAGTA-
