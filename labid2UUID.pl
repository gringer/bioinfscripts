#!/usr/bin/perl

use warnings;
use strict;

use Text::CSV;

open(my $idLookupFile, "<", "/home/gringer/bioinf/GU-2012-Apr-01-RLMB/Miles/individuals/NI_UUID_Ped_2012-Oct-23.csv") or die("Cannot open lookup file");

my $csv = Text::CSV->new ({ binary => 1, eol => $/ });

my %col = ();
my $headerRow = $csv->getline($idLookupFile);
my $fieldCount = 0;
foreach my $field (@$headerRow){
  $col{$field} = $fieldCount++;
}

my %dataReplacement = ();
my $fileType = "tfam";

while (my $row = $csv->getline($idLookupFile)){
  my @fields = @$row;
  if($fields[$col{"LAB_ID"}] ne "NA"){
    my $labid = $fields[$col{"LAB_ID"}];
    my $uuid = $fields[$col{"UUID"}];
    my $gender = $fields[$col{"Gender"}];
    my $genderVal = ($gender eq "Male")?1:($gender eq "Female")?2:0;
    if(($fileType eq "tfam") || ($fileType eq "ped")){
      $dataReplacement{$labid} = 
          join(" ",1,$uuid,
               ($fields[$col{"patID"}] eq "NA")?0:$fields[$col{"patID"}],
               ($fields[$col{"matID"}] eq "NA")?0:$fields[$col{"matID"}],
               $genderVal);
    } else {
      $dataReplacement{$labid} = $uuid;
    }
  }
}
close($idLookupFile);

while(<>){
  if(($fileType eq "tfam") || ($fileType eq "ped")){
    my $output = $_;
    my @fields = split(/\s+/,$output,6);
    if(exists($dataReplacement{$fields[1]})){
      # try the 'individual' column first
      $output = $dataReplacement{$fields[1]}." ".$fields[5];
    } elsif(exists($dataReplacement{$fields[0]})){
      # try the 'family' column
      $output = $dataReplacement{$fields[0]}." ".$fields[5];
    }
    print($output);
  }
}
