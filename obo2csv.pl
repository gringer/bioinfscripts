#!/usr/bin/perl

use warnings;
use strict;

my $category = "";
my %data = ();

print "id,is_a,name,def\n";

while(<>){
  chomp;
  if(/^\[(.*?)\]/){
    if($category eq "Term"){
      printf("%s,%s,\"%s\",%s\n",
             $data{"id"},
             $data{"is_a"} ? $data{"is_a"} : "",
             $data{"name"} ? $data{"name"} : "",
             $data{"def"} ? $data{"def"} : "");
    }
    %data = ();
    $category = $1;
  } elsif(/^([^:]+): (.*)$/){
    my $field = $1;
    my $value = $2;
    $value =~ s/\[.*?\]//g;
    if($field eq "is_a"){
      $value =~ s/ .*$//;
    }
    if($category eq "Term"){
      $data{$field} = ($data{$field}) ? ($data{$field}.";".$value) : $value;
    }
  }
}
