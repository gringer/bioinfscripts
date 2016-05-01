#!/usr/bin/perl

my @lineBuffer = ();
my @ids = ();

while(<>){
  my $line = $_;
  if(/^([a-z])\s(.*?)\s/){
    my ($flag, $id) = ($1, $2);
    if($flag eq "s"){
        push(@ids, $id);
    } elsif($flag eq "a"){ # new alignment
      if((scalar(@ids) != 2) || ($ids[0] ne $ids[1])){
        print(join("",@lineBuffer));
      }
      # if(@ids){
      #   printf("a '%s' '%s' %d %d\n", $ids[0], $ids[1],
      #          scalar(@lineBuffer), scalar(@ids));
      # } else {
      #   print("a\n");
      # }
      @lineBuffer = ();
      @ids = ();
    }
  }
  push(@lineBuffer, $line);
}

if((scalar(@ids) != 2) || ($ids[0] ne $ids[1])){
  print(join("",@lineBuffer));
}
