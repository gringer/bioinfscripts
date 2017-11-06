#!/usr/bin/perl

## sam2mapped -- created mapped sequences (as fasta files) for aligned regions

my $pos = -1;
my $seqName = "";
my $bestFlags = "";
my $bestID = "";
my $bestSeq = "";
my $bestQual = "";
my $bestLine = "";
my $seenCount = 0;

my $output = "fa"; # can be "fa"

sub printFQ {
  my ($id, $seq, $qual) = @_;
  if($id){
    printf("@%s\n%s\n+\n%s\n", $id, $seq, $qual);
  }
}

sub printFA {
  my ($id, $seq, $qual) = @_;
  if($id){
    printf(">%s\n%s\n", $id, $seq);
  }
}

while(<>){
  if(/^@/){
    next;
  }
  my $line = $_;
  chomp;
  my @F = split(/\t/);
  my $refPos = $F[3];
  my $cigar = $F[5];
  $cigar =~ s/[0-9]S$//;
  my $seq = $F[9];
  my $qual = $F[10];
  my $startTrim = 0;
  my $matchLen = 0;
  while($cigar =~ s/^([0-9]+)([MIDNSHP=X])//){
    my $subLen = $1;
    my $op = $2;
    if($op eq "S"){
      $seq = substr($seq, $subLen);
      $qual = substr($seq, $subLen);
    }
    if($op =~ /[M=XI]/){
      $matchLen += $subLen;
    }
  }
  $seq = substr($seq, 0, $matchLen);
  $qual = substr($qual, 0, $matchLen);
  if($seq){
    printFA($F[0], $seq, $qual);
  }
}
if($output eq "fastq"){
  printSeq($bestID, $bestSeq, $bestQual);
} elsif($output eq "sam"){
  print($bestLine);
}
