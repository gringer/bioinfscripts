#!/usr/bin/perl

# Author: David Eccles (gringer) 2012 [based on vcf2fq by lh3]

use strict;
use warnings;
use Getopt::Std;

&vcf2fq();
exit;

sub vcf2fq {
  my %opts = (d=>3, D=>100000, Q=>10, L=>100, f=>"");
  getopts('d:D:Q:L:f:', \%opts);
  die(qq/
Usage:   vcfutils.pl vcf2fq [options] <all-site.vcf>

Options: -d INT    minimum depth                   [$opts{d}]
         -D INT    maximum depth                   [$opts{D}]
         -Q INT    min RMS mapQ                    [$opts{Q}]
         -L INT    min INDEL Qual                  [$opts{L}]
         -f FASTA  file with reference sequence(s) [$opts{f}]
\n/) if (@ARGV == 0 && -t STDIN);

  my ($last_chr, $seq, $qual, $last_pos, @gaps);
  my $_Q = $opts{Q};
  my $_L = $opts{L};
  my $_d = $opts{d};
  my $_D = $opts{D};

  my %baseSeqs = ();

  if($opts{f}){
    open(my $fastaFile, "<", $opts{f}) or die("unable to open reference sequence file: ".$opts{f});
    my $seqID = "";
    $seq = "";
    while(<$fastaFile>){
      chomp;
      if(/^>([^ ]+)/){
        $baseSeqs{$seqID} = $seq unless ($seqID eq "");
        $seqID = $1;
        $seq = "";
      } else {
        $seq .= $_;
      }
    }
    $baseSeqs{$seqID} = $seq unless ($seqID eq "");
    close($fastaFile);
  }

  my %het = (AC=>'M', AG=>'R', AT=>'W', CA=>'M', CG=>'S', CT=>'Y',
			 GA=>'R', GC=>'S', GT=>'K', TA=>'W', TC=>'Y', TG=>'K');

  $last_chr = '';
  while (<>) {
	next if (/^#/);
	my @t = split;
	if ($last_chr ne $t[0]) {
	  &v2q_post_process($last_chr, \$seq, \$qual, \@gaps) if ($last_chr);
	  ($last_chr, $last_pos) = ($t[0], 0);
	  $seq = $qual = '';
          if(defined($baseSeqs{$t[0]})){
            $seq .= $baseSeqs{$t[0]};
            $qual .= '~' x length($seq);
          }
	  @gaps = ();
	}
	die("[vcf2fq] unsorted input\n") if ($t[1] - $last_pos < 0);
	if ($t[1] > length($seq)) {
          $seq .=  'n' x ($t[1] - $last_pos);
          $qual .= '!' x ($t[1] - $last_pos);
	}
	if (length($t[3]) == 1 && $t[7] !~ /INDEL/ && $t[4] =~ /^([A-Za-z.])(,[A-Za-z])*$/) { # a SNP or reference
	  my ($ref, $alt) = ($t[3], $1);
	  my ($b, $q);
	  $q = $1 if ($t[7] =~ /FQ=(-?[\d\.]+)/);
	  if ($q < 0) {
		$_ = ($t[7] =~ /AF1=([\d\.]+)/)? $1 : 0;
		$b = ($_ < .5 || $alt eq '.')? $ref : $alt;
		$q = -$q;
	  } else {
		$b = $het{"$ref$alt"};
		$b ||= 'N';
	  }
	  $b = lc($b);
	  $b = uc($b) if (($t[7] =~ /MQ=(\d+)/ && $1 >= $_Q) && ($t[7] =~ /DP=(\d+)/ && $1 >= $_d && $1 <= $_D));
	  $q = int($q + 33 + .499);
	  $q = chr($q <= 126? $q : 126);
	  substr($seq,$t[1]-1,1) = $b;
	  substr($qual,$t[1]-1,1) = $q;
	} elsif (($t[4] ne '.') && ($t[7] =~ /MQ=(\d+)/ && $1 >= $_Q) &&
            ($t[5] > $_L)) {
          # an INDEL
          my $fq = 126;
          if($t[7] =~ /MQ=(\d+)/){
            $fq = $1;
          }
	  push(@gaps, [$t[1], length($t[3]), (split(/,/,$t[4],2))[0], $fq]);
	}
	$last_pos = $t[1];
  }
  &v2q_post_process($last_chr, \$seq, \$qual, \@gaps);
}

sub v2q_post_process {
  my ($chr, $seq, $qual, $gaps) = @_;
  for my $g (@$gaps) {
    # blank out replaced sequence (there was an off-by-1 error in the original code)
    substr($$seq, ($g->[0]-1), $g->[1]) = '-' x ($g->[1]);
    substr($$qual, ($g->[0]-1), $g->[1]) = ' ' x ($g->[1]);
  }
  my $newSeq = '';
  my $newQual = '';
  my $seqPos = 1;
  my $indelOffset = 0;
  for my $g (@$gaps) {
    $newSeq .= substr($$seq,$seqPos-1,($g->[0])-$seqPos);
    $newQual .= substr($$qual,$seqPos-1,($g->[0])-$seqPos);
    $indelOffset = $indelOffset - (($g->[0])-$seqPos);
    $indelOffset = 0 if ($indelOffset < 0);
    $seqPos = $g->[0];
    $newSeq .= lc(substr($g->[2],$indelOffset)) unless ($indelOffset > length($g->[2]));
    my $q = int($g->[3] + 33 + .499);
    $q = chr($q <= 126? $q : 126);
    $newQual .= $q x (length($g->[2]) - $indelOffset) unless ($indelOffset > length($g->[2]));
    $indelOffset = $g->[1] unless ($indelOffset > $g->[1]);
  }
  if($seqPos < length($$seq)){
    $newSeq .= substr($$seq,$seqPos);
    $newQual .= substr($$qual,$seqPos);
  }
  $newSeq =~ tr/\-//d;
  $newQual =~ tr/ //d;
  print "\@$chr\n"; &v2q_print_str(\$newSeq);
  print "+\n"; &v2q_print_str(\$newQual);
}

sub v2q_print_str {
  my ($s) = @_;
  my $l = length($$s);
  for (my $i = 0; $i < $l; $i += 60) {
	print substr($$s, $i, 60), "\n";
  }
}
