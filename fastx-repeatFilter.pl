#!/usr/bin/perl

use warnings;
use strict;

use Pod::Usage; ## uses pod documentation in usage code
use Getopt::Long qw(:config auto_version auto_help);
use POSIX; ## for ceil

our $VERSION = "0.1";
our $DEBUG = 0;

=head1 NAME

fastx-repeatFilter.pl -- Repeat Match Algorithm for Long-read Evaluation

=head1 SYNOPSIS

./fastx-repeatFilter.pl <reads.fq> [options]

=head2 Options

=over 2

=item B<-help>

Only display this help message

=item B<-lmin>

Minimum repeat length

=item B<-lmax>

Maximum repeat length

=item B<-length>

Repeat length; sets I<lmin> and I<lmax> to the same value

=item B<-skew>

Maximum skew between repeats for matching lines

=item B<-threshold>

Threshold score for identifying repeat reads

=item B<-fraction>

Fraction of read that must pass threshold score

=item B<-reverse>

Remove repeat-containing reads instead of extracting them

=back

=head1 DESCRIPTION

Identifies (filters) repeat-containing reads in a long-read dataset.

=head1 METHODS

something

=cut

=head2 processSeq(id, seq, qual, lmin, lmax, skew, threshold, fraction, reverse)

Analyses the sequence I<seq> to work out if it is likely to contain a
substantial proportion of repeats of length in the range I<lmin> to
I<lmax>. Depending on the value of I<reverse> and whether a repeat
read was detected, either return an empty string or the sequence.

=cut

sub processSeq{
  my ($id, $seq, $qual, $lmin, $lmax, $skew,
      $threshold, $fraction, $reverse) = @_;
  if(!$id || (length($seq) < $lmin)){
    if($reverse){
      if($qual){
        print("@${id}\n${seq}\n+\n${qual}\n");
      } else {
        print(">${id}\n${seq}\n");
      }
      return;
    }
  }
  foreach my $len ($lmin..$lmax){
    my @scores = ();
    for(my $spos = $skew; ($spos+2*$len+$skew+1) < length($seq); $spos += $len){
      my $maxScore = 0;
      foreach my $ofs (-$skew..$skew){
        my $score=0;
        foreach my $c (0..($len-1)){
          if(substr($seq,$c+$spos,1) eq substr($seq,$spos+$c+$len+$ofs,1)){
            $score++;
          }
        }
        if($score > $maxScore){
          $maxScore = $score;
        }
        # printf("---\n%s\n%s\n--- [%d+%d, %0.2f]\n",
        #        substr($seq,$spos,$len),
        #        substr($seq,$spos+$len+$ofs,$len),
        #        $spos, $ofs, $score);
      }
      push(@scores, $maxScore);
      #    printf("%3d %0.2f\n", $spos, $maxScore / $len);
    }
    @scores = (sort {$b <=> $a} (@scores))[0..($#scores * $fraction)];
    print("${len},$scores[0],$scores[$#scores]\n");
  }
  #printf("Length %d, min score with fraction %0.2f: %0.2f\n",
  #      $len, $fraction, $scores[$#scores]/$len);
}


####################################################
# Command line parsing and verification starts here
####################################################

my $argLine = join(" ",@ARGV);

my $options =
  {
   "length" => 100,
   "lmin" => -1,
   "lmax" => -1,
   "skew" => 3,
   "threshold" => 0.8,
   "fraction" => 0.5,
   "reverse" => 0
};

GetOptions($options,
           'length|l=i',
           'lmin=i',
           'lmax=i',
           'skew=i',
           'threshold=f',
           'fraction=i',
           'reverse|v!',
           'debug!' => \$DEBUG,
) or pod2usage(1);

if($options->{"skew"} >= $options->{"length"}){
  die("Skew must be less than repeat length");
}

if(($options->{lmin} == -1) || ($options->{lmax} == -1)){
  $options->{lmin} = $options->{length};
  $options->{lmax} = $options->{length};
}

my $inQual = 0; # false
my $seqID = "";
my $qualID = "";
my $seq = "";
my $qual = "";
while(<>){
  chomp;
  chomp;
  if(!$inQual){
    if(/^(>|@)((.+?)( .*?\s*)?)$/){
      my $newSeqID = $2;
      my $newShortID = $3;
      if($seqID){
        processSeq($seqID, $seq, $qual,
                   $options->{"lmin"}, $options->{"lmax"},
                   $options->{"skew"},
                   $options->{"threshold"}, $options->{"fraction"},
                   $options->{"reverse"});
      }
      $seq = "";
      $qual = "";
      $seqID = $newSeqID;
    } elsif(/^\+(.*)$/) {
      $inQual = 1; # true
      $qualID = $1;
    } else {
      $seq .= $_;
    }
  } else {
    $qual .= $_;
    if(length($qual) >= length($seq)){
      $inQual = 0; # false
    }
  }
}

if($seqID){
  processSeq($seqID, $seq, $qual,
             $options->{"lmin"}, $options->{"lmax"}, $options->{"skew"},
             $options->{"threshold"}, $options->{"fraction"},
             $options->{"reverse"});
}
