#!/usr/bin/perl

# simplegt2tped.pl -- converts from a simplegt formatted file into a
# transposed PLINK pedigree (tped) file.

# Author: David Eccles (gringer), 2015 <bioinformatics@gringene.org>

use warnings;
use strict;
use FileHandle;
use Pod::Usage;
use Getopt::Long qw(:config auto_version auto_help pass_through);

=head1 DESCRIPTION

converts from a simplegt formatted file into a transposed PLINK
pedigree (tped) file.

=head1 SYNOPSIS

./simplegt2tped.pl [options] -m <marker location file> <genotype file>

=head2 Options

=over 2

=item B<-phenotype> I<file>

Define phenotype information file (used for sex)

=item B<-marker> I<file>

Define marker location file (used for position)

=item B<-outbase> I<string>

Define base file name for output files

=back

=cut

my $phenoFile = "";

my %genogroup = ();

my $numgts = 0;
my $phenoFileName = "";
my $markerFileName = "";
my $outFileBase = "data";
my %markerInfo = ();
my %sexInfo = ();

my @inputFiles = ();

my $uafChange = 0; # false
my $affChange = 0; # false

my $markerColumn = 0;
my $chrColumn = 1;
my $posColumn = 2;
my $idColumn = 0;
my $sexColumn = 1;

GetOptions(
           'marker=s' => \$markerFileName,
           'phenotype=s' => \$phenoFileName,
           'outbase=s' => \$outFileBase,
          );

if(!$markerFileName){
    print(STDERR "Error: no marker files specified on command line\n");
    pod2usage(2);
}

print(STDERR "Reading in marker file... ");
open(my $markerFile, "< $markerFileName")
  or die("cannot open $markerFileName for reading");
while(<$markerFile>){
  my @data = split(/\s+/, $_);
  if(scalar(@data) < 3){
    next;
  }
  if($data[0] !~ /^#/){
    $markerInfo{$data[$markerColumn]} =
      sprintf("%s %s 0 %s",
              @data[($chrColumn,$markerColumn,$posColumn)]);
  }
}
close($markerFile);
print(STDERR "done!\n");

if($phenoFileName){
  print(STDERR "Reading in phenotype file... ");
  open(my $phenoFile, "< $phenoFileName")
    or die("cannot open $phenoFileName for reading");
  while(<$phenoFile>){
    my @data = split(/\s+/, $_);
    if($data[0] !~ /^#/){
      $sexInfo{$data[$idColumn]} = $data[$sexColumn];
    }
  }
  print(STDERR "done!\n");
}

while(@ARGV){
  my $arg = shift(@ARGV);
  if(-f $arg){
    push(@inputFiles, $arg);
  } else {
    print(STDERR "Unknown argument '$arg'");
    pod2usage(1);
  }
}

@ARGV = @inputFiles;

my @indLabels = ();
my @phenoVals = ();
my $phenoNum = 1;
my $skippedMarkers = 0;
my $includedMarkers = 0;
my $lineCounter = 0;

print(STDERR "Reading in input and writing TPED file (${outFileBase}.tped)...");
open(my $outFile, "> ${outFileBase}.tped");
while (<>){
  s/\s+$//;
  my $line = $_;
  if($lineCounter++ > 5000){
    printf(STDERR ".");
  }
  if($line =~ /^##/){
    ## Determine individual labels
    ## This works even if more than one <ID> region is present in the
    ## header line, as might be the case in a 'join'ed file
    if($line =~ /IDs:\s+(.*?)\s*>/){
      @indLabels = ();
      @phenoVals = ();
    }
    while($line =~ /IDs:\s+(.*?)\s*>/){
      my @lineData = split(/\s+/, $1);
      push(@indLabels, @lineData);
      ## generate an array of (@lineData) copies of $phenoNum
      push(@phenoVals, (($phenoNum) x scalar(@lineData)) );
      $phenoNum++;
      $line =~ s/^.*?>//;
    }
  } else{
    my ($marker, $rest) = split(/\s+/, $line, 2);
    if(!$numgts){
      my @genotypes = split(/\s+/,$rest);
      $numgts = (@genotypes);
    }
    $rest = uc($rest);
    $rest =~ s/\s+/ /g; #replace whitespace with spaces
    $rest =~ s/([ACGTN])([ACGTN])/ $1 $2/g; #put spaces in for genotypes
    if(defined($markerInfo{$marker})){
      printf($outFile "%s %s\n", $markerInfo{$marker}, $rest);
      $includedMarkers++;
    } else {
      $skippedMarkers++;
    }
  }
}
close($outFile);
print(STDERR " done");
if($skippedMarkers){
  printf(STDERR "(%d markers skipped)!\n", $skippedMarkers);
} else {
  print(STDERR "!\n");
}

if(scalar(keys(%markerInfo)) > $includedMarkers){
    printf(STDERR "Warning: map file contains more markers (%d) ".
           "than have been retrieved from simplegt file (%d).\n",
        scalar(keys(%markerInfo)), $includedMarkers);
}

print(STDERR "Creating individual labels... ");
if(@indLabels){
  if(scalar(@indLabels) != $numgts){
    die(sprintf("Number of genotypes (%d) does not match labels in ".
                "header line (%d)", $numgts, scalar(@indLabels)));
  }
} else {
  @indLabels = (1 .. $numgts);
}
my $rest = scalar(@indLabels) - ($affChange + $uafChange);
if($affChange){
  @phenoVals = ((2) x $affChange);
  push(@phenoVals, (1) x $rest);
} elsif($uafChange){
  @phenoVals = ((1) x $uafChange);
  push(@phenoVals, (2) x $rest);
} else {
  @phenoVals = ((0) x $rest);
}
print(STDERR "done!\n");

print(STDERR "Writing TFAM file (${outFileBase}.tfam)..");
open(my $famFile, "> ${outFileBase}.tfam");
print($famFile join("\n",(keys(%sexInfo))[0..100])."\n");
for(my $line = 0; $line < @indLabels; $line++){
  my $indLabel = $indLabels[$line];
  my $sex = $sexInfo{$indLabel} ? $sexInfo{$indLabel} : 0;
  #Family ID,Individual ID,Paternal ID,Maternal ID,Sex,Phenotype
  #Sex: (1=male; 2=female; other=unknown)
  printf($famFile "%s %s 0 0 %s %s\n",
         ($line+1, $indLabel, $sex, $phenoVals[$line]));
}
close($famFile);
print(STDERR "done!\n");
