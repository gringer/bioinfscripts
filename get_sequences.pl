#!/usr/bin/perl

# get_sequences.pl -- retrieves sequences from NCBI, using a file list

# Author: David Eccles (gringer) 2011 <david.eccles@mpi-muenster.mpg.de>

# use like this: 
# ./get_sequences.pl names.txt nuccore > out.fasta

use warnings;
use strict;

use LWP::Simple;
use LWP::UserAgent;
use XML::Parser;
use XML::Simple;
use Cwd;

sub usage(){
    print "usage : get_sequences.pl <accession file> <database>\n";
    print "output: a fasta file with sequence data\n";
    print "\nOther options:\n";
    print " -help       : Show this screen\n";
    print " -v          : increase verbosity of output\n";
}

our @processContent = ();
our $indentLevel = 0;
our %currentEntry = ();
our $currentName = "";
our $currentType = "";
our $headersPrinted = 0; # false

our $verbose = 0;
our $resultsReceived = 0;
our $totalReceived = 0;

our @validHeadings = ();

our $xParse = XML::Parser->new('Style' => 'Stream');
$xParse->setHandlers (
    Start => \&StartTag,
    End   => \&EndTag,
    Char  => \&Characters);
our $xPStream;

sub StartTag {
    my ($e, $name, %lookup) = @_;
    if($name eq 'DocSum'){
        %currentEntry = ();
    }
    if($name eq 'Id'){
        $currentName = 'Id';
    } elsif(($name eq 'Item') && $lookup{'Name'}){
        if(($lookup{'Type'}) && ($lookup{'Type'} eq 'String')){
            $currentType = 'String';
        } else {
            $currentType = '';
        }
        $currentName = $lookup{'Name'};
    } else {
        $currentName = '';
    }
    if($verbose > 1){
        printf(STDERR "%sStarting %s... ", (' ' x $indentLevel), $name);
        foreach my $key (keys(%lookup)){
            printf(STDERR "%s=%s;", $key, $lookup{$key});
        }
        print(STDERR "\n");
    }
    $indentLevel++;
}
    
sub EndTag {
    my ($e, $name) = @_;
    # do something with end tags
    $indentLevel--;
    if($name eq 'DocSum'){
        $resultsReceived++;
        # foreach my $key (keys(%currentEntry)){
        #     print(STDERR "$key=".$currentEntry{$key}.":");
        # }
        # print(STDERR "\n");
        my @headings = ();
        my @values = ();
        if(!@validHeadings){
            @headings = ('Id', 'Gi','NomenclatureSymbol', 
                         'NomenclatureName', 'TaxId', 
                         'ChrAccVer', 'Chromosome', 'ChrStart', 'ChrStop',
                         'Caption', 'Title', 'Extra');
            @validHeadings = 
                grep {defined($currentEntry{$_})?$_:0} @headings;
            # print column headers
            print(join(',',@validHeadings)."\n");
            print("# ".join(',',keys(%currentEntry))."\n");
            $headersPrinted = 1; # true
        }
        @values = map{
            defined($currentEntry{$_})?$currentEntry{$_}:"";
        } @validHeadings;
        printf("%s\n", join(',', @values));
        printf(STDERR "%s\n", join(";", %currentEntry)) if ($verbose > 1);
        if(($verbose > 0) && ($resultsReceived % 1000 == 0)){
            print(STDERR ".");
        }
    }
    if(($name eq 'Item') || ($name eq 'Id')){
        $currentType = '';
        $currentName = '';
    }
    if($verbose > 1){
        printf(STDERR "%sEnding %s\n", (' ' x $indentLevel), $name);
    }
}
    
sub Characters {
    my ($e, $data) = @_;
    if($verbose > 1){
        if($data !~ /^\s*$/){
            printf(STDERR "%sData: %s\n", (' ' x $indentLevel), $data);
        }
    }
    if($currentName eq 'Id'){
        $currentEntry{'Id'} = $data;
        printf(STDERR "%sAdding Id=%s\n", (' ' x $indentLevel),
               $data) if ($verbose > 1);
    } elsif($currentName){
        # remove trailing spaces
        $data =~ s/\s+$//;
        if(($currentType eq 'String') && ($data)){
            # surround string values with quotes
            $data = "\"$data\"";
        }
        printf(STDERR "%sAdding %s=%s\n", (' ' x $indentLevel),
               $currentName, $data) if ($verbose > 1);
        $currentEntry{$currentName} = $data;
    }
}

sub parseText {
    my ($data, $resp, @other) = @_;
    $xPStream->parse_more($data);
}

my $esearchURL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?';
my $efetchURL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?';

my $searchParams = 'usehistory=y&tool=get_sequences.pl';
my $fetchParams = 'retmode=fasta&rettype=text&tool=get_sequences.pl';

my $dbName = 'nuccore';
my $fileExtension = "_$dbName.fa";

my $fetchFasta = 1;

my @files = ();
my @queryNames = ();

my $database = "";

while(@ARGV){
    my $argument = shift(@ARGV);
    if(-f $argument){ # argument is a file that exists
        push(@files, $argument);
    } elsif ($argument eq "-help"){
        usage();
        exit(0);
    } elsif ($argument eq "-v"){
        $verbose++;
        print(STDERR "increasing verbosity\n");
    } else {
        if(($argument =~ /^([a-z_]+)$/i) && !$database){
            $database = $argument;
        } else {
            print(STDERR "Error: unknown option '$argument'\n");
            usage();
            exit(2);
        }
    }
}

@ARGV = @files;

if(!$database){
    print(STDERR "Error: necessary options were not specified. ".
          "Please indicate the database that names are from\n");
    usage();
    exit(2);
}

my $dbquery = 'db='.$database;

while(<>){
    chomp;
    if(!(/^\s*$/)){
        push(@queryNames, $_);
    }
}

if(!@queryNames){
    print(STDERR "Error: no query terms provided\n");
    usage();
    exit(2);
}


# convert to Entrez query
my @formattedTerms = map{$_.'[Accession]'} @queryNames;
my $query = 'term='.join (" OR ", @formattedTerms);

my $ua = LWP::UserAgent->new;

# use equery to get list of IDs
my $shortenedQuery = $query;
$shortenedQuery =~ s/^(.{10}).*(.{10})$/$1...$2/;
printf(STDERR "querying IDs for '%s' in database '%s'... ", 
       $query, $database);
#       $shortenedQuery, $database);
my $searchURL = $esearchURL.join("&",$dbquery,$query,$searchParams);
printf(STDERR "URL[%s]... ", $searchURL) if ($verbose > 1);
my %searchResult = %{XMLin(get($searchURL))};
my %fetchResult = ();
if(exists($searchResult{'WebEnv'})){
    printf(STDERR "done\n");
    # retrieve environment/query for easier retrieval of results
    my $webEnv = 'WebEnv='.$searchResult{'WebEnv'};
    my $queryKey = 'query_key='.$searchResult{'QueryKey'};
    while($totalReceived < $searchResult{'Count'}){
        # use efetch to get actual results
        printf(STDERR "retrieving %d results for '%s' in database '%s' ...",
               ($searchResult{'Count'} - $totalReceived), $shortenedQuery, 
               $database);
        my $resultStart = sprintf("retstart=%d",$totalReceived);
        my $fetchURL = $efetchURL.join("&",$dbquery,$webEnv,$queryKey,
                                           $fetchParams,$resultStart);
        printf(STDERR "URL[%s]...", $fetchURL) if ($verbose > 1);
        my $res = $ua->get($fetchURL);
        if($res->is_success){
            my $content = $res->content;
            if($content =~ />/){
                for (split /^/, $content){
                    if(/^>/){
                        $resultsReceived++;
                    }
                    if(!/^\s*$/){
                        print;
                    }
                }
            }
            printf(STDERR " done (%d results received)\n", $resultsReceived);
        } else {
            printf(STDERR "Error: failed to fetch results\n");
            exit(3);
        }
        if($resultsReceived == 0){
            printf(STDERR "Warning: no more results received, ".
                   "but more expected\n");
            $totalReceived = $searchResult{'Count'};
        }
        $totalReceived += $resultsReceived;
        $resultsReceived = 0;
    }
} elsif (exists($searchResult{'ERROR'})) {
    printf(STDERR "failed to find any results [%s]\n", 
           $searchResult{'ERROR'});
} else {
    printf(STDERR "failed to find any results\n");
}
