#!/usr/bin/perl

use warnings;
use strict;

use Pod::Usage;
use Getopt::Long qw(:config auto_version auto_help pass_through);
use Time::localtime;

use HTML::TreeBuilder::XPath;

=head1 DESCRIPTION

Retrieve exchange rates for foreign currencies, either a specified
date in the past (from Oanda), or for today (from TSB Bank)

=head1 SYNOPSIS

exchange_rate.pl [<date>] <currency>

Today's date will be used if the date is not specified

=cut

my %monthNums = ( Jan=>1, Feb=>2, Mar=>3, Apr=>4,
                  May=>5, Jun=>6, Jul=>7, Aug=>8,
                  Sep=>9, Oct=>10, Nov=>11, Dec=>12);
my @numMonths = qw(   Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );

my $todayStr = # today
  sprintf("%04d-%03s-%02d",
          localtime->year + 1900, $numMonths[localtime->mon], localtime->mday);

my %options = (date => $todayStr);

GetOptions(\%options, "date=s") or pod2usage();

if((@ARGV) && ($ARGV[0] =~ /20..-([A-Z][a-z]{2}|\d{2})-\d{2}/)){
  $options{"date"} = shift(@ARGV);
}
if(!@ARGV){
  pod2usage("Currency must be specified");
}

my $currency = shift(@ARGV);

$options{"date"} =~ s/(20..-)(\d{2})(-\d{2})/$1.$numMonths[$2-1].$3/e;

my $todayDate = $options{"date"} eq $todayStr;

##printf("TSB Bank NZ;%s\n","0.9629");
##exit(0);

printf(STDERR "Getting exchange rate information for %s%s\n", $options{"date"},
       ( $todayDate ? " [today]" : ""));

if($todayDate){
  printf(STDERR "[Fetching NZD/%s rate for %s from xe.com]\n",
         $currency, $options{"date"});
  my $tree = HTML::TreeBuilder->
    new_from_url('http://www.xe.com/currencyconverter/convert/?Amount=1&From=NZD&To='.$currency)
    or die("Cannot load URL");
  for my $elt ($tree->findnodes('//span[@class="uccResultUnit"]')){
    printf("xe.com,%0.4f\n", $elt->attr("data-amount"));
  }
} else {
  my $date = $options{"date"};
  $date =~ s/20(..)-(...)-(\d{2})/$1."\/".$monthNums{$2}."\/".$3/e;
  my $url = sprintf("http://www.oanda.com/currency/".
                    "historical-rates-classic?date_fmt=jp&".
                    "date=%s&date1=%s&exch=NZD&expr=%s&".
                    "margin_fixed=0&format=CSV&redirected=1",
                   $date,$date,$currency);
  printf(STDERR "[Fetching NZD/%s rate for %s from oanda.com]\n",
         $currency, $options{"date"});
  printf(STDERR "[$url]\n");
  # print(STDERR "[$url]\n");
  my $tree = HTML::TreeBuilder->new_from_url($url);
  my $rate = $tree->findvalue('//pre');
  $rate =~ s/^.*?,//;
  print("oanda.com;$rate\n");
}
