#!/usr/bin/perl

# job2tex.pl -- Create job invoice based on comma-separated job file

# Author: David Eccles (gringer), 2010-2015 <programming@gringer.org>

# use diagnostics -verbose;
use strict;
use warnings;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Time::localtime;
use Text::CSV;
use FindBin qw($Bin);
use File::Temp qw(tempfile);

sub usage {
  print("usage: ./job2tex.pl <input file> [options]\n");
  print("\nOther Options:\n");
  print("-help             : Only display this help message\n");
  print("-output           : Output base file name\n");
	print("-advance          : Invoice is prior to job start\n");
  print("-rate <price>     : Alter per-unit job rate\n");
  print("-fixed <price>    : Fixed-price job\n");
  print("-month <number>   : Paid per month (assumes fixed)\n");
  print("-gstinc           : Total amount includes GST\n");
  print("-desc <string>    : Append <string> to job description\n");
  print("-PO <string>      : Include purchase order ID\n");
  print("-date             : Show date on job summary details\n");
  print("-noupdate,-draft  : Don't update invoice number\n");
  print("-amend <n> <date> : Amend invoice #n, produced on date\n");
  print("                    Can also use <date>,<source>;<rate>\n");
  print("-paid             : Add a 'paid' message on the invoice\n");
  print("\n");
}

sub getCurrency {
  my ($date, $currency) = @_;
  open(my $fh, '-|', "perl $Bin/exchange_rate.pl $date $currency") or die $!;
  my ($source,$rate) = split(/;/, <$fh>);
  return ($source,$rate);
}

my $jobInFilename = 0; # false

my $outputBaseName = 0; # false
my $updateInvoiceNum = 1; # true
my $amendedInvoice = 0; # false
my $jobType = "waged";
my $invoicePaid = 0; # false
my $gstRate = 0.15;
my $gstExclusive = 1; # true
my $gstDefined = 0; # false
my $showDate = 0; # false
my $rateKnown = 0; # false
my %templateFields = ();
my %templateRemoveFlags = ();

$templateFields{"unitRate"} = ""; # false
$templateFields{"advanceText"} = "";
$templateFields{"totalUnits"} = 0;
$templateFields{"payType"} = "hour";
$templateFields{"purchaseOrder"} = "";

$templateFields{"jobDesc"} = "Bioinformatics Services";


my @months = qw(   Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
my @monthdays = qw( 31  28  31  30  31  30  31  31  30  31  30  31 );

# extract command line arguments
while(@ARGV){
  my $argument = shift @ARGV;
  if(-f $argument){ # file existence check
    $jobInFilename = $argument;
  } else {
    if($argument eq "-help"){
      usage();
      exit(0);
    }
    elsif($argument eq "-fixed"){
      $jobType = "fixed";
      $templateFields{"unitRate"} = shift @ARGV;
      if($templateFields{"totalUnits"} == 0){
        $templateFields{"totalUnits"} = 1;
      }
    }
    elsif($argument eq "-advance"){
      $templateFields{"advanceText"} = " proposed";
    }
    elsif($argument eq "-paid"){
      $invoicePaid = 1; # true
    }
    elsif($argument eq "-date"){
      $showDate = 1; # true
    }
    elsif($argument eq "-desc"){
      $templateFields{"jobDesc"} .= " [" . (shift @ARGV) ."]";
    }
    elsif($argument eq "-PO"){
      $templateFields{"purchaseOrder"} = shift @ARGV;
    }
    elsif($argument eq "-gstinc"){
      $gstExclusive = 0; # false
      $gstDefined = 1;
    }
    elsif($argument eq "-rate"){
      $templateFields{"unitRate"} = shift @ARGV;
    }
    elsif($argument eq "-month"){
      $jobType = "fixed";
      $templateFields{"payType"} = "month";
      $templateFields{"totalUnits"} = shift @ARGV;
    }
    elsif($argument eq "-output"){
      $outputBaseName = shift @ARGV;
    }
    elsif($argument eq "-amend"){
      $amendedInvoice = 1;
      $templateFields{"invNumber"} = shift(@ARGV);
      ($templateFields{"priorInvDate"}, $rateKnown) =
        split(/,/, shift(@ARGV));
      printf(STDERR "Amending invoice #%d on %s\n",
             $templateFields{"invNumber"},
             $templateFields{"priorInvDate"});
      $updateInvoiceNum = 0; # false
    }
    elsif(($argument eq "-noupdate") || ($argument eq "-draft")){
      $updateInvoiceNum = 0; # false
    }
    else{
      printf(STDERR "Error: Unknown command-line option, '$argument'\n");
      usage();
      exit(4);
    }
  }
}

if(!$jobInFilename){
    print(STDERR "Error: No valid job input file given\n");
    usage();
    exit(1);
}

if(!$templateFields{"unitRate"}){
    print(STDERR "Error: No unit rate specified\n");
    usage();
    exit(1);
}

if(!$outputBaseName){
  if($jobInFilename =~ m/^(.*?)\./){
    # truncate at location of first '.'
    $outputBaseName = $1;
  } else {
    $outputBaseName = $jobInFilename;
  }
}

if(!$templateFields{"purchaseOrder"}){
  $templateRemoveFlags{"PURCORDER"} = 1;
}

if($templateFields{"priorInvDate"}){
  $templateFields{"amendDate"} =
    sprintf("%04d-%03s-%02d",
            localtime->year + 1900, $months[localtime->mon], localtime->mday);
  $templateFields{"invDate"} = $templateFields{"priorInvDate"};
  $templateRemoveFlags{"ORIGINALINVOICE"} = 1;
  $templateRemoveFlags{"DRAFTINVOICE"} = 1;
} else {
  $templateFields{"invDate"} =
    sprintf("%04d-%03s-%02d",
            localtime->year + 1900, $months[localtime->mon], localtime->mday);
  $templateRemoveFlags{"AMENDEDINVOICE"} = 1;
}

if($showDate){
  $templateRemoveFlags{"NODATELINE"} = 1;
} else {
  $templateRemoveFlags{"DATELINE"} = 1;
}

print(STDERR "Reading job file ($jobInFilename)... ");
# Note: this works with both gzipped files and plain text files
my $jobInFile = new IO::Uncompress::Gunzip "$jobInFilename" or
    die "Unable to open $jobInFilename\n";

my $csv = Text::CSV->new({binary => 1, auto_diag => 1});
my @lineData = ();
my %colHeadings = ();
my $firstLine = 1; # true

my %subjobUnits = ();
my %subjobDates = ();
my @subjobOrder = ();
my $subjobID = 0;

while(<$jobInFile>){
  if((/^[^#]/) && (!/^\s+$/)){
    my $units = 0;
    if ($csv->parse($_)) {
      @lineData = $csv->fields();
    } else {
      my $err = $csv->error_input;
      die "Failed to parse line: $err";
    }
    if($firstLine){
      $firstLine = 0; # false
      foreach my $i (0 .. (scalar(@lineData) - 1)){
        $colHeadings{$lineData[$i]} = $i;
      }
    } else {
      my $timeStart = $lineData[$colHeadings{"TimeStart"}];
      my $timeEnd = $lineData[$colHeadings{"TimeEnd"}];
      $timeStart =~ s#^(..)(..)$#$1 + $2 / 60#e;
      $timeEnd =~ s#^(..)(..)$#$1 + $2 / 60#e;
      #           printf(STDERR "%s  - %s\n", $timeEnd, $timeStart);
      if($lineData[$colHeadings{"DateStart"}] eq
         $lineData[$colHeadings{"DateEnd"}]){
        $units = $timeEnd - $timeStart;
      } else {
        my ($yearS, $monthS, $dayS) =
          split(/-/, $lineData[$colHeadings{"DateStart"}]);
        my ($yearE, $monthE, $dayE) =
          split(/-/, $lineData[$colHeadings{"DateEnd"}]);
        if(($yearS eq $yearE) && ($monthS eq $monthE)){
          $units = $timeEnd - $timeStart + ($dayE - $dayS) * 24;
        } else {
          print(STDERR "Start and end dates for '".
                $lineData[$colHeadings{"Description"}].
                "' do not match. " .
                "Units will be set to 0 for this subjob.");
        }
      }
      if($units < 0){
        print(STDERR "Calculated cost is negative for '".
              $lineData[$colHeadings{"Description"}].
              "'. " .
              "Units will be set to 0 for this subjob.\n");
        $units = 0;
      }
      if($units == 0){
        print(STDERR "Calculated cost is 0 hours for '".
              $lineData[$colHeadings{"Description"}].
              "'. " .
              "Units will remain at $units for this subjob.\n");
      }
      if($units > 24){
        print(STDERR "Calculated cost is greater than 24 hours for '".
              $lineData[$colHeadings{"Description"}].
              "'. " .
              "Units will remain at $units for this subjob.\n");
      }
      if(defined($subjobUnits{$colHeadings{"Description"}})){
        printf(STDERR "WARNING: Description '%s' already exists ".
               "and will be replaced\n", $colHeadings{"Description"});
      }
      $lineData[$colHeadings{"Description"}] .=
        sprintf('_%02d',$subjobID);
      push(@subjobOrder, $lineData[$colHeadings{"Description"}]);
      $subjobUnits{$lineData[$colHeadings{"Description"}]} = $units;
      $subjobDates{$lineData[$colHeadings{"Description"}]} =
        $lineData[$colHeadings{"DateStart"}];
      if($jobType ne "fixed"){
        $templateFields{"totalUnits"} += $units;
      }
      $subjobID++;
    }
  }
}
close($jobInFile);
print(STDERR "done!\n");

$templateFields{"logoLocation"} =
  "/bioinf/Work/images/gringene_8_logo_imageonly.pdf";
if(!$templateFields{"invNumber"}){
  my $invoicePath = "../invnumber.txt";
  $templateFields{"invNumber"} = qx{cat ${invoicePath}};
  chomp($templateFields{"invNumber"});
  $templateFields{"invNumber"}++;
  if($updateInvoiceNum){
    open(INVOICEFILE, "> ${invoicePath}")
      or die("Cannot open invoice number file for writing");
    print(INVOICEFILE ($templateFields{"invNumber"})."\n");
    close(INVOICEFILE);
  } else {
    $templateFields{"invNumber"} = "DRAFT";
  }
}
my $templatePath = $Bin."/jobTemplate.tex";

$templateFields{"jobID"} = "";
my $codeName = "";
$templateFields{"clientBusinessName"} = "";
$templateFields{"clientAddress"} = "";
$templateFields{"accountDetails"} = <<EOT;
Name:   D A Eccles \\& J A D Eccles
Bank:   TSB Bank
Branch: TSB Bank Direct, New Plymouth
Number: 15-3959-0522435-00
EOT
$templateFields{"GSTNumber"} = "73-199-131"; # GST registration number

$templateFields{"currency"} = "NZD";
my $country = "nz";

my $detailsHeader = "\\textbf{Sub-job Details} & \\textbf{Hours}\\\\";

$templateFields{"bankMessage"} = "Please transfer payment to the following bank account:";

if(!$templateFields{"invNumber"}){
  my $invoicePath = "../invnumber.txt";
  $templateFields{"invNumber"} = qx{cat ${invoicePath}};
  chomp($templateFields{"invNumber"});
  $templateFields{"invNumber"}++;
  if($updateInvoiceNum){
    open(INVOICEFILE, "> ${invoicePath}")
      or die("Cannot open invoice number file for writing");
    print(INVOICEFILE ($templateFields{"invNumber"})."\n");
    close(INVOICEFILE);
  } else {
    $templateFields{"invNumber"} = "DRAFT";
  }
}

my %fieldLookup = ();
my %clientLookup = ();
{
  my $clientLookupPath = "clientlookup.csv";
  open(my $clientFile, "< ${clientLookupPath}");
  while(my $row = $csv->getline($clientFile)){
    if($row->[0] eq "clientCode"){
      my $fieldIndex = 0;
      foreach my $field (@{$row}){
        $fieldLookup{$field} = $fieldIndex++;
      }
    } else {
      grep {s/\\n/\n/g} @{$row};
      #print(STDERR "Adding ".join(";",@{$row})."\n");
      foreach my $field (keys(%fieldLookup)){
        my $data = $row->[$fieldLookup{$field}];
        $clientLookup{$row->[$fieldLookup{clientCode}]}{$field} =
          $data;
      }
    }
  }
  close($clientFile);
}

my %codeLookup = ();
{
  my $codeLookupPath = "codelookup.csv";
  my %fieldLookup = ();
  open(my $codeFile, "< ${codeLookupPath}");
  while(my $row = $csv->getline($codeFile)){
    if($row->[0] ne "clientCode"){
      $codeLookup{$row->[0]}{$row->[1]} = $row->[2];
    }
  }
  close($codeFile);
}

if($outputBaseName =~
   m#(([A-Z]+)-([0-9]{6}|20[0-9]{2}-[A-Z][a-z]{2}-[0-9]{2})-([A-Z]+))/#){
  printf(STDERR "Found an appropriate job code\n");
  $templateFields{"jobID"} = $1;
  my $clientCode = $2;
  $templateFields{"clientBusinessName"} = $2;
  $codeName = $4;
  if(defined($clientLookup{$clientCode})){
    foreach my $field (keys(%fieldLookup)){
      if($clientLookup{$clientCode}{$field}){
        $templateFields{$field} = $clientLookup{$clientCode}{$field};
      }
    }
  }
  if(defined($codeLookup{$clientCode}{$codeName})){
    $templateFields{clientName} = $codeLookup{$clientCode}{$codeName};
  }
  if ($country eq "de") {
    $templateFields{"currency"} = "EUR";
    $templateFields{"accountDetails"} = <<EOT;
Name:         David Eccles \\& Jessica Eccles
Bank:         Sparkasse M\\"unsterland Ost
IBAN:         DE35 4005 0150 0135 6447 14
BIC:          WELADED1MST
Bankleitzahl: 400 501 50
Kontonummer:  0135644714
EOT
  }
  if($country eq "nz"){
    printf(STDERR "Country: New Zealand\n");
    $templateFields{"currency"} = "NZD";
    if(!$gstDefined){
      $gstExclusive = 1; # true
      $gstDefined = 1;
    }
  }
  if($country eq "us"){
    printf(STDERR "Country: United States\n");
    $templateFields{"currency"} = "USD";
    if(!$gstDefined){
      $gstExclusive = 0; # false
      $gstDefined = 1;
    }
  }
  if($country eq "au"){
    printf(STDERR "Country: Australia\n");
    if(!$gstDefined){
      $gstExclusive = 1; # true
      $gstDefined = 1;
    }
    $templateFields{"currency"} = "AUD";
    $templateFields{"bankMessage"} = "Please use the following details ".
      "for international money transfer:";
    $templateFields{"accountDetails"} = <<EOT;
Bank:   HSBC Bank
Branch: 9th Floor, 1 Queen St, Auckland NZ
SWIFT:  HSBCNZ2A

*For futher credit to*
Bank:   TSB Bank NZ Ltd, Bank Direct
Branch: 120 Devon Street East, New Plymouth

*For credit to*
Name:   D A Eccles \\& J A D Eccles
Bank:   TSB Bank
Branch: TSB Bank Direct, New Plymouth
Number: 15-3959-0522435-00
EOT
  }

  if(!$gstDefined){
    printf(STDERR "Value is exclusive of GST\n");
    $gstExclusive = 1;
    $gstDefined = 1;
  }

  $templateFields{"clientAddress"} =~ s/\n/\\\\\n/g;
  $templateFields{"accountDetails"} =~ s/^\s*(.+)$/\\textbf{$1/mg;
  $templateFields{"accountDetails"} =~ s/^\\textbf\{\*/\\emph{/mg;
  $templateFields{"accountDetails"} =~ s/:/} & /g;
  $templateFields{"accountDetails"} =~ s/\*/:} & /g;
  $templateFields{"accountDetails"} =~ s/\n/\\\\\n/g;
  my $suff = ($amendedInvoice)?"_amended":"";
  $suff .= (!$updateInvoiceNum && !$amendedInvoice)?"_draft":"";
  my $outFileName =
    sprintf("%s_invoice_%s_%s%s.tex",
            $outputBaseName, $templateFields{"invNumber"},
            $templateFields{"invDate"}, $suff);
  my $outPDFName =
    sprintf("%s_invoice_%s_%s%s.pdf",
            $outputBaseName, $templateFields{"invNumber"},
            $templateFields{"invDate"}, $suff);
  if($jobType eq "fixed"){
    $detailsHeader = "\\textbf{Sub-job Details}\\\\";
  }
  $templateFields{"jobLines"} = "";
  foreach my $subJob (@subjobOrder){
    my $subJobDesc = $subJob;
    $subJobDesc =~ s/_[0-9]+$//;
    $subJobDesc =~ s/&/\\&/g;
    $subJobDesc =~ s/#/\\#/g;
    $subJobDesc =~ s/_/\\_/g;
    $subJobDesc =~ s/\\n/\n/g;
    if($jobType ne "fixed"){
      if(!$showDate){
        $templateFields{"jobLines"} .=
          "\\multicolumn{2}{l}{$subJobDesc} & ".
            sprintf("%-1.2f", $subjobUnits{$subJob})."\\\\\n";
      } else {
        $templateFields{"jobLines"} .=
          "$subJobDesc & ".
          $subjobDates{$subJob} . " & ".
            sprintf("%-1.2f", $subjobUnits{$subJob})."\\\\\n";
      }
    } else {
      $templateFields{"jobLines"} .=
        "\\multicolumn{4}{l}{$subJobDesc}\\\\"
    }
  }
  $templateFields{"tAmt"} =
    sprintf("%-1.2f",$templateFields{"totalUnits"} *
            $templateFields{"unitRate"});
  if(($jobType eq "fixed") && ($templateFields{"totalUnits"} == 1)){
    $templateRemoveFlags{"UNITLINE"} = 1;
  } else {
    $templateRemoveFlags{"FIXEDLINE"} = 1;
  }
  if($templateFields{"payType"} ne "month"){
    $templateFields{"totalUnits"} =
      sprintf("%-1.2f",$templateFields{"totalUnits"});
  }
  if($gstExclusive){
    printf(STDERR "Calculating additional GST\n");
    my $gstComponent = $templateFields{"tAmt"} * $gstRate;
    my $amountDue = $templateFields{"tAmt"} +
      $templateFields{"tAmt"} * $gstRate;
    $templateFields{"GSTAmt"} = sprintf("%0.2f", $gstComponent);
    $templateFields{"dAmt"} = sprintf("%0.2f", $amountDue);
    $templateRemoveFlags{"GSTINC"} = 1;
  } else {
    printf(STDERR "Calculating GST component\n");
    my $gstComponent = ($templateFields{"tAmt"} * 3) / 23;
    $templateFields{"GSTAmt"} = sprintf("%0.2f", $gstComponent);
    $templateRemoveFlags{"GSTADD"} = 1;
    $templateRemoveFlags{"DUEAMOUNT"} = 1;
    $templateRemoveFlags{DUENZDLINE} = 1;
  }

  if($templateFields{"currency"} ne "NZD"){
    $templateFields{"exDate"} = $templateFields{"invDate"};
    if($rateKnown){
      ($templateFields{"exSource"},$templateFields{"exRate"}) =
        split(';', $rateKnown);
      printf(STDERR "Parsing known rate: %s [from %s]\n",
            $templateFields{"exRate"},$templateFields{"exSource"});
    } else {
      ($templateFields{"exSource"},$templateFields{"exRate"}) =
        getCurrency($templateFields{"exDate"}, $templateFields{"currency"});
    }
    $templateFields{"tNZDAmt"} =
      sprintf("%0.2f",$templateFields{"tAmt"} / $templateFields{"exRate"});
    $templateFields{"GSTNZDAmt"} =
      sprintf("%0.2f",$templateFields{"GSTAmt"} / $templateFields{"exRate"});
    if(exists($templateFields{"dAmt"})){
      $templateFields{"dNZDAmt"} =
        sprintf("%0.2f",$templateFields{"dAmt"} / $templateFields{"exRate"});
    }
  } else {
    $templateRemoveFlags{GSTNZDLINE} = 1;
    $templateRemoveFlags{DUENZDLINE} = 1;
    $templateRemoveFlags{TOTNZDLINE} = 1;
    $templateRemoveFlags{NZDLINE} = 1;
  }

  if(!$invoicePaid){
    $templateRemoveFlags{"PAIDINVOICE"} = 1;
  } else {
    $templateRemoveFlags{"DUEINVOICE"} = 1;
  }

  if(!$templateFields{"clientName"}){
    $templateRemoveFlags{"CLIENTNAME"} = 1;
  }

  if(!$updateInvoiceNum){
    $templateRemoveFlags{"ORIGINALINVOICE"} = 1;
    $templateRemoveFlags{"INVNUM"} = 1;
  } else {
    $templateRemoveFlags{"DRAFTINVOICE"} = 1;
  }

  print (STDERR "Creating invoice file ($outFileName)... ");
  open(TEMPLATEFILE, "< ${templatePath}")
    or die ("Cannot open template invoice file for reading");
  open(TEXFILE, "> $outFileName");
  while(<TEMPLATEFILE>){
    if(/%%([A-Z]+)$/){
      if(exists($templateRemoveFlags{$1})){
        $_ = "";
      }
    }
    while(/@@@(.*?)@@@/){
      if(exists($templateFields{$1})){
        s/@@@(.*?)@@@/$templateFields{$1}/e;
      } else {
        printf(STDERR "Warning: template field '%s' not found\n",
               $1);
        s/(@@@)(.*?)(@@@)/@@ $2 @@/;
      }
    }
    print(TEXFILE $_);
  }
  close(TEXFILE);
  print(STDERR "done!\n");
  print(STDERR "Compiling invoice file to PDF file...\n");
  my $outDir = $outFileName;
  $outDir =~ s#/[^/]*$##;
  $outFileName =~ s#^.*/([^/]*)$#$1#;
  my @args = ("-output-directory",$outDir,
              "-interaction","batchmode",
              $outFileName);
  my $retVal = system("pdflatex", "-draftmode",
                      @args, "2>&1");
  if($retVal == 0){
    # run twice again to correct indexes
    system("pdflatex",@args);
    system("pdflatex",@args);
    printf(STDERR "done, '(%s)' should now exist!\n", $outPDFName);
  } else {
    printf(STDERR "**** Error running script to generate file: %s ****\n",
           $outPDFName);
    print(STDERR "Command run was:\n");
    print(STDERR "pdflatex ".join(" ",@args)."\n");
    print(STDERR "\n##############################\n\n");
    @args = ("-output-directory",$outDir,
                "-interaction","nonstopmode",
                $outFileName);
    my $retVal = system("pdflatex", "-draftmode",
                        @args, "2>&1");
    print(STDERR "\n##############################\n\n");
  }
} else {
  printf(STDERR "Error: output base name (%s) does not match ".
         "expected format\n".
         "(([A-Z]+)-([0-9]{6}|20[0-9]{2}-[A-Z][a-z]{2}-[0-9]{2})-([A-Z]+))\n",
         $outputBaseName);
}
