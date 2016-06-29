#!/usr/bin/perl

my $pos = -1;
my $seqName = "";
my $bestFlags = "";
my $bestID = "";
my $bestSeq = "";
my $bestQual = "";
my $bestLine = "";

my $output = "sam"; # can be "sam" or "fastq"

sub printSeq {
    my ($id, $seq, $qual) = @_;
    if($id){
	printf("@%s\n%s\n+\n%s\n", $id, $seq, $qual);
    }
}


while(<>){
    if(/^@/){
	print;
	next;
    }
    my $line = $_;
    chomp;
    my @F = split(/\t/);
    if(($F[2] ne $seqName) || ($F[3] != $pos) || (length($bestSeq) < length($F[9]))){
	if(($F[2] ne $seqName) || ($F[3] != $pos)){
	    if($output eq "fastq"){
		printSeq($bestID, $bestSeq, $bestQual);
	    } elsif($output eq "sam"){
		print($bestLine);
	    }
	}
	$seqName = $F[2];
	$pos = $F[3];
	$bestLine = $line;
	$bestID = $F[0];
	$bestFlags = $F[1];
	$bestSeq = $F[9];
	$bestQual = $F[10];
    }
}
if($output eq "fastq"){
    printSeq($bestID, $bestSeq, $bestQual);
} elsif($output eq "sam"){
    print($bestLine);
}
