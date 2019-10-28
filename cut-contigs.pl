#!/usr/bin/perl
use strict;
use warnings;
use integer;

if(not defined $ARGV[0]){
	print "Cutting multifasta files without extra dependencies.\n";
	print "Usage: $0 <input.fasta> [window size] > cutted.fasta\n";
	print "By default window size equals 10000.\n";
	print "\*.\* Aleksei Korzhenkov. 2017-2018 \*.\*\n";
	exit;	
}

my $filename = $ARGV[0];

#Setting default maximum fragment length
my $maxl = 10000;
if(defined $ARGV[1]){
	$maxl = $ARGV[1];
}

#Short contigs removal
my $minl = 1000;

#Reading data
local $/ = '>';
open(FILE,"<$filename") or die;
my @fasta = <FILE>;
close FILE;
shift @fasta;

#Line formating and printing
sub printseq{
	my $line = 60;
	my $inseq = shift;
	while(length($inseq) > $line){
		print(substr($inseq,0,$line,'')."\n");
	}
	print("$inseq\n");
}

foreach (@fasta){
	my ($id, $seq) = split /\n/, $_, 2;
	$seq =~ s/\W//g;
	my $len = length($seq);
	next if $len < $minl;
	if($len < $maxl*2){
		print(">$id\n");
		printseq($seq);
	}else{
		my $localmax = $len / ($len/$maxl);
		my $i = 0;
		while(length($seq) >= 2*$localmax){
			print(">$id"."_[$i]\n");
			printseq(substr($seq,0,$localmax));
			substr($seq,0,$localmax,'');
			$i++;
		}
		print(">$id"."_[$i]\n");
		printseq($seq)
	}
}
