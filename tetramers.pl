#!/usr/bin/env perl
use strict;
use warnings;

usage() if $#ARGV != 0;

#Reading data
local $/ = ">";
open(F,'<',"$ARGV[0]") || die $!;
my @data = <F>;
close F;
shift @data;

#Generating tetramers
my @nucl = ("A", "C", "G", "T");
my @tetramers;

foreach my $i (@nucl){
	foreach my $j (@nucl){
		foreach my $k (@nucl){
			push @tetramers,$i.$j.$k.$_ foreach(@nucl)
		}
	}
}

#Printing header
print "#SequenceName";
print map {"\t" . $_} @tetramers;
print "\n";

#Iterating through data
foreach(@data) {
	my @result = ();
	my ($title, $seq) = split /\n/,$_,2;
	$seq =~ s/\n//;
	my $length = length($seq) - 3;
	foreach(@tetramers){
		my $count = () = ($seq =~ /$_/gi);
		push @result,sprintf "%.4f", 100*$count/$length;
	}
	print "$title";
	print map {"\t" . $_} @result;
	print "\n";	
}

sub usage{
	print "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n";
	print "\#Tetranucleotide statistics for (meta)genomic DNA sequences \#\n";
	print "\#\#\#\#\#\#\# Usage: tetramers.pl <input.fasta> > output.csv \#\#\#\#\#\#\n";
	print "\#\#\#\#\#\# The output is tab-separated and includes header \#\#\#\#\#\#\n";
	print "\#\#\#\#\#\#\#\#\#\#\# The values are presented as percents \#\#\#\#\#\#\#\#\#\#\#\#\n";
	print "\#\#\#\#\#\#\#\#\#\#\# Licensed under The MIT License (MIT) \#\#\#\#\#\#\#\#\#\#\#\#\n";
	print "\#\#\#\#\#\# Copyright (c) Aleksei A. Korzhenkov, 2016-2018 \#\#\#\#\#\#\#\n";
	print "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n";
	exit;
}