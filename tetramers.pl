#!/usr/bin/env perl
use strict;
use warnings;

usage() if $#ARGV != 0;

# Reading data
local $/ = ">";
open(F,'<',"$ARGV[0]") || die $!;
my @data = <F>;
close F;
shift @data;

# Generating nonredundant tetramers
my @nucl = ("A", "C", "G", "T");
my @tetramers;

foreach my $i (@nucl){
	foreach my $j (@nucl){
		foreach my $k (@nucl){
			foreach(@nucl){
				my $t = $i.$j.$k.$_;
				push @tetramers,$t if (not grep { revcomp($t) eq $_} @tetramers);
			}

		}
	}
}

# Printing header
print "#SequenceName";
print map {"\t" . $_} @tetramers;
print "\n";

# Iterating through the data
foreach(@data) {
	my @result = ();
	my ($title, $seq) = split /\n/,$_,2;
	$seq =~ s/\n//;
	my $length = length($seq) - 3;
	foreach(@tetramers){
		my $count = () = ($seq =~ /$_/gi);
		$_ = revcomp($_);
		my $count_r = () = ($seq =~ /$_/gi);
		push @result,sprintf "%.4f", 100*($count+$count_r)/$length;
	}
	
	print "$title";
	print map {"\t" . $_} @result;
	print "\n";	
}

#
sub revcomp{
	my $str = shift;
	my @out = ();
	push @out,comp($_) foreach(reverse split '',$str);
	return join '',@out;
}

# Return complement nucleotide
sub comp{
	my $i = uc shift;
	return "T" if($i eq "A");
	return "G" if($i eq "C");
	return "C" if($i eq "G");
	return "A" if($i eq "T");
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
