#!/usr/bin/perl

use strict;
use warnings;

my $basename = 'aln';		#HACK
if (@ARGV) {
	$basename = shift;
}

die "Unrecognised arguments." if @ARGV;

$_ = <STDIN>;
chomp;
my ($nTaxa, $nSites) = /(\d+)\s+(\d+)/ or die;

my @data;
while (<STDIN>) {
	tr/\r//d;		# Convert to UNIX line endings
	push @data, $_;
}

die "Read the wrong number of sequence data lines!" if @data != $nTaxa;

for (my $i = 3; $i <= $nTaxa; ++$i) {
	my $fn = "$basename.$i.phy";
	print STDERR "Creating '$fn'...\n";
	open my $f, ">", $fn or die "Couldn't open '$fn': $!";
	print $f "$i $nSites\n";
	for (my $j = 0; $j < $i; ++$j) {
		print $f $data[$j];
	}

	close $f;
}
