#!/usr/bin/perl

use strict;
use warnings;

binmode STDIN;
binmode STDOUT;

while (<>) {
	tr/\r//d;
	if (s/\A\tTREE FASTDNAMP_\d+ = \[&U\] //) {
		print;
	}
}
