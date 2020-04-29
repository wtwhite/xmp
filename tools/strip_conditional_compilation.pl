#!/usr/bin/perl

use strict;
use warnings;

my $token = shift;
my $sense = shift;
my $keep = 1;
my $level = 0;
my $activelevel = -1;

while (<>) {
	if (/^#\s*if(n?def)?\b/) {
		if ($activelevel == -1 && m!^#\s*if(n?)def\s+\Q$token\E\b\s*(//.*)?$!) {
			$keep = ($sense + 0) ^ (($1 eq 'n') + 0);
			$activelevel = $level;
		} else {
			if (/^#\s*if\s+.*\b\Q$token\E\b/) {
				print STDERR "Warning: Token '$token' appears on unparseable #if line on line $. of $ARGV!\n";
			}
			print if $keep;
		}

		++$level;
	} elsif (/^#\s*else\b/) {
		--$level;
		if ($level == $activelevel) {
			$keep = !$keep;
		} else {
			print if $keep;
		}

		++$level;
	} elsif (/^#\s*endif\b/) {
		--$level;
		if ($level == $activelevel) {
			$keep = 1;
			$activelevel = -1;
		} else {
			print if $keep;
		}
	} else {
		print if $keep;
	}
}

die "Final level == $level!" if $level != 0;
