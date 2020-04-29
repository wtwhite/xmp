#!/usr/bin/perl

# FastdnampPerfTesting.pm
# =======================
# Created by: WTJW
# Created on: 22/03/2010

package FastdnampPerfTesting;

use strict;
use warnings;
use Exporter;

our @ISA = 'Exporter';
our @EXPORT_OK = qw(
	@files
	&boundFor
);
our %EXPORT_TAGS = (all => [ @EXPORT_OK ]);

#HACK: Yes, it's a bit gross populating this this way, but who cares.
our @files;
if (@ARGV) {
	@files = @ARGV;
} else {
	@files = glob("*.phy");
}

print STDERR "Data files to be run:\n", join("\n", @files), "\n";

my %bounds = qw(
	its36			233
	53humans.32		269
);

sub boundFor($;$) {
	my ($bn, $plain) = @_;
	
	if (exists $bounds{$bn}) {
		if ($plain) {
			return $bounds{$bn};
		} else {
			return "-b $bounds{$_[0]}";
		}
	} else {
		if ($plain) {
			return 0;
		} else {
			return "";
		}
	}
}

1;
