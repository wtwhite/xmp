#!/usr/bin/perl

use strict;
use warnings;
use FindBin;
use lib $FindBin::Bin;
use FastdnampPerfTesting ':all';

my $FASTDNAMP = "../build-orac/fastdnamp";
my $mpichanneltype = "sock";
my $hostname = "orac";

#"Time to read input:" => 
#"Time to order taxa:");
#"Time to determine initial upper bound:");
#"Time to determine lower bounds:");
#"Time to perform branch and bound search:");
#"Time to write output:");
#"Total time taken:");
#sub grabTimes($) {
#	my ($bn) = @_;
#	
#	my $fn = "$bn.phy.err";		#HACK
#	open my $f, "<", $fn or die "Couldn't open '$fn' for input.";
#	local $_;
#	while (<$f>) {
#		
#	}
#	close $f;
#}

for (my $ncpus = 9; $ncpus >= 2; --$ncpus) {
	foreach my $lb ("-Bp", "-Bd") {
		foreach my $fn (@files) {
			my $bn = $fn;
			$bn =~ s/\.phy$//i;
			if (!boundFor($bn, 1) || $lb eq '-Bp') {		# Don't try -Bd on huge datasets...
				my $obn = "$bn.n$ncpus.$mpichanneltype.$lb";

				my $cmd = "mpiexec -n $ncpus \"$FASTDNAMP\" -o \"$obn.tre\" $lb " . boundFor($bn) . " -m 100000 --seed=1 --displaynewbounds=n --displaynewtrees=n --updatesecs=0 \"$fn\" > \"$obn.out\" 2> \"$obn.err\"";
				print STDERR "$cmd\n";
				
				if (system $cmd) {
					print STDERR "Command failed with exit code $?!\n";
				} else {
					print join(",", "fastdnamp", "fastc64", $ncpus, $mpichanneltype, $hostname, $bn, boundFor($bn, 1)), "\n";
				}
			}
		}
	}
}
