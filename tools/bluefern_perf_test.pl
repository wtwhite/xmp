#!/usr/bin/perl

use strict;
use warnings;
use FindBin;
use lib $FindBin::Bin;
use FastdnampPerfTesting ':all';

my $FASTDNAMP = "../build-bluegene/fastdnamp";
my $mpichanneltype = "bgl";

sub get_class_for($) {
	my ($nodes) = @_;

	if ($nodes <= 32) {
		return "bg32_100";
	} elsif ($nodes <= 128) {
		return "bg128_100";
	} elsif ($nodes <= 512) {
		return "bg512_100";
	} else {
		die;
	}
}

for (my $nodes = 1; $nodes <= 512; $nodes *= 2) {
	foreach my $lb ("-Bp", "-Bd") {
		foreach my $fn (@files) {
			my $bn = $fn;
			$bn =~ s/\.phy$//i;
			
			my $ncpus = $nodes * 2;			# Cos we use VN mode
			my $class = get_class_for($nodes);
			
			next if $ncpus < 64 && $bn eq '53humans.32';	# Need lotsa CPUs for this one
			if (!boundFor($bn, 1) || $lb eq '-Bp') {		# Don't try -Bd on huge datasets...
				my $obn = "$bn.n$ncpus.$mpichanneltype.$lb";
				
				my $bound = boundFor($bn);
				
				my $x = <<'THE_END';
#!/bin/sh
# @ output = $(Executable).$(Cluster).out
# @ error = $(Executable).$(Cluster).err
# @ wall_clock_limit = 03:00:00
# @ notification = complete
# @ notify_user = w.t.white@massey.ac.nz
# @ job_type = bluegene
# Set your shell type if it's not ksh or bash:
# @ shell = /bin/sh
# @ group = bg01
THE_END
				
				$x .= <<THE_END;
# @ class = $class
# @ bg_size = $nodes
# @ queue
mpirun -np $ncpus -verbose 1 -mode VN -env "TMPDIR=/hpc/bluefern/wtw13/tmp" -cwd /hpc/home/wtw13/parallel_parsimony/measure "$FASTDNAMP" -o "$obn.tre" $lb $bound -m 100000 --seed=1 --displaynewbounds=n --displaynewtrees=n --updatesecs=0 "$fn"
THE_END
				
				my $fname = "submit.$obn.sh";
				print STDERR "Writing $fname...\n";
				open my $f, '>', $fname or die;
				print $f $x;
				close $f;
			}
		}
	}
}
