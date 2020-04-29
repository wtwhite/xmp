#!/usr/bin/perl

# convert_intel_asm_to_gcc.pl
# ===========================
# Created by: WTJW
# Created on: 16/02/2005

use strict;

my $I2G_PATH = "/home/wtwhite/intel2gas/intel2gas-1.3.3";

my $orig_fname = $ARGV[0] or die "Specify filename!";
#my $int_fname = $orig_fname;
#$int_fname =~ s/\.c$/-gcc-raw.c/ or die;
my $int_fname = "$0.tmp";
#my $new_fname = $orig_fname;
#$new_fname =~ s/\.c$/-gcc-cooked.c/ or die;
my ($basename) = $orig_fname =~ m%([^/]+)\..*?$% or die;		# Used for making unique labels

# First, run the file through the intel2gas tool:
my $cmd = "I2G_DATA=$I2G_PATH/ $I2G_PATH/intel2gas -m -I $orig_fname -o $int_fname";
print STDERR "$cmd\n";
system $cmd;

# Now, go through and make various "manual" (well, slightly less automated)
# changes.
open IN, $int_fname or die;
#open OUT, ">$new_fname" or die;

TOP:
while (<IN>) {
	if (/\b__asm__ \(/) {
		# Entering ASM mode
		print;
		
		my $input_opno = 0;
		my $output_opno = 0;
		my @input_ops = ();
		my @output_ops = ();
		my @lines = ();
		my @clobbered = ();
		my @unclobber = ();
		my %labels = ();
		my $something = 0;
		while (<IN>) {
			if (/^\s*:/) {
				# Entering output operand list.  Should be empty.
				my @end_lines = ();
				push @end_lines, $_;
				
				# Go through all the lines, replacing "%INPUTnn"
				# pseudo-params with actual param numbers.
				foreach my $line (@lines) {
					$line =~ s/%INPUT(\d+)/'%' . ($1 + @output_ops)/ge;
					
					# Replace any labels with their "global" name
					foreach my $label (keys %labels) {
						$line =~ s/\b$label\b/$labels{$label}/g;
					}
					
					print $line;
				}
				
				# List output params
				++$something if @output_ops;
				push @end_lines, join(", ", map(defined($_->[0]) ? "[$_->[0]] \"$_->[1]\" ($_->[2])" : "\"$_->[1]\" ($_->[2])", @output_ops)), "\n";
				
				while (<IN>) {
					push @end_lines, $_;
					if (/^\s*:/) {
						# Entering input operand list.  Should be empty.
						++$something if @input_ops;
						push @end_lines, join(", ", map(defined($_->[0]) ? "[$_->[0]] \"$_->[1]\" ($_->[2])" : "\"$_->[1]\" ($_->[2])", @input_ops)), "\n";
						
						my $clobber_mode = 0;
						while (<IN>) {
							if (s/^\s*://) {	# Yes, destroy this
								$clobber_mode = 1;
							}
							
							if (/^\s*\);$/) {
								push @end_lines, ": ", join(", ", map("\"$_\"", @clobbered)), "\n";
								
								# At this point, we decide whether to print anything at all.
								if ($something) {
									foreach (@end_lines) {
										print;
									}
								}
								
								print;
								next TOP;
							}
							
							if ($clobber_mode) {
								# This line contains a list of regs that intel2gas thinks are clobbered.
								chomp;
								while (s/^\s*"(.*?)",?//) {
									my $reg = $1;
									
									# Check whether we should no longer mark this reg as "clobbered"
									if (!grep $reg eq $_, @unclobber) {
										push @clobbered, $1;
										++$something;
									} else {
										print "//AUTOUNCLOBBERED: $reg\n";
									}
								}
							}
						}
					}
				}
##			} elsif (m%^MISMATCH:\s*"\s*movq\s+(mm\d),\[QWORD PTR \1val\]"%) {
#			} elsif (m%^MISMATCH:\s*"\s*movq\s+(mm\d),\[QWORD PTR (\w+)\]\s*"%) {
#				# We patch up this special case
#				push @lines, "//AUTOPATCHED:\n";
#				push @lines, "\"movq	%INPUT$input_opno,%%$1\\n\\t\"\n";
##				push @input_ops, [ "m", "*((long long *) $1val)" ];
#				push @input_ops, [ "m", "*((long long *) $2)" ];
#				push @clobbered, "%$1";
#				++$input_opno;
#			} elsif (m%^MISMATCH:\s*"\s*movdqa\s+(xmm\d),\[QWORD PTR (\w+)\]\s*"%) {
#				# We patch up this special case
#				push @lines, "//AUTOPATCHED:\n";
#				push @lines, "\"movdqa	%INPUT$input_opno,%%$1\\n\\t\"\n";
##				push @input_ops, [ "m", "*((long long *) $1val)" ];
#				push @input_ops, [ "m", "*((long long *) $2)" ];
#				push @clobbered, "%$1";
#				++$input_opno;
#			} elsif (m%^MISMATCH:\s*"\s*movdqa\s+(xmm\d),(xmm\d)\s*"%) {
#				# We patch up this special case
#				push @lines, "//AUTOPATCHED:\n";
#				push @lines, "\"movdqa	%%$2,%%$1\\n\\t\"\n";
#				push @clobbered, "%$1";
#			} elsif (m%^MISMATCH:\s*"\s*movdqa\s+(xmm\d),\[(\w+)\+(\w+)([-+]\d+)?\]\s*"%) {
#				# We patch up this special case
#				push @lines, "//AUTOPATCHED:\n";
#				push @lines, "\"movdqa	$4(%%$2,%%$3),%%$1\\n\\t\"\n";
#				push @clobbered, "%$1";
#			} elsif (m%^MISMATCH:\s*"\s*movdqa\s+\[(\w+)\+(\w+)([-+]\d+)?\],(xmm\d)\s*"%) {
#				# We patch up this special case
#				push @lines, "//AUTOPATCHED:\n";
#				push @lines, "\"movdqa	%%$4,$3(%%$1,%%$2)\\n\\t\"\n";
#			} elsif (m%^MISMATCH:\s*"\s*pshufd\s+(xmm\d),(xmm\d),(\w+)\s*"%) {
#				# We patch up this special case
#				my ($r1, $r2, $immed) = ($1, $2, $3);
#				$immed =~ s/^([0-9][0-9a-fA-F]*)h$/0x$1/;
#				push @lines, "//AUTOPATCHED:\n";
#				push @lines, "\"pshufd	\$$immed,%%$r2,%%$r1\\n\\t\"\n";
#				push @clobbered, "%$r1";
			} else {
				# First, extract out any {{GCC...}} directives.
				my $should_print = 1;
				#HACK: Avoid acting on {{{GCC...}} directives that have been commented "more than once."  Yeech.
#				while (s%\{\{GCC\s+(.*?)\s*\}\}%%) {
				while (!m%^.*//.*//.*\{\{GCC\s+.*?\s*\}\}% && s%\{\{GCC\s+(.*?)\s*\}\}%%) {
					my $directive = $1;
#					if ($directive =~ /^([IO]) "(.*?)" \((.*)\)$/) {
					if ($directive =~ /^([IO]) (?:\[(\w+)\] )?"(.*?)" \((.*)\)$/) {		# Allow optional "[name]" field
						# This indicates a line that should be replaced with an operand constraint.
						my ($mode, $name, $constraint, $expr) = ($1, $2, $3, $4);
						if ($mode eq "I") {
							push @input_ops, [ $name, $constraint, $expr ];
							++$input_opno;
						} else {
							push @output_ops, [ $name, $constraint, $expr ];
							++$output_opno;
						}
						
						# "Autounclobber" feature: If we introduce a constraint on a
						# general-purpose reg, remove this reg from the list of
						# clobbered regs that intel2gas has prepared.
						while ($constraint =~ s/([abcd])//) {
							push @unclobber, "%$1l", "%$1x", "%e$1x";
						}
						
						push @unclobber, "%si", "%esi" if $constraint =~ /S/;
						push @unclobber, "%di", "%edi" if $constraint =~ /D/;
					} elsif ($directive eq 'DEL') {
						# Usually used with the above directive.
						push @lines, "//AUTOREMOVED:$_";
						$should_print = 0;
					} elsif ($directive =~ /^INS (".*")$/) {
						push @lines, "$1\n";
					} else {
						die "Unrecognised directive '$directive'.";
					}
				}
				
				# Is this a label?  If so, we need to rename it so that it
				# is globally unique.
				if (/^\s*"\s*(\w+):/) {
					my $label = $1;
					$labels{$label} = "__${basename}_L$._$label";		# Including line number should make this label globally unique
					print STDERR "Saw label '$label', will rename to '$labels{$label}'.\n";		#DEBUG
				}
				
				s/%/%%/g;
				push @lines, $_ if $should_print;
			}
		}
	} else {
		print;
	}
}

close IN;
#close OUT;
