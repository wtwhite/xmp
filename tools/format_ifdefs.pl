#!/usr/bin/perl

while (<>) {
	if (/^#if(|n?def) +(.*?)$/) {
		my ($directive, $cond) = ($1, $2);
		if ($directive eq 'ndef') {
			$cond = "not $cond";
		}
		push @ifdefs, $cond;
		print;
	} elsif (/^#endif\b/) {
		die "Unmatched #endif on line $.!\n" if !@ifdefs;
		my ($cond) = pop @ifdefs;
		print "#endif\t// $cond\n";
	} elsif (/^#else\b/) {
		die "Unmatched #else on line $.!\n" if !@ifdefs;
		my ($cond) = pop @ifdefs;
		$cond = "not $cond";
		$cond =~ s/^not not //;
		print "#else\t// $cond\n";
		push @ifdefs, $cond;		# The #endif at the end of an #if... that contains an #else will repeat the condition stated for the #else.
	} else {
		print;
	}
}
