#!/bin/sh
for f in *.c; do perl -lne 'BEGIN { $fn = $ARGV[0]; } $seen = 1 if /^#include .(common|switches)\.h./; if (/^#include .dbgprint\.h./ && !$seen) { print $fn; };' $f; done
