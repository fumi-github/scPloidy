#!/bin/env perl

use strict;
use warnings;

my $maxlength = 2000;

my @prev = (); # previous line
my @curr = (); # current line

while (<>) {
    chomp;
    if (scalar(@prev) == 0) {
	@prev = split(/\t/);
	next;
    }
    @curr = split(/\t/);
    if (# not paired reads
	$prev[0] ne $curr[0]) {
	@prev = @curr;
	next;
    }
    if (# Reference sequence name
	$prev[2] eq $curr[2] and
        # CIGAR; skip soft clip, hard clip etc.
	$prev[5] =~ m/^\d+M$/ and
	$curr[5] =~ m/^\d+M$/ and
	abs($prev[8]) <= $maxlength) {
	if ($prev[8] > 0) {
	    print "$prev[2]\t$prev[3]\t" . ($prev[3] + $prev[8] - 1) . "\t$prev[0]\n";
	} elsif ($curr[8] > 0) {
	    print "$curr[2]\t$curr[3]\t" . ($curr[3] + $curr[8] - 1) . "\t$curr[0]\n";
	}
    }
    @prev = ();
}
