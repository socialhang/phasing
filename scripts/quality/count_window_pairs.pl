#!/usr/bin/env perl
# Counter used in the source workflow: count complete input lines.
use strict;
use warnings;
my $file = shift;
my %counts;
open my $input, '<', $file or die "$file: $!";
while (my $line = <$input>) {
    chomp $line;
    $counts{$line}++;
}
close $input;
for my $key (keys %counts) {
    print "$key\t$counts{$key}\n";
}
