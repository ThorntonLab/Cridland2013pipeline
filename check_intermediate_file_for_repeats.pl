#!perl -w

use strict;

my $input = shift(@ARGV) or die;
my $output = shift(@ARGV) or die;
unlink(qq{$output});

my %hash = ();

open(A, "<$input");
while(my $line = <A>) {

    chomp $line;

    $hash{$line} = 1;

}
close A;
open(B, ">>$output");
while((my $key, my $value) = each(%hash)) {

    print B $key, "\n";

}
close B;
