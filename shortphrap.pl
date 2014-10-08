#!perl -w

use strict;

my $input = shift(@ARGV) or die;
my $output = shift(@ARGV) or die;
my $qual = shift(@ARGV) or die;

my %hash = (); 
my %hash2 = ();
my %qualhash = ();

open (A, "<$input");

my $name = "";

while(my $line = <A>) {

    chomp $line;

    if($line =~ m/>/) {

	$name = $line;

    }else {

	my @a = split(//, $line);

	my $a = 0;
	my $c = 0;
	my $t = 0;
	my $g = 0;
	my $count = 0;

	foreach my $element (@a) {

	    $count++;

	    if ($element eq 'A') {

		$a++;

	    }
	    if ($element eq 'C') {

		$c++;

	    }
	    if ($element eq 'T') {

		$t++;

	    }
	    if ($element eq 'G') {

		$g++;

	    }
	}

	
	if ((($a / $count) < 0.5) and (($c / $count) < 0.5) and (($g / $count) < 0.5) and (($t / $count) < 0.5)){
	    
	    if (!(exists($hash{$line}))) {

		$hash{$line} = $name; 
		$hash2{$name} = "";
	    
	    }
	}
    }
}

close A;

my $qualfile = $input . ".qual";

open(C, "<$qualfile");

my $next = 0;
my $current = 0;

LOOP:while(my $line2 = <C>) {

    chomp $line2;

    if($line2 =~ m/>/) {
	
	if(exists($hash2{$line2})) {

	    $next = 1;
	    $current = $line2;
	    next LOOP;

	}
	
    }elsif($line2 !~ m/>/) {
	
	if($next == 1) {
	    
	    $qualhash{$current} = $line2;
	    $next = 0;
	    next LOOP;
	}
    }
}

close C;

unlink(qq{$output});
open(B, ">>$output");
unlink(qq{$qual});
open(D, ">>$qual");


while((my $key, my $value) = each(%hash)) {

    print B $value, "\n", $key, "\n";
    print D $value, "\n", $qualhash{$value}, "\n";
}
close B;
close D;
