#!perl -w

use strict;

my $line = shift(@ARGV);
my $tepos = shift(@ARGV) or die;
my $p99 = shift(@ARGV) or die;
my $bam = shift(@ARGV) or die;
my $intermediate = shift(@ARGV) or die;
my $flags = shift(@ARGV) or die;

open(F, "<$flags");
my %flag = ();
while(my $l = <F>){

    chomp $l;

    my @f = split(/\t/, $l);

    $flag{$f[0]} = $f[1];
}
close F;
my $left = 0;
my $leftminus = 0;
my $right = 0;
my $rightplus = 0;
my $rangeleft = "";
my $rangeright = "";
my $sam = "samtools view ";
my $outtemp = "";

open(A, "<$tepos");
open(C, ">>$intermediate");
while(my $in = <A>) {

    chomp $in;
    my @a = split(/\t/, $in);

    $left = $a[2];
    $right = $a[3];
    $leftminus = $a[2] - $p99;
    $rightplus = $a[3] + $p99;

    if($leftminus < 1) {

	$leftminus = 1;
    }

    $rangeleft = $a[0] . ":" . $leftminus . "-" . $left;
    $rangeright = $a[0] . ":" . $right . "-" . $rightplus;

    $outtemp = $a[0] . "." . $a[1] . "." . "left.bam";

    system(qq{$sam $bam $rangeleft -o $outtemp});

    open(B, "<$outtemp");
    
    while(my $lineleft = <B>) {

	chomp $lineleft;
	
	my @b = split(/\t/, $lineleft);
	my @c = split(/:/, $b[0]);

	if($flag{$b[1]} == 0) {
	
	    print C $b[2], "\t", $a[1], "\t", $c[1], "\t", $c[2], "\t", "0", "\n";
	}
    }
    close B;

    system(qq{rm -f $outtemp});

    $outtemp = $a[0] . "." . $a[1] . "." . "right.bam";

    system(qq{$sam $bam $rangeright -o $outtemp});

    open(E, "<$outtemp");

    while(my $lineright = <E>) {

	chomp $lineright;

	my @d = split(/\t/, $lineright);
	my @e = split(/:/, $d[0]);

	if($flag{$d[1]} == 1) {
	    
	    print C $d[2], "\t", $a[1], "\t", $e[1], "\t", $e[2], "\t", "1", "\n";

	}
    }
    close E;

    system(qq{rm -f $outtemp}); 

}
close C;
