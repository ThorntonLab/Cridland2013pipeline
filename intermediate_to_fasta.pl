#!perl -w

use strict;

my $input = shift(@ARGV) or die;
my $reads = shift(@ARGV) or die;
my $value = shift(@ARGV) or die;
my $output = shift(@ARGV) or die;
my %left = ();
my %right = ();
open(A, "<$input");#need to read in input file and make a hash of where all the files should go
while(my $line = <A>) {

    chomp $line;
    my @a = split(/\t/, $line);
    if($a[4] == 0) {
	my $f = $a[0] . "." . $a[1] . ".left.fasta";
	push @{$left{$a[2] . ":" . $a[3]}}, $f;
	
    }elsif($a[4] == 1) {
	my $f2 = $a[0] . "." . $a[1] . ".right.fasta";
	push @{$right{$a[2] . ":" . $a[3]}}, $f2;
    }else{
	die;
    }
}
close A;

open(B, "gunzip -c $reads |") || die "can't open pipe to $input";
while(my $line2 = <B>) {
    
    chomp $line2;
    $line2 =~ s/\s+/\t/go;
    my @b = split(/\t/, $line2);
    if(exists($left{$b[1] . ":" . $b[2]})) {
	foreach my $e (@{$left{$b[1] . ":" . $b[2]}}) {
	   # print $e, "\n";
	    #fix the qual score
	    my $qual = $b[6];
	    #print $qual, "\n";
	    $qual =~ s/\\\\/\\/go;
	    my $ffile = $output . $e;
	    my $qfile = $output . $e . ".qual";
	    open(C, ">>$ffile");
	    open(D, ">>$qfile");
	    
	    print C ">" . $b[0] . "." . $b[1] . "." . $b[2] . "." . $b[3] . "." . $e, "\n";
	    print C $b[5], "\n";	
	    print D ">" . $b[0] . "." . $b[1] . "." . $b[2] . "." . $b[3] . "." . $e, "\n";
	    my @c = split(//, $qual);
	    foreach my $element (@c ){
		print D (ord($element) - $value), " ";
	    }
	    print D "\n";
	    close C;
	    close D;
	}
    }
    if(exists($right{$b[1] . ":" . $b[2]})) {
	foreach my $e2 (@{$right{$b[1] . ":" . $b[2]}}){
	    #fix the qual score
	    my $qual = $b[6];
	    #print $qual, "\n";
	    $qual =~ s/\\\\/\\/go;
	    my $ffile = $output . $e2;
	    my $qfile = $output . $e2 . ".qual";
	    open(E, ">>$ffile");
	    open(F, ">>$qfile");
	    
	    print E ">" . $b[0] . "." . $b[1] . "." . $b[2] . "." . $b[3] . "." . $e2, "\n";
	    print E $b[5], "\n";
	    
	    print F ">" . $b[0] . "." . $b[1] . "." . $b[2] . "." . $b[3] . "." . $e2, "\n";
	    my @d = split(//, $qual);
	    foreach my $element2 (@d ){
		print F (ord($element2) - $value), " ";
	    }
	    print F "\n";
	    
	    close E;
	    close F;
	}
    }    
}
close B;
