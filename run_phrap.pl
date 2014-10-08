#!perl -w

use strict;
use warnings;
use Parallel::ForkManager;
my $forks = shift(@ARGV) or die;
my $manager = new Parallel::ForkManager($forks); #the number of forks must be adjusted by the user

#these may need to be adjusted by the user to specify the correct directory
my $phrap = 'phrap';
my $longfasta = 'shortphrap.pl'; 

#run phrap

my @out = glob('*.out');

foreach my $out (@out){

    unlink(qq{$out});

}

my @short = glob('short*');

foreach my $short (@short){

    unlink(qq{$short});

}

my @output = glob('*.output');

foreach my $output (@output) {

    unlink(qq{$output});

}

my @contigs = glob('*.contigs');

foreach my $contigs (@contigs) {

    unlink(qq{$contigs});

}

my @ace = glob('*.ace');

foreach my $ace (@ace) {

    unlink(qq{$ace});

}

my @singlet = glob('*.singlets');

foreach my $singlet (@singlet) {

    unlink(qq{$singlet});

}

my @problems = glob('*.problems');

foreach my $problems (@problems) {

    unlink(qq{$problems});

}

my @log = glob('*.log');

foreach my $log (@log) {

    unlink(qq{$log});

}

print "Running Phrap...\n";

my @fasta2 = glob('*.fasta');

foreach my $fasta2 (@fasta2) {

    $manager -> start and next;

    my $linecount = 0;

    open(FILE, "<$fasta2");

    $linecount++ while <FILE>;

    close FILE;

    if ($linecount < 10000) {

    system(qq{$phrap $fasta2 -vector_bound 0 -forcelevel 10 -minscore 10 -minmatch 10 -new_ace});

    }elsif ($linecount >= 10000) {

        my $newfile = "short." . $fasta2;
        my $newqual = $newfile . ".qual";

        system(qq{perl $longfasta $fasta2 $newfile $newqual});

        system(qq{$phrap $newfile -vector_bound 0 -forcelevel 10 -minscore 10 -minmatch 10 -new_ace});

    }
    $manager -> finish;
}
$manager -> wait_all_children;
print "Done! \n";
