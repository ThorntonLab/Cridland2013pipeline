#!perl -w
use strict;
use warnings;
use Parallel::ForkManager;

##################################################################################
#This pipeline detects Transposable elements by processing BWA output
##################################################################################
#This pipeline requires that te_alignment_processing, samtools and phrap are installed and in the users path
##################################################################################

#values that may need to be changed by the user
##################################################################################
my $highch = 14; #The highest number chromosome for a species - 14 is for Drosophila melanogaster - This is because we identify chromosomes by number - see README file
my $mdist = 1000;
my $filternum = 3; #require this many read pairs to call an event
my $check_num = 0;
##################################################################################

#Command line input
##################################################################################
#my $forks = shift(@ARGB) or die; #the number of processors to use
#my $manager = new Parallel::ForkManager($forks); #the number of forks must be adjusted by the user - this is primarily used by Phrap (www.phrap.org)

my $linenum = shift(@ARGV); #the line number

#This pipeline requires several output files from te_alignment_processing - user must specify the full path to these files
#the mdist file - must specify the full path
my $mapout = shift(@ARGV) or die;

#The umu file - must specify the full path
my $umu = shift(@ARGV) or die;

#The bam file - must specify the full path - the user must provide a single merged bam file
my $bam = shift(@ARGV) or die;

#The file with TE positions in the reference sequence - the user must provide the full path
#this file is tab delimited with three columns
#chromosomenumber start stop
my $TE_file = shift(@ARGV) or die; #for D. melanogaster - TE_positions - provided with scripts

#The file with the positions that are NOT TEs - provided for D. melanogaster
my $nottefile = shift(@ARGV) or die; #NOT_TE_positions

#The formatted file with reference TEs - for calling TEs known in the reference
my $ref_teclust_output = shift(@ARGV) or die; #a formatted list of known TE positions - provided for D. melanogaster - reference_teclust.output

#This indicates which strand a read is on based on the sam flag  - provided sam_flags
my $flag = shift(@ARGV) or die; # a file with sam flags

#The user must specify and output directory
my $outdir = shift(@ARGV) or die;#full path of the output directory

#script commands
##################################################################################
my $teclust = 'teclust';
my $sam = "samtools";
my $view = "view";
my $L = "-L";
my $filter = 'filter_edit';
my $pair = 'get_pair_ids_bwa';
my $pairmapping = 'get_pair_ids_mapping_bwa_2.pl';
my $line_reprocess = 1; #this reprocess fastq files if they have been formatted to be uploaded as sql tables - this is necessary to create fastq files for phrap that won't crash 1 is won't 0 is will
##################################################################################

##################################################################################

#file definitions to leave alone - these are all output files - some of which are temporary files
##################################################################################
my $output = $outdir . "line" . $linenum . "_te_match.sam";#first sam file - reads that overlap TE regions
unlink(qq{$output});
my $output2 = $outdir . "line" . $linenum . "_te_match_ranges";#file of positions for making the second bam file - save line/lane/pair for later matching and store positions of mates of TE reads to use samtools view to grab reads in these positions
unlink(qq{$output2});
my $output3 = $outdir . "line" . $linenum . "_mates.sam";#second sam file - of mates of TE reads
unlink(qq{$output3});
my $output4 = $outdir . "line" . $linenum . "_mates.bam";#output3 sam file to bam format
unlink(qq{$output4});
my $output5 = $outdir . "line" . $linenum . "_unique.sam";#samtools view to grab reads that are in NOT TE ranges
unlink(qq{$output5});
my $output6 = $outdir . "line" . $linenum . "_teclust_input";#teclust input - verify that reads have line/lane/pair of a TE read
unlink(qq{$output6});
my $ummoutfile = $outdir . "line" . $linenum . "_umu_teclust_input";#teclust input - straight from umu file
unlink(qq{$ummoutfile});
my $teclustout = $outdir . "line" . $linenum . ".teclust_output.gz";
unlink(qq{$teclustout});
my $filtered = $outdir . "line" . $linenum . "_teclust_filtered.gz";
unlink(qq{$filtered});
my $filtered_cat = $outdir. "line" . $linenum . "_teclust_filtered_all";
unlink(qq{$filtered_cat});
my $filtered_cat2 = $filtered_cat . ".gz";
unlink(qq{$filtered_cat2});
my $pairout = $outdir . "line" . $linenum . "_intermediate";
unlink(qq{$pairout});
my $pairout2 = $outdir . "line" . $linenum . "_te_pos_estimate";
unlink(qq{$pairout2});
##################################################################################

#read .mdist file to get the 99% percentile of mapping fragment distances
##################################################################################
my $p99 = 0;
open(A, "gunzip -c $mapout |") || die "can't open pipe to $mapout";
 LOOP:while(my $line2 = <A>) {     
     if($line2 =~ m/#/) {         
         next LOOP;
     }    
     if($line2 =~ m/distance/) {        
         next LOOP;
     }
     $line2 =~ s/\s+/\t/go;     
     my @a = split("\t", $line2);
     my @b = split("e-", $a[2]);
     if($b[1] != 01){
         next LOOP;
     }
     if (($b[1] == 01) and ($b[0] <= 9.9)) {
         $p99 = $a[0];
     }
     if (($b[1] == 01) and ($b[0] > 9.9)) {
         last LOOP;
     }
}
close A;
print "99% of fragment lengths is ", $p99, "\n";
##################################################################################
#store the conversion for sam flags - this indicates the strand of the read
my %flag = ();
$flag{65} =  0;
$flag{69} =  2;
$flag{73} =  0;
$flag{77} =  2;
$flag{81} =  1;
$flag{83} =  1;
$flag{89} =  1;
$flag{97} =  0;
$flag{99} =  0;
$flag{113} =  1;
$flag{117} =  2;
$flag{129} =  0;
$flag{133} =  2;
$flag{137} =  0;
$flag{141} =  2;
$flag{145} =  1;
$flag{147} =  1;
$flag{153} =  1;
$flag{161} =  0;
$flag{163} =  0;
$flag{165} =  2;
$flag{177} =  1;      
$flag{181} =  2;
##################################################################################
system(qq{$sam $view $L $TE_file $bam > $output});#first sam file - reads that overlap TE regions
#run samtools view searching for all reads that overlap positions of known TEs 
open(B, "<$output"); #readthrough the output of the first samtools search - store line/lane/pair for reads matching the search
my %hash = ();
my %llp = ();
while(my $line2 = <B>) {
    chomp $line2;
    my @b = split(/\s+/, $line2);
    my @llp = split(/:/, $b[0]);
    if($b[6] eq "=") {
	$hash{$b[2] . "\t" . $b[7]} = 1; #mates of TE reads - positions
	$llp{$llp[0] . ":" . $llp[1] . ":" . $llp[2]} = 1;
    }else {
	$hash{$b[6] . "\t" . $b[7]} = 1; #mates of TE reads - positions
	$llp{$llp[0] . ":" . $llp[1] . ":" . $llp[2]} = 1;
    }
}
close B;
my %chroms = ();#push hashes of positions into arrays 
while((my $key, my $value) = each(%hash)) {
    my @s = split(/\t/, $key);
    push @{$chroms{$s[0]}}, $s[1];
}
open(C, ">>$output2"); #here we are condensing and printing out all of the positions of mates of TEs
for(my $t = 0; $t <=$highch; $t++){ 
    if(exists($chroms{$t})) {
        my @sorted = sort {$a <=> $b} @{$chroms{$t}};        
        push(@sorted, -1);        
        my $currstart = -1;
        my $currstop = -1;
        my $new = 0;
      LOOP:foreach my $element (@sorted) {            
	  if($element != -1) {                
	      if(($new == 0) and ($currstart == -1)) {		  
		  $currstart = $element;
		  $currstop = $element;
		  $new = 1;
		  next LOOP;		  
	      }elsif($new == 1){                       
		  if(($element - $currstop) == 1) {                          
		      $currstop = $element;
		      next LOOP;                          
		  }else{                          
		      print C $t, "\t", $currstart, "\t", $currstop, "\n";
		      $currstart = $element;
		      $currstop = $element;
		      $new = 0;
		      next LOOP;                           
		  }
	      }elsif(($new == 0) and ($currstart != -1)) {
		  if(($element - $currstop) == 1) {
		      $currstop = $element;
		      $new = 1;
		      next LOOP;                           
		  }else{                           
		      print C $t, "\t", $currstart, "\t", $currstop, "\n";
		      $currstart = $element;
		      $currstop = $element;
		      $new = 0;
		      next LOOP;                         
		  } 
	      }
	  }else{                
	      print C $t, "\t", $currstart, "\t", $currstop, "\n";               
	  }
      }
        @sorted = ();
    }
}
close C;
#pull out all reads in the ranges of the mates
system(qq{$sam $view -h $L $output2 $bam > $output3});

#compress sam to bam
system(qq{$sam $view -bS $output3 > $output4});

system(qq{$sam $view $L $nottefile $output4 > $output5});

open(D, "<$output5");
open(OUT, ">>$output6");
my %uniquehash = ();
LOOP:while(my $final = <D>) {
    chomp $final;
    my @f = split(/\t/, $final);
    if($f[2] =~ m/\*/) {
	next LOOP;
    }
    my @inner = split(/:/, $f[0]);
    if(exists($llp{$inner[0] . ":" . $inner[1] . ":" . $inner[2]})) {
	if($flag{$f[1]} != 2) {
	    $uniquehash{$f[0]} = 1;
	    print OUT $f[3], "\t", $f[2], "\t", $flag{$f[1]}, "\n";
	}
    }
}
close OUT;
close D;
open(UMU, "gunzip -c $umu |") || die;
open(UMUOUT, ">>$ummoutfile");
while(my $line2 = <UMU>) {
    chomp $line2;
    my @b = split(/\t/, $line2);
    if(!(exists ($uniquehash{$b[0] . ":" . $b[1] . ":" . $b[2] . ":" . $b[3]}))) {
        print UMUOUT $b[6], "\t", $b[5], "\t", $b[8], "\n";
	delete $uniquehash{$b[0] . ":" . $b[1] . ":" . $b[2] . ":" . $b[3]};
    }
}
close UMU;
close UMUOUT;
##################################################################################

#run teclust
print "Clustering reads...\n";
system(qq{$teclust $TE_file $p99 $mdist $teclustout $output6 $ummoutfile});
#run filter
print "Filtering clusters...\n";
system(qq{$filter $teclustout $filtered $filternum});
#concatenate shared ranges with novel ranges
system(qq{gunzip $filtered});
$filtered =~ s/\.gz//go;
system(qq{cat $filtered $ref_teclust_output > $filtered_cat});
system(qq{gzip $filtered});
system(qq{gzip $filtered_cat});
$filtered = $filtered . ".gz";
$filtered_cat = $filtered_cat . ".gz";

#run get_pair_ids
print "Finding read pairs...\n";

my @pairscommand = ($pair, " " , $filtered_cat, " ", $p99, " ", $pairout, " ", $pairout2);

push(@pairscommand, " ");
push(@pairscommand, $umu);

system(qq{@pairscommand});

system(qq{perl $pairmapping $linenum $pairout2 $p99 $bam $pairout $flag});

unlink(qq{$output});
unlink(qq{$output2});
unlink(qq{$output3});
unlink(qq{$output4});
unlink(qq{$output5});
print "Done! \n";
