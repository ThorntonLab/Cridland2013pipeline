This code is used to detect structural variation from paired-end Illumina data for systems with a high-quality reference genome.

The code is written by Kevin Thornton (krthornt@uci.edu) and is property of the University of California.

To compile:
   run "make"

   You need libsequence (www.molpopgen.org) and the boost libraries (www.boost.org) installed.  libsequence also depends on the GNU Scientific Library for installation (www.gnu.org/software/gsl).

   Please note that we rely heavily in boost_iostreams, which is a compiled library.  Therefore you need a full boost installation, and not simply the header files.

   NOTE: please refer to the above web pages for instructions on installing libsequence, boost, and the GSL.

   The software is intended to be run on a POSIX-style Unix machine.  All testing and running has been done on Linux, and we guarantee nothing about Apple's OS X, especially if you use third-party means of installing libraries such that they end up in directories other than /usr/local/[include|lib]. 

   We offer no support for Windows machines, even those running Cygwin.  Nor is there any guarantee that this code will work on such a system.

   Any reports of compilation issues and/or bugs must include a copy/paste of the error messages, otherwise we cannot help you.  You must also provide your operating system and version, and your compiler system and version.

Note: see documentation for bwa and samtools for those steps.

Input files:
      1. ref.fa or ref.fa.gz = your reference genome, in FASTA, or gzipped FASTA
      2. left.fastq.gz and right.fastq.gz = left and right read of a lane of PE data
      3. left2.fastq.gz and right2.fastq.gz = left and right read of a second lane of PE data (so that the examples of using > 1 lane make sense)

Pipeline outline:


	 The pipeline starts with FASTQ files, not bam files.  So, if you have already aligned your data, then you cannot easily use those bam files without first renaming all the reads within them.  We do not support that use.

	 1. rename_reference ref.fa[.gz] ref.renamed.fa.gz name_table.csv seq_table.csv
	    Output: ref.renamed.fa.gz contains the reference with sequence names replaced with integers.  name_table.csv is a tab-delimited table mapping integers to original names, and seq_table.csv is a tab-delimited table with the new names and the sequences.  This last table can be deleted--we used to upload it to MySQL, but no longer do so.
	 2. Process ref.renamed.fa.gz with bwa build as you normally would.
	 3. fastq_to_table line_id lane_id read_dir fastqfile outfile.csv.gz outfile.fastq.gz
	    Input: line_id = an integer representing the sample (>= 0)
	    	   lane_id = an integer representing the lane (>= 0)
		   read_dir = 0 for left read, 1 for right read.
		   outfile.csv.gz = name of a gzipped file in tab-delimited format containing line_id lane_id pair_id read_dir followed by the original info from the fastq file.  [NOTE: can be deleted]
		   outfile.fastq.gz = a new fastq file with each read name replaced with line_id:lane_id:pair_id:read_dir, where pair_id is auto-filled during execution.  This read renaming is critical for douwnstream programs to work
		   
		   Example: fastq_to_table 0 0 0 left.fastq.gz left.csv.gz left.renamed.fastq.gz
		   	    fastq_to_table 0 0 1 right.fastq.gz right.csv.gz right.renamed.fastq.gz
			    fastq_to_table 0 1 0 left2.fastq.gz left2.csv.gz left2.renamed.fastq.gz
			    fastq_to_table 0 1 1 right2.fastq.gz right2.csv.gz right2.renamed.fastq.gz

			    And now you have the two lanes of sequence from sample 0 all renamed.

			    
	4.  Align all 4 fastq files to ref.renamed.fa.gz as per normal.  Use "bwa sampe" to resolve PE mappings.  See Cridland et al. for mapping parameters that we used.  
	5.  Merge the bam files across lanes.
	6.  Use samtools to create two different sorted bam files.  One is sorted by position (sample0.sorted.bam), and the other by read name (sample0.readsorted.bam).  The first is the "normal" bam file, and the latter is used by downstream programs
	7.  Run bwa_bam_to_mapfiles to collect "unusual/interesting" read pairs for further analysis:
	    samtools view -f 1 sample0.readsorted.bam | bwa_bam_to_mapfiles sample0.structural sample0.um

	    The output will be several files:
	    sample0.structural_left.gz
	    sample0.structural_right.gz
	    sample0.structural.sam.gz
	    sample0.um_u.gz
	    sample0.um_m.gz
	    sample0.um_sam.gz

	    The *.structural_[left|right].gz files contain information about read pairs in either parallel (PAR) or divergent (DIV) orientation, or pairs that map uniquely to different chromosomes (UL = unlinked).  These files are used for calling structural variants (See Cridland and Thornton (2010) Genome Biology and Evolution 2010: 83-101 for a discussion of how these reads may be used).  The format is tab-delimited:

	    line_id lane_id pair_id read_dir mapping_qual chrom_id start stop strand[0/1=fwd/rev] mismatches num_gaps type[=UL|PAR|DIV]

	    The *.structural.sam file contains the original read data in simple "sam" format as compressed ASCI.  This is really for validation/debugging.

	    In the structural file, each read is uniquely-mapping.

	    The um_u/um_m files are use for detection of TE insertions.  The um_u file contains information for a uniquely-mapping read, and the um_m files containg the mappings for the multiply-mapping reads.

	    The um_u/um_um files have the following format:

	    line_id lane_id pair_id read_dir mapping_qual chrom_id start stop strand[0/1=fwd/rev] mismatches num_gaps 

	    The difference between these and the "structural" files is that the um_u files have 1 line per read because the reads map uniquely, while the um_m file has > 1 line per read.

	8.  Get the distribution of mapping distances between unique pairs in the proper orientation:

	    samtools view -f 1 sample0.readsorted.bam | bwa_mapdistance sample0.mdist.gz

	    The output is the empirical cumulative distribution function of mapping distances, and is easily processed in R to find any quantile.

	    For example:

	    zcat /home/krthornt/new_yakuba_data/alignments_bwa_analysis/cnv2/line01_CY20A.mdist.gz | head -n 10
	    distance	number	cprob
	    74	7	5.914142e-08
	    75	5	1.013853e-07
	    76	25	3.126047e-07
	    77	1389	1.204795e-05
	    78	1506	2.477181e-05
	    79	1665	3.883902e-05
	    80	1635	5.265276e-05
	    81	1488	6.522454e-05
	    82	1464	7.759355e-05    

	    distance = distance b/w where reads map (mapping start position) in bp
	    number = number of reads at that distance
	    cprob = cumulative probability

	    To find the 99.9th quantile of the mapping distance in R:

	    > x=read.table("/home/krthornt/new_yakuba_data/alignments_bwa_analysis/cnv2/line01_CY20A.mdist.gz",header=T)
	    > z=which(x$cprob >= 0.999)
	    > x[z[1],]
	        distance number     cprob
		925      998    279 0.9990015

	    And therefore 925bp is the mapping distance corresponding to the 99th quantile of mapping distances


       9.  Cluster the um_u/um_m files into putative TE calls:
	
		umm_te_finder sample0.um_u.gz sample0.um_m.gz min_mapping_qual sample0.tefinder_out.gz refdata.txt

		Where:
			sample0.[um_u|um_m.gz] are the output of step 7.  
			min_mapping_qual is the minimum mapping quality to allow for unique reads [ >= 0 ]
			sample0.tefinder_out.gz = the name of the output file (gzipped text)
			refdata.txt = tab-delimited uncompressed text file with 3 columns (all integers): chrom start stop, representing the positions of TEs in ref.renamed.fa.gz


		The output will be three columns (integers):

		position   chrom   strand

		The positions correspond to the mapping start positions of uniquely-mapping reads whose multiply-mapping partners map to intervals specified in refdata.txt


		Then:

		teclust  refdata.txt insert_size max_dist outfile.gz


		Where:
			refdata.txt = same file that goes into umm_te_finder describing TEs in ref.renamed.fa.gz
			insert_size = An estimate of the upper bound on insert sizes of your library.  For example, the 99.9th quantile obtained in step 8.  This is used
				      as a cutoff for merging reads into clusters
		        max_dist = maximum distance for merging up- and down- stream clusters into the same "event"
			outfile.gz = output file name [gzipped]

		The contents of outfile.gz are the following columns:
