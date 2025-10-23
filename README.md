cd # SaugerParentage
Created by Sam Johnson, 10/25/24

sjohn208@uwyo.edu

# Project Summary

Below is the code and notes for the Wind River Sauger Parentage Analysis bioinformatics workflow. This began on 10/25/24 with file transfers and establishment of the Sauger Parentage github repository, followed by unzipping files (10/29/24) and parsing (11/01/24).

Problems arose in the alignment phase, where vastly different datasets were produced by aligning to either the yellow perch or walleye reference genome. Information regarding these alignments are deonted by: 
Trial 1: Alignment to the Yellow Perch genome (Perca flavescens) produced very few loci after initial filtering. 
Trial 2: Alignment to the Walleye refernce genome (Sander vitreus)

# Workflow Outline

* Demultiplex
   * Unzip Raw .fastq's

   * Count Raw Reads

   * Fastqc Analysis

   * Splitting Raw Files

   * Parsing

   * Debugging perl on Medbow

* Split .fastq

   * Split to 1184 fastq files

   * Count Raw Reads Per Individual

* Alignment and Choice of Reference Genome (Yellow Perch and Walleye, Trials 1 and 2, respectively)

   * Remove 60 line endings

   * Change Reference names

   * Indexing

   * Run bwa

   * Count Aligned Reads per Individual

* Variant Calling (Yellow Perch and Walleye, Trials 1 and 2, respectively)

   * .sam to .bam

   * Make bam list

   * Call Variants

* Filtering

   * Make ID File for reheadering

   * Reheader

   * First Filter investigation (Both References)

* Data Quality Investigations
   
   * Percent Missing Data Per Locus Investigation (Both References)

   * Generate RadMapReport (Both References)

* BWA MEM Investigation

   * Running BWA MEM

   * Variant Calling BWA MEM, Reheadering

   * SFS Analysis for all combinations of reference genomes and alignment algorithms

* Sequoia

   * Generating input for Sequoia
   
   * Sequoia using F0's & Test F1's (BWA ALN to Walleye Reference)

   * Desperate Times...

   * Picking Up Speed on Sequoia...

   * Expanding Beyond the Test Group....

* Principal Component Analysis

* Simulation Work
   
   * Adding Missing Data

* Bedtools Intersect (aln/samse and mem overlap)

* HipHop

   * randomForest

* Contaminant Filtering

   * fastp filtering

* Entropy

# Demultiplex

### Unzip raw .fastq's

Requires an ssh key for your hpc server to connect it with the repository. Then you'll git clone that repository to your server's account. Unzip using this script from JPJ.

```{bash}
sbatch slurm_zip.sh
```

### Count Raw Reads

See code below. Establish an allocated interactive job, then you're checking for n raw reads using the @ symbols that separate reads in the .fastq file structure.

```{bash}
salloc --account=ysctrout --time=3:00:00
grep -c "^@" 1SaugEvens.fastq
   # 1,111,224,597
grep -c "^@" 1SaugOdds.fastq
   # 1,156,307,046
```

### Fastqc Analysis

Fastqc runs various analyses on your raw data to get you a cursory look before you expand to other analyses.

Start by running fastqc on your raw data files (the large ones from each libarary) using this slurm script . In this case, those are 1SaugEvens.fastq and 1SaugOdds.fastq.

```{bash}
sbatch slurm_fastqc.sh
```

Outputs will be returned as zipped directories. Unzip them using: 

```{bash}
salloc --account=ysctrout --time=3:00:00
unzip 1SaugEvens_fastqc.zip
unzip 1SaugOdds_fastqc.zip
```

Notice that these are not .gz files, so we do not use gunzip.
We'll then need to go into these output directories and scp the .htmls to our local machine to open in Chrome.
```{bash}
# (Run locally, from /Users/samjohnson/Documents/Sauger_102824/Fastqc_Output)

scp -r sjohn208@medicinebow.arcc.uwyo.edu:/project/ysctrout/hatchsauger/1Saug/rawreads/out_fastqc/1SaugEvens_fastqc .

scp -r sjohn208@medicinebow.arcc.uwyo.edu:/project/ysctrout/hatchsauger/1Saug/rawreads/out_fastqc/1SaugOdds_fastqc .
```


### Splitting Raw Files
Since it would take incrdibly long to parse the 1SaugEvens/Odds.fastqs, we need to make them smaller to parallelize the parsing process. In this case, we're splitting the big files into smaller ones that are each 55 million lines long.

```{bash}
split -d -l 55000000 --verbose --additional-suffix \.fastq 1SaugEvens.fastq even_fq_
split -d -l 55000000 --verbose --additional-suffix \.fastq 1SaugOdds.fastq odd_fq_
```

### Parsing

Created new scripts to parallelize the parsing process across all split files in evens and odds that start with "even_fq" or "odd_fq". Those scripts are run_parse_evens.pl and run_parse_odds.pl.


```{bash}
perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/run_parse_evens.pl even_fq_*
perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/run_parse_odds.pl odd_fq_*
```
Parsed files contain reads that start with the individual sample id (SAR_YY_XXXX). However, the files are not sorted by sample id. That's the split .fastq step (see below). 

The parsing process also generates reports, and these reports describe, for the parsed files, how many good mids there are among all of the reads.

Count the number of good mids within each parsed report. To do this, you'll need to load R as a module.

```{bash}
## how many good mids?

grep "Good mids" parsereport_* | cut -f 4 -d " " > good_mids_count.txt
   #for all of these parsed report files, take the fourth column of the line that has "Good mids"

module load arcc/1.0 gcc/14.2.0 r/4.4.0
R

dat <- read.delim("good_mids_count.txt", header = FALSE)
sum(dat[,1])
	## 166 files, 1,394,271,363 good mids

#exit R with 
quit()
```
This stayed the same after fixing the demux error (from 11/27) on 12/18/24, and rerunning it on 12/19/24. No surprises there. Number of good mids should stay the same.

Next step is to concatenate all parsed files into a single file called all_parsed.fastq using slurm_cat.sh
Now, this file contains ALL of the reads that correspond to sample id's, which theoretically, should be the same number of all of the good mids from the above step. In other words, n reads in all_parsed.fastq should = 1,394,271,363. 

```{bash}
sbatch slurm_cat.sh
```

And when you run this...
```{bash}
salloc --account=ysctrout --time=3:00:00
grep -c "^@" all_parsed.fastq
# 1,394,271,363
```
You do in fact get the number of good mids!

### Debugging perl on Medbow

After the transition, perl parsing script no longer worked. Required this command to activate perl module.

```{bash}
 /usr/bin/perl -MCPAN -e'install Text::Levenshtein::XS'
```
No longer required to run the perl module at the beginning of this script. Reason unknown. 

# Split .fastq
### Splitting to 1184 .fastq files
Make ID's file and check nrow using...

```{bash}
grep -h "_" *Demux.csv | cut -f 3 -d "," > sauger_ids.txt
wc -l sauger_ids.txt
# 1184 (12 plates * 96) + 32 other samples
```

Splitting the concatenated all_parsed.fastq into individual .fastq files that correspond to each individual Sample ids using slurm_sauger_split.sh 

First, the script makes empty files in your rawreads directory (SAR_YY_XXXX.fastq, n = 1184), then assigns each read in all_parsed.fastq to those sample-specific .fastq's.

```{bash}
sbatch slurm_sauger_split.sh
```
Keep in mind, this slurm script is running an embedded perl script, splitFastq_universal_regex.pl. In the last few lines of that slurm script. You're setting the working directory to the raw reads, then providing the file path to the perl script, then the id's text file, and the all_parsed.fastq. 

### Count Raw Reads Per Individual
We're concerned about this lack of sites, and whether or not its a structural problem with the data. Katie's asked to see a histogram of nreads per individual.

To get the number of raw reads per individual (from the raw .fastq files), and store them in a .txt file, we'll use this code (run from the rawfastqs directory).
```{bash}
grep -c "^@" *.fastq > nreads_per_ind.txt

# That .txt file is located in this directory
/Users/samjohnson/Documents/Sauger_102824/GeneticData/AssembledReads/n_rawreads_per_indiv
```

We'll then load this .txt into R and generate a histogram. (see nrawreads.R)

# Alignment and Choice of Reference Genome

### Yellow Perch Reference (Trial 1)

#### Remove 60 line endings
First step is to remove the line endings after 60 characters, using remove60_fasta.pl
```{bash}
## remove 60 line endings
perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/remove60_fasta.pl Perca_flavescens.fasta
```

#### Change Reference Names
Second step is to rename the reference names into something a bit more meaningful/reasonable.
```{bash}
## change reference names to something shorter

perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/rename_scaff.pl no60_Perca_flavescens.fasta

mv renamed_no60_Perca_flavescens.fasta.txt yellowperch_genome.fna
```
#### Indexing
Third step is to make an index for the reference genome bwa (Burrow Wheeler Align).
```{bash}
sbatch slurm_sauger_bwa_index.sh
```

#### Run bwa
Now you're ready to run Burrow Wheeler Align on /project/ysctrout/hatchsauger/1Saug/rawfastqs/*.fastq

Generates for you a sam_sai directory that houses the output from this step, aligned reads with assoc. alignment scores. For each individual's split .fastq, you get a .sam and .sai file (standard alignment map, index, respectively). Can less into .sam but not .sai
```{bash}
sbatch slurm_sauger_runbwa.sh
```

#### Count Aligned/Assembled Reads (Yellow Perch)
Run this slurm script from your sam_sai_pflav directory. It will output a .txt file that you can load into R and use to plot a histogram of the number of aligned reads per individual and compare it to that of the walleye aligned reads. (See /project/ysctrout/hatchsauger/SaugerParentage/r_scripts/assembledreads/nalignedreads.R)

```{bash}
sbatch /project/ysctrout/hatchsauger/SaugerParentage/slurms/count_assembled.sh
```

### Walleye Reference (Trial 2)
Download Walleye Reference Genome to new directory in ysc/reference_genomes/ called Sander_vitreus and unzip.
```{bash}
rsync --copy-links --recursive --times --verbose rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/031/162/955/GCA_031162955.1_sanVit1/*.fna.gz .
gunzip GCA_031162955.1_sanVit1/*.fna.gz
```
#### Remove 60 line endings
First step is to remove the line endings after 60 characters, using remove60_fasta.pl
```{bash}
## remove 60 line endings
perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/remove60_fasta.pl GCA_031162955.1_sanVit1_genomic.fna 
```
#### Change Reference Names
Second step is to rename the reference names into something a bit more meaningful/reasonable.
```{bash}
## change reference names to something shorter

perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/rename_scaff.pl no60_GCA_031162955.1_sanVit1_genomic.fna

mv renamed_no60_GCA_031162955.1_sanVit1_genomic.fna.txt walleye_genome.fna
```

#### Indexing
Third step is to make an index for the reference genome bwa (Burrow Wheeler Align).
```{bash}
sbatch slurm_sauger_bwa_index.sh
```

#### Run bwa
Run Burrow Wheeler Align on /project/ysctrout/hatchsauger/1Saug/rawfastqs/*.fastq

Generates for you a sam_sai directory that houses the output from this step, aligned reads with assoc. alignment scores. For each individual's split .fastq, you get a .sam and .sai file (standard alignment map, index, respectively). Can less into .sam but not .sai
```{bash}
sbatch slurm_sauger_runbwa.sh
```

#### Count Aligned/Assembled Reads Per Individual (Walleye)
Run this slurm script from your sam_sai_svit directory. It will output a .txt file that you can load into R and use to plot a histogram of the number of aligned reads per individual and compare it to that of the yellow perch aligned reads (See nalignedreads.R).

```{bash}
sbatch /project/ysctrout/hatchsauger/SaugerParentage/slurms/count_assembled.sh
```

# Variant Calling and Filtering

## Variant Calling (Yellow Perch)

### .sam to .bam (Yellow Perch)
```{bash}
sbatch slurm_sauger_sam2bam.sh
```
Note: Delete .sam's and .sai's after converting to .bam's and .bai's to save space. Navigate to directory, rm *.sam, rm *.sai

### Make bam list (Yellow Perch)
Variants slurm script requires a text file with all of the sorted.bam file names. This generates that file.
```{bash}
ls *.sorted.bam > bam_list.txt
```

### Call variants (Yellow Perch)
This script relies on samtools, bcftools, and paths to both the reference genome and to your sam_sai_pflav directory (which should also house your bam_list.txt from above).
```{bash}
sbatch slurm_sauger_variants.sh
```
Outputs a .vcf, aligned reads and variants for each scaffold.

To see how many variants were called, run: 
```{bash}
grep -c "^scaff" variants_rawfiltered_012325.vcf
79,272 variants called
```

## Filtering (Yellow Perch)
All taking place inside sam_sai_pflav.

```{bash}
salloc --account=ysctrout --time=3:00:00
```

### Making ID file for reheadering:
At this stage, all of the reads are assigned to names aln_Sample_ID.sorted.bam.txt. This command takes those names, cuts off the "aln_" and "sorted.bam", and stores new names in sauger_ids_col.txt that are just the Sample_ID (SAR_YY_XXXX).
```{bash}
sed -s "s/aln_//" bam_list.txt | sed -s "s/.sorted.bam//" > sauger_ids_col.txt
```

### Reheader 
This "reheader"ing step now takes those polished names and assigns the reads in variants_rawfiltered_012325.vcf to those names.
```{bash}
module load arcc/1.0 gcc/14.2.0 bcftools/1.20

bcftools reheader -s sauger_ids_col.txt variants_rawfiltered_012325.vcf -o rehead_variants_rawfiltered_012325.vcf
```

### First filter investigation 
Worked with Maria for a while to get the run_first_filter.pl script to run on rehead_variants_rawfiltered_012325.vcf.

Because you're running this from the sam_sai_pflav directory, you need to provide the whole path to the perl script.

```{bash}
perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/run_first_filter_MPR.pl rehead_variants_rawfiltered_012325.vcf
```
Now we have (in sam_sai_pflav/first_filter_out) a series of standard output files for maf(1,2,3,4,5) and miss(6,7,8,9). 

```{bash}
grep "Sites" first_filter_out/*
first_filter_out/stdout_maf1_miss6:After filtering, kept 218 out of a possible 79272 Sites
first_filter_out/stdout_maf1_miss7:After filtering, kept 195 out of a possible 79272 Sites
first_filter_out/stdout_maf1_miss8:After filtering, kept 169 out of a possible 79272 Sites
first_filter_out/stdout_maf1_miss9:After filtering, kept 144 out of a possible 79272 Sites
first_filter_out/stdout_maf2_miss6:After filtering, kept 185 out of a possible 79272 Sites
first_filter_out/stdout_maf2_miss7:After filtering, kept 167 out of a possible 79272 Sites
first_filter_out/stdout_maf2_miss8:After filtering, kept 144 out of a possible 79272 Sites
first_filter_out/stdout_maf2_miss9:After filtering, kept 125 out of a possible 79272 Sites
first_filter_out/stdout_maf3_miss6:After filtering, kept 161 out of a possible 79272 Sites
first_filter_out/stdout_maf3_miss7:After filtering, kept 150 out of a possible 79272 Sites
first_filter_out/stdout_maf3_miss8:After filtering, kept 128 out of a possible 79272 Sites
first_filter_out/stdout_maf3_miss9:After filtering, kept 112 out of a possible 79272 Sites
first_filter_out/stdout_maf4_miss6:After filtering, kept 143 out of a possible 79272 Sites
first_filter_out/stdout_maf4_miss7:After filtering, kept 133 out of a possible 79272 Sites
first_filter_out/stdout_maf4_miss8:After filtering, kept 113 out of a possible 79272 Sites
first_filter_out/stdout_maf4_miss9:After filtering, kept 99 out of a possible 79272 Sites
first_filter_out/stdout_maf5_miss6:After filtering, kept 135 out of a possible 79272 Sites
first_filter_out/stdout_maf5_miss7:After filtering, kept 125 out of a possible 79272 Sites
first_filter_out/stdout_maf5_miss8:After filtering, kept 108 out of a possible 79272 Sites
first_filter_out/stdout_maf5_miss9:After filtering, kept 94 out of a possible 79272 Sites
```
Went back into run_first_filter_MPR.pl and killed the flag in the vcftools line (57) that tells it to remove sites that are within 100 bp of one another. Slight improvement, but minimal.
```{bash}
grep "Sites" first_filter_out/*
first_filter_out/stdout_maf1_miss6:After filtering, kept 300 out of a possible 79272 Sites
first_filter_out/stdout_maf1_miss7:After filtering, kept 258 out of a possible 79272 Sites
first_filter_out/stdout_maf1_miss8:After filtering, kept 199 out of a possible 79272 Sites
first_filter_out/stdout_maf1_miss9:After filtering, kept 166 out of a possible 79272 Sites
first_filter_out/stdout_maf2_miss6:After filtering, kept 254 out of a possible 79272 Sites
first_filter_out/stdout_maf2_miss7:After filtering, kept 220 out of a possible 79272 Sites
first_filter_out/stdout_maf2_miss8:After filtering, kept 167 out of a possible 79272 Sites
first_filter_out/stdout_maf2_miss9:After filtering, kept 142 out of a possible 79272 Sites
first_filter_out/stdout_maf3_miss6:After filtering, kept 228 out of a possible 79272 Sites
first_filter_out/stdout_maf3_miss7:After filtering, kept 202 out of a possible 79272 Sites
first_filter_out/stdout_maf3_miss8:After filtering, kept 151 out of a possible 79272 Sites
first_filter_out/stdout_maf3_miss9:After filtering, kept 129 out of a possible 79272 Sites
first_filter_out/stdout_maf4_miss6:After filtering, kept 208 out of a possible 79272 Sites
first_filter_out/stdout_maf4_miss7:After filtering, kept 183 out of a possible 79272 Sites
first_filter_out/stdout_maf4_miss8:After filtering, kept 135 out of a possible 79272 Sites
first_filter_out/stdout_maf4_miss9:After filtering, kept 115 out of a possible 79272 Sites
first_filter_out/stdout_maf5_miss6:After filtering, kept 198 out of a possible 79272 Sites
first_filter_out/stdout_maf5_miss7:After filtering, kept 174 out of a possible 79272 Sites
first_filter_out/stdout_maf5_miss8:After filtering, kept 129 out of a possible 79272 Sites
first_filter_out/stdout_maf5_miss9:After filtering, kept 109 out of a possible 79272 Sites
```

### First filter investigation (no maf, missing data up to 60%)
Solo work on 02/05/25 to develop run_first_filter_noMAF.pl script to run on rehead_variants_rawfiltered_012325.vcf.

Should have no MAF filter, and run for the following missing data values. 

```{bash}
my @misses = ('9', '8', '7', '6', '5', '4');
# 10 to 60% missing data
```

Because you're running this from the sam_sai_pflav directory, you need to provide the whole path to the perl script.

```{bash}
perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/run_first_filter_noMAF.pl rehead_variants_rawfiltered_012325.vcf
```

Now we have (in sam_sai_pflav/first_filter_out_noMAF) a series of standard output files for miss(4,5,6,7,8,9). 

```{bash}
grep "Sites" first_filter_out_noMAF/*
first_filter_out_noMAF/stdout_miss4:After filtering, kept 5461 out of a possible 79272 Sites
first_filter_out_noMAF/stdout_miss5:After filtering, kept 4772 out of a possible 79272 Sites
first_filter_out_noMAF/stdout_miss6:After filtering, kept 4388 out of a possible 79272 Sites
first_filter_out_noMAF/stdout_miss7:After filtering, kept 4058 out of a possible 79272 Sites
first_filter_out_noMAF/stdout_miss8:After filtering, kept 3685 out of a possible 79272 Sites
first_filter_out_noMAF/stdout_miss9:After filtering, kept 3318 out of a possible 79272 Sites
```


## Variant Calling (Walleye)
All taking place inside sam_sai_svit.

### .sam to .bam (Walleye)

```{bash}
sbatch slurm_sauger_sam2bam.sh
```
Note: Delete .sam's and .sai's after converting to .bam's and .bai's to save space. Navigate to directory, rm *.sam, rm *.sai

### Make bam list (Walleye)
The variant calling slurm script requires a text file with all of the sorted.bam file names. This generates that file.
```{bash}
ls *.sorted.bam > bam_list.txt
```

### Call variants (Walleye)
This script relies on samtools, bcftools, and paths to both the reference genome and to your sam_sai directory (which should also house your bam_list.txt from above).
```{bash}
sbatch slurm_sauger_variants.sh
```
Outputs a .vcf, aligned reads and variants for each scaffold.

To see how many variants were called, run: 
```{bash}
grep -c "^scaff" variants_rawfiltered_svit_020625.vcf
# 366,797 variants called
```

## Filtering (Walleye)
All taking place inside sam_sai_svit.

```{bash}
salloc --account=ysctrout --time=3:00:00
```

### Making ID file for reheadering:
At this stage, all of the reads are assigned to names aln_Sample_ID.sorted.bam.txt. This command takes those names, cuts off the "aln_" and "sorted.bam", and stores new names in sauger_ids_col.txt that are just the Sample_ID (SAR_YY_XXXX).
```{bash}
sed -s "s/aln_//" bam_list.txt | sed -s "s/.sorted.bam//" > sauger_ids_col.txt
```

### Reheader 
This "reheader"ing step now takes those polished names and assigns the reads in variants_rawfiltered_svit_020625.vcf to those names.
```{bash}
module load arcc/1.0 gcc/14.2.0 bcftools/1.20

bcftools reheader -s sauger_ids_col.txt variants_rawfiltered_svit_020625.vcf -o rehead_variants_rawfiltered_svit_020625.vcf
```

### First filter investigation 
Solo work on 02/07/25 to use run_first_filter_MPR.pl script on rehead_variants_rawfiltered_svit_020625.vcf.

Because you're running this from the sam_sai_svit directory, you need to provide the whole path to the perl script.

```{bash}
perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/run_first_filter_MPR.pl rehead_variants_rawfiltered_svit_020625.vcf
```
Now we have (in sam_sai_svit/first_filter_out) a series of standard output files for maf(1,2,3,4,5) and miss(4,5,6,7,8,9).

```{bash}
grep "Sites" first_filter_out/*
first_filter_out/stdout_maf1_miss4:After filtering, kept 10491 out of a possible 366797 Sites
first_filter_out/stdout_maf1_miss5:After filtering, kept 9151 out of a possible 366797 Sites
first_filter_out/stdout_maf1_miss6:After filtering, kept 8396 out of a possible 366797 Sites
first_filter_out/stdout_maf1_miss7:After filtering, kept 7775 out of a possible 366797 Sites
first_filter_out/stdout_maf1_miss8:After filtering, kept 7225 out of a possible 366797 Sites
first_filter_out/stdout_maf1_miss9:After filtering, kept 6488 out of a possible 366797 Sites
first_filter_out/stdout_maf2_miss4:After filtering, kept 8772 out of a possible 366797 Sites
first_filter_out/stdout_maf2_miss5:After filtering, kept 7753 out of a possible 366797 Sites
first_filter_out/stdout_maf2_miss6:After filtering, kept 7222 out of a possible 366797 Sites
first_filter_out/stdout_maf2_miss7:After filtering, kept 6750 out of a possible 366797 Sites
first_filter_out/stdout_maf2_miss8:After filtering, kept 6331 out of a possible 366797 Sites
first_filter_out/stdout_maf2_miss9:After filtering, kept 5736 out of a possible 366797 Sites
first_filter_out/stdout_maf3_miss4:After filtering, kept 7881 out of a possible 366797 Sites
first_filter_out/stdout_maf3_miss5:After filtering, kept 6984 out of a possible 366797 Sites
first_filter_out/stdout_maf3_miss6:After filtering, kept 6519 out of a possible 366797 Sites
first_filter_out/stdout_maf3_miss7:After filtering, kept 6113 out of a possible 366797 Sites
first_filter_out/stdout_maf3_miss8:After filtering, kept 5757 out of a possible 366797 Sites
first_filter_out/stdout_maf3_miss9:After filtering, kept 5228 out of a possible 366797 Sites
first_filter_out/stdout_maf4_miss4:After filtering, kept 7204 out of a possible 366797 Sites
first_filter_out/stdout_maf4_miss5:After filtering, kept 6396 out of a possible 366797 Sites
first_filter_out/stdout_maf4_miss6:After filtering, kept 5988 out of a possible 366797 Sites
first_filter_out/stdout_maf4_miss7:After filtering, kept 5625 out of a possible 366797 Sites
first_filter_out/stdout_maf4_miss8:After filtering, kept 5310 out of a possible 366797 Sites
first_filter_out/stdout_maf4_miss9:After filtering, kept 4849 out of a possible 366797 Sites
first_filter_out/stdout_maf5_miss4:After filtering, kept 6702 out of a possible 366797 Sites
first_filter_out/stdout_maf5_miss5:After filtering, kept 5956 out of a possible 366797 Sites
first_filter_out/stdout_maf5_miss6:After filtering, kept 5597 out of a possible 366797 Sites
first_filter_out/stdout_maf5_miss7:After filtering, kept 5268 out of a possible 366797 Sites
first_filter_out/stdout_maf5_miss8:After filtering, kept 4978 out of a possible 366797 Sites
first_filter_out/stdout_maf5_miss9:After filtering, kept 4554 out of a possible 366797 Sites
```

### First filter investigation (no maf, missing data up to 60%)
Solo work on 02/07/25 to use run_first_filter_noMAF.pl script on rehead_variants_rawfiltered_svit_020625.vcf.

Should have no MAF filter, and run for the following missing data values. 

```{bash}
my @misses = ('9', '8', '7', '6', '5', '4');
# 10 to 60% missing data
```

Because you're running this from the sam_sai_svit directory, you need to provide the whole path to the perl script.

```{bash}
perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/run_first_filter_noMAF.pl rehead_variants_rawfiltered_svit_020625.vcf
```

Now we have (in sam_sai_svit/first_filter_out_noMAF) a series of standard output files for miss(4,5,6,7,8,9). 

```{bash}
grep "Sites" first_filter_out_noMAF/*
first_filter_out_noMAF/stdout_miss4:After filtering, kept 87350 out of a possible 366797 Sites
first_filter_out_noMAF/stdout_miss5:After filtering, kept 78726 out of a possible 366797 Sites
first_filter_out_noMAF/stdout_miss6:After filtering, kept 74268 out of a possible 366797 Sites
first_filter_out_noMAF/stdout_miss7:After filtering, kept 70270 out of a possible 366797 Sites
first_filter_out_noMAF/stdout_miss8:After filtering, kept 66084 out of a possible 366797 Sites
first_filter_out_noMAF/stdout_miss9:After filtering, kept 59169 out of a possible 366797 Sites
```

# Data Quality Investigations

## Percent Missing Data Per Locus Investigation

We are concerned about the low number of reads that were retained for the yellow perch allignment. Is it possible that there's a lot of missing data for each site? This script investigates that.

(From slurm scripts directory)

This script generates a summary table for each site on each scaffold that describes, among other things, the percent missing data for that site. You need to give it the .vcf of interest and the path to that .vcf.

```{bash}
sbatch slurm_sauger_permd.sh
```

This script outputs a .missing.lmiss tab delimited file that we can read into R and use to plot histograms of the percent missing data per site.

## RadMap Report Generation

Generating a RadMapReport is another way that we can investigate the quality of our data. This report gives you information (for each individual) on the number of raw reads, aligned reads, number of loci, mean coverage, number of loci with 10x coverage or greater, and how many loci have that good coverage. This information can be plotted in a variety of ways.

I modifited the script by Joana Meier at this link to create slurm_radmapreport.sh: https://speciationgenomics.github.io/allelicBalance/

 Run the script on your directory of .sorted.bam files (See line 19). (e.g., /project/ysctrout/hatchsauger/sam_sai_svit/sorted.bams)

```{bash}
sbatch slurm_radmapreport_fix.sh
```

We then must rename the output .txt files to describe species using mv. I then scp'd these files to my local machine to read into R.

See instructions in this R Script to transform output from slurm_radmapreport_fix.sh to usable radmapreports, or see reports themselves. 

(From this Repository)
```{bash}
cd /project/ysctrout/hatchsauger/SaugerParentage/r_scripts
open RadMapReport.R
```
OR from local
```{bash}
cd /Documents/Sauger_102824/GeneticData/RadMapReport/RadMap_svit
cd /Documents/Sauger_102824/GeneticData/RadMapReport/RadMap_pflav
```
to see complete reports as they are intended to be used in plotting.


# BWA MEM Investigation

We found out that Will's data were aligned using BWA MEM rather that BWA ALN and SAMSE. Josh and I have generated runbwa_mem.pl and have run it on the walleye and yellow perch reference genomes.

The outputs are stored in /project/ysctrout/hatchsauger/sam_sai_pflav_mem and /project/ysctrout/hatchsauger/sam_sai_svit_mem. These sorted bams were then used to count the number and percentage of aligned reads for each using count_assembled_pflav_mem.sh and count_assembled_svit_mem.sh.

There is a difference in format of mem-aligned reads than aln/samse-aligned reads, so the RadMapReport script does not work for mem-aligned reads. Will need to circle back to this.

Upon looking at the number of assembled reads for each species, they were the same for both bwa_mem outputs. I will go try them again for both species and double check the reference genome indices.

I have rerun bwa mem for each species from empty directories, and I will revisit the results tomorrow. I anticipate that this flaw was due to a failure to git push/pull before running the perl script for the other reference genome.

The r script to plot the histograms of percent aligned reads (nalignedreads.R) was updated and pushed. 

That concludes work 3/12/25. BWA MEM's are currently running, script is ready for when they're done.

BWA MEM's are finishing, and I am starting count_assembled_svit_mem.sh and count_assembled_pflav_mem.sh.

Count_assembled analyses have finished on the bwa mem aligned reads vs the bwa aln samse aligned reads (for both walleye and yellow perch) and the results are shocking. See local directory /GeneticData/AssembledReads/n_aln_per_indiv for plots and /SaugerParentage/r_scripts/assembledreads for the R script.

## Variant Calling (BWA MEM) and SFS Analysis for BWA MEM (Each step was done for svit_mem and pflav_mem)
Because BWA MEM generates .bam, .sorted.bam, and .sorted.bam.bai, I have bypassed using slurm_sauger_sam2bam.sh and have moved to slurm_sauger_variants.sh for both sam_sai_svit_mem and sam_sai_pflav_mem. Variant calling for each of those sets was initiated at 11:30pm on 03/13/25.

03/14/25 - All four directories aln/mem pflav/svit now have rawfiltered .vcfs in them. 

I want to take each of those .vcfs and filter them for MAF 0.01 and Miss 0.9. Then, turn those into genotype  matrices, and plot an SFS for each that I can send to them in a 2x2 panel.

First step is to create reheadered .vcfs for both mem directories:

### Making ID file for reheadering (both directories):
At this stage, all of the reads are assigned to names mem_Sample_ID.sorted.bam in bam_list.txt. This command takes those names, cuts off the "mem" and "sorted.bam", and stores new names in sauger_ids_col.txt that are just the Sample_ID (SAR_YY_XXXX).

```{bash}
sed -s "s/mem_//" bam_list.txt | sed -s "s/.sorted.bam//" > sauger_ids_col.txt
```

### Reheader 
This "reheader"ing step now takes those polished names and assigns the reads in variants_rawfiltered_svit_mem_031325.vcf to those names.

```{bash}
module load arcc/1.0 gcc/14.2.0 bcftools/1.20

bcftools reheader -s sauger_ids_col.txt variants_rawfiltered_svit_mem_031325.vcf -o rehead_variants_rawfiltered_svit_mem_031325.vcf
bcftools reheader -s sauger_ids_col.txt variants_rawfiltered_pflav_mem_031325.vcf -o rehead_variants_rawfiltered_pflav_mem_031325.vcf
```

Josh showed me this amazing line to filter vcfs without dealing with the script.

```{bash}
/project/ysctrout/software/vcftools/bin/vcftools --vcf rehead_variants_rawfiltered_svit_mem_031325.vcf --out variants_maf1_miss9 --remove-filtered-all --maf 0.01 --max-missing 0.9 --recode
/project/ysctrout/software/vcftools/bin/vcftools --vcf rehead_variants_rawfiltered_pflav_mem_031325.vcf --out variants_maf1_miss9 --remove-filtered-all --maf 0.01 --max-missing 0.9 --recode

/project/ysctrout/software/vcftools/bin/vcftools --vcf rehead_variants_rawfiltered_svit_mem_031325.vcf --out variants_maf30_miss9 --remove-filtered-all --maf 0.3 --max-missing 0.9 --recode
```

Filtering output (svit_mem): MAF 0.01 Miss 0.9
```{bash}
After filtering, kept 1184 out of 1184 Individuals
After filtering, kept 7390 out of a possible 404644 Sites
```

Filtering output (svit_mem): MAF 0.30 Miss 0.9
```{bash}
After filtering, kept 1184 out of 1184 Individuals
After filtering, kept 1540 out of a possible 404644 Sites
```

Filtering output (pflav_mem): MAF 0.01 Miss 0.9
```{bash}
After filtering, kept 1184 out of 1184 Individuals
After filtering, kept 204 out of a possible 136652 Sites
```

Just to be sure, I'm also going to do this same filtering with this line in the aln directories, sam_sai_svit and sam_sai_pflav

```{bash}
/project/ysctrout/software/vcftools/bin/vcftools --vcf rehead_variants_rawfiltered_svit_020625.vcf --out variants_maf1_miss9 --remove-filtered-all --maf 0.01 --max-missing 0.9 --recode
/project/ysctrout/software/vcftools/bin/vcftools --vcf rehead_variants_rawfiltered_pflav_012325.vcf --out variants_maf1_miss9 --remove-filtered-all --maf 0.01 --max-missing 0.9 --recode
```

Filtering output (svit aln/samse): 
```{bash}
After filtering, kept 1184 out of 1184 Individuals
After filtering, kept 6488 out of a possible 366797 Sites
```

Filtering output (pflav aln/samse): 
```{bash}
After filtering, kept 1184 out of 1184 Individuals
After filtering, kept 166 out of a possible 79272 Sites
```

Now that we have filtered .vcfs to MAF 0.01 and Miss 0.90 for all four combinations, let's convert them to genotype matrices and prepare for SFS. I'm also going to convert all of them for Sequoia, since I imagine I'll put them in there at some point. Here's the example for pflav_mem

```{bash}
salloc --account=ysctrout --time=3:00:00
module load arcc/1.0 gcc/14.2.0 vcftools/0.1.17
cd /project/ysctrout/hatchsauger/sam_sai_pflav_mem/
vcftools --vcf variants_maf1_miss9.recode.vcf --012 --out variants_maf1_miss9_mem_pflav
sed "s/-1/-9/g" variants_maf1_miss9_mem_pflav.012 > variants_maf1_miss9_mem_pflav.012_conv
# conv = converted to missing data = "-9"
```

SFS generated in SaugerParentage/r_scripts/sfs_all.R on 03/14/25

## radmap Reports from BWA MEM

Remember that there are also a set of scripts to generate radmapreports from data aligned using bwa mem. We were having trouble with this at first. See here:

```
/project/ysctrout/hatchsauger/SaugerParentage/slurms
ls -ltrh *rad*
```

Again, svit mem outperforms everything else (species and alignment tools), and there do not appear to be differences among our two libraries. Great news.

# Sequoia

## Generating Input for Sequoia
It appears that the alignment to the walleye reference generated a much higher quality dataset. Now we are ready to filter start with a standard filter of MAF ≥ 1 and Missing Data < 10% per site.

We're going to use a command in vcftools (--012) to generate a file that lists the genotypes, for each individual, for each locus in terms of 0,1,2 (copies of the reference/nonreference allele), and -1 for missing. Sequoia takes missing data as -9, so we'll need to find and replace every instance of -1 with -9.

```{bash}
salloc --account=ysctrout --time=3:00:00
module load arcc/1.0 gcc/14.2.0 vcftools/0.1.17
cd /project/ysctrout/hatchsauger/sam_sai_svit
vcftools --vcf variants_maf1_miss9.recode.vcf --012 --out variants_maf1_miss9
sed "s/-1/-9/g" variants_maf1_miss9.012 > variants_maf1_miss9.012_conv
# conv = converted to missing data = "-9"
```

## Sequoia F0/Test F1 Analysis (BWA ALN to Walleye Reference)

Sequoia recommends more than 90% genotype call rate, and a MAF of 30% or higher. I have already filtered the rawfiltered .vcf,  rehead_variants_rawfiltered_svit_020625.vcf, using: 

```{bash}
perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/run_first_filter_MPR.pl rehead_variants_rawfiltered_svit_020625.vcf
```

So today, I will start with: 

```{bash}
/project/ysctrout/hatchsauger/sam_sai_svit/first_filter_out/variants_maf3_miss9.recode.vcf
```

...which has 5228 out of a possible 366797 Sites.

I will now generate the genotype matrix for this filtered .vcf using: 

We're going to use a command in vcftools (--012) to generate a file that lists the genotypes, for each individual, for each locus in terms of 0,1,2 (copies of the reference/nonreference allele), and -1 for missing. Sequoia takes missing data as -9, so we'll need to find and replace every instance of -1 with -9.

```{bash}
salloc --account=ysctrout --time=3:00:00
module load arcc/1.0 gcc/14.2.0 vcftools/0.1.17
cd /project/ysctrout/hatchsauger/sam_sai_svit/first_filter_out/
vcftools --vcf variants_maf3_miss9.recode.vcf --012 --out variants_maf3_miss9
sed "s/-1/-9/g" variants_maf3_miss9.012 > variants_maf3_miss9.012_conv
# conv = converted to missing data = "-9"
```
 
 Here's the plan: First, let's review the sequoia scripts we have written so far and take a minute to understand what's going on. Then, take this genotype matrix and load it into R. Take your sauger_ids_col.txt and load it in, along with a .csv of JUST the sample id's of the 2015 F0 fish and their Test F1 progeny. Then, filter the genotype matrix to include rows whose sample id is in the vector of individuals that's Test F1's and 2015 F0's, then edit the Sequoia_Full.R script and run it using slurm_sequoia_first.sh. 

 This was done in the afternoon on 3/13 and finished in the evening. Tomorrow, 3/14 I will scp the data to my local machine and investigate relationships.

No relationships found on that data, perhaps because we filtered with maf of 0.03 above instead of 0.3. Ugh. Resume after break and try again. The filtered .vcf is already there in sam_sai_svit.

03/24/25 Work: Took the maf30_miss9 .vcf and turned it into the genotype matrix. prepared Sequoia_Full.R for a run on that. There are still 1400 sites. May be too many still. We shall see.

Results: Only one parent offspring duo was produced by this run.

03/26/25 Work: Trying maf30_misss9 on svit_mem. Prepared and ran Sequoia_Full.R on that converted matrix. 1540 sites. This should let us know if the number of sites is the problem, or if it's an alignment issue.

Results: NO RELATIONSHIPS RETURNED. 

Weeks of 06/02/25:
I returned to sequoia after learning that hiphop requires the specification of sexes for the putative parents. This makes the use of hiphop impossible for assigning to F2 -> F1's.

My first step upon returning to sequoia as the primary candidate for parentage was to see how it would perform with more strict filtering, and thus with fewer, QUALITY loci.

This still did not return adequate numbers of PO relationships for the Test F1 -> F0 assignment. See results from:

```{bash}
/project/ysctrout/hatchsauger/SaugerParentage/r_scripts/sequoia/Sequoia_ReducedNLoci.R
```

I was then concerned with the strictness of the OH and ME (Opposite Homozygosity and Mendelian Error) filters that are used to filter potential relationships in GetMaybeRel(). I then used this script:

```{bash}
/project/ysctrout/hatchsauger/SaugerParentage/r_scripts/sequoia/gmr_custom.R
```

to customize the GetMaybeRel() function and attempt to loosen those filters. No significant improvements were made to the ability to generate quality assignments for the Test F1's -> F0's.

I then became skeptical of the input structure (of the genotype matrix generated by vcftools --012). I asked Will what he used to create genotype matrices for sequoia. He replied, “I used the vcf2geno function from the LEA package to convert my .vcf into a matrix and then replaced all of the "9" (missing data) with "-9" (what sequoia said it wanted for missing data)”.

I did this and compared the matrix to that created by vcftools –012. There were no discrepancies. This is good, but it doesn’t put me any closer to figuring out what’s going on. For code, see:

```{bash}
/project/ysctrout/hatchsauger/SaugerParentage/r_scripts/sequoia/check_geno_matrices.R
```

### Desperate times...
I have also tried sequoia with a few "hard" filtering appraoches that I originally used for hiphop (see below). As of July 2025 I have been unsuccessful in returning adequate numbers of relationships, even with strict filters for things like depth and quality, lots of turning knobs and tweaking function arguments, and using my custom functions. Eventually I called in the cavalry and emailed Dr. Huisman. Her advice was fantastic. See Sauger Meeting Notes F24/S25/Summer25. Email came 06/30/25.

The first step is to thin our sites for LD, in case we had not yet done this. I believe this happened somewhere upstream in the raw filtering step for fresh .vcf's, but used vcftools --thin to achieve this. pflav_mem and svit_mem .vcfs had now been thinned. svit_mem has yet again emerged with many more sites, see:

```{bash}
/project/ysctrout/hatchsauger/sam_sai_svit_mem/thin_filtered
```

Next step is to try sequoia on these thinned .vcf's (especially svit_mem), then increase the Tassign (use default functions), and alter the complex argument (default functions as well).


### Picking Up Speed on Sequoia...

Lots of knob turning occurred in the summer and fall of 2025. All of it's documented in:

```
cd /project/ysctrout/hatchsauger/SaugerParentage/r_scripts/sequoia
```

I won't document it all here, but big changes include trying earlier sequoia program versions, adding the contaminant filtering described below, incorporating (unsuccessfully) genotype likelihood and correlation-based LD filtering using vcftools. Since I was concerned about the quality of the genotype calls, I also upped min mean read depth to 4 and then 8, and messed with md and thinning to see how many snps I would retain, then created the sequoia_params tables to see how many assignments and accurate assignments we could make for the test group with various combinations of error and thinning. Also created the sequoia_params_plots to visualize.

I also tried leaving out loci that were out of HWE using sequoia's SNPStats() function, but that didn't improve things either. This was implemented in

```
/project/ysctrout/hatchsauger/SaugerParentage/r_scripts/sequoia/Sequoia_ContamFilt_mindep8_md5.R
```

and the output is stored in sequoia_params_snpchip.xlsx (Genetic Data folder)

A breakthrough finally happened when I shared those results with Dr. Jisca Huisman, who explained that the error matrix that informs the pedigree reconstruction has been built on SNP chip data, which has drastically different error rates and structures than RADseq/GBS data. She pointed me in the direction of this sequoia info page,

https://jiscah.github.io/articles/howto_RADseq.html

, as well as this paper by Bresadola et al, (2020, Mol. Ecol. Resources). I was then compelled to switch error types and generate the sequoia_params_radseq tables. Using this radseq error bumped up my compsite scores (Assignment x Accuracy) for the test group, to just over 83%. This was for a few different degrees of thinning, and with a few combinations of E0 and E1.

### Expanding beyond the Test Group...

#### Final Checks (SFS)

Once that breakthrough happened, I was ready to expand to the full dataset. I decided to go with this one:

```
/project/ysctrout/hatchsauger/sam_sai_contam_fastp_svit_mem/vcfs/q40/mindep8/rehead_variants_rawfiltered_svit_mem_contam_fastp_bial_noindels_q40_mindep8_maxdep75_maf30_miss95_thin100K.recode.vcf
```

because it had all of the up front filering for the raw reads, high site quality, high min depth, the right maf for sequoia, and low md per site, all criteria of a good dataset according to the documentation. This dataset got a composite of 83.158% in sequoia_params_radseq.xlsx, and had the most snps of any dataset that got that score. Sequoia says you want lots of snps for small, potentially inbred populations, but that in some cases, fewer snps are better because the dataset are error prone and fewer snps means fewer OH mismatches, which means higher assignment rates. Since this population is small, low diversity, and because I was about to add a lot more individuals by expanding, I wanted to keep as many loci as possible to retain as much signal as possible.

We wanted to undergo a few final checks before expanding, just to be safe. The first of which was to check the SFS for the parent dataset before maf filtering, as well as the data as they were going into sequoia (with final maf and miss filters). The SFS were generated using 

```
/project/ysctrout/hatchsauger/SaugerParentage/r_scripts/sfs_sequoia_ready_datasets
```

and were deemed reasonable, so we pushed forward.

#### Final Checks (AB, Allelic Balance)

During the process of changing error matrices to resemble RADseq rates rather than SNPchip error rates, I became increasingly concerned with allelic dropout, where, because of PCR issues or mutations at cut sites, true heterozygotes only have one or their alleles amplified and genotyped. In the error matrix that I chose, the E1 (error rate for heterozygotes) was much less than E0 (for homozygotes), which seems backwards given the expectation of high rates of ADO in this kind of dataset. Because of that, I wanted to analyze AB at heterozygous sites and for each individual to see if the allelic balances looked off, even though I am using strict filtering.

See these two scrips, modified from MPR, who modified Jessi Rick's approach from an earlier Journal of Heredity paper.

```
/project/ysctrout/hatchsauger/SaugerParentage/slurms/allelicBalance-Site.sh 
```

and

```
/project/ysctrout/hatchsauger/SaugerParentage/slurms/allelicBalance-Samples.sh     
```

The first script takes each site where there are heterozygotes, and takes a mean AB among all individuals genotyped at that site. The second takes each individual, and gives me a mean AB for that individual across all heterozygous sites. I completed these calculations and plotted histograms for these using this r script.

```
/project/ysctrout/hatchsauger/SaugerParentage/r_scripts/AllelicBalancePlots.R
```

Both showed normal distributions centered right around 0.5, but I want to filter out sites/inds with poor AB's. That's coming soon.

#### Running sequoia() and GetMaybeRel() on the full dataset

First trials are stored in:

```
/project/ysctrout/hatchsauger/SaugerParentage/r_scripts/Sequoia_ContamFilt_mindep8_md5_RADseqErr_ALLINDS.R
```

which was created on 102225. First attempts poor. Improvements incoming.




## Principal Component Analysis

Upon learning that sequoia returned 0 relationships between the Test F1's and the 2015 F0 Parents, we wanted to be sure that we didn't get the incorrect plate of test samples from Montana. We performed this validation by running a PCA. If we got the correct samples, they should stack on top of or very close to the rest of our samples. If not, they may be the wrong samples. 

Results: They stack nicely. Few fish from Willie's Lab may be WAE, as well as a few of the F2's. That warrants furtuer investigation. See the results in the local GeneticData/hiphop or GeneticData/pca directories, and see the script for making the pca in: 

```{bash}
/project/ysctrout/hatchsauger/SaugerParentage/r_scripts/pca
```

## A Note on Nomenclature

03/14/25 - Changed all of the .bam, .sorted.bam, and .sorted.bam.bai in both mem directories to start with the prefix "mem_" rather than "aln_". The runbwa_mem.pl script was also updated to reflect this change, should we have to run that alignment again. Keep that in mind for any discrepancies in documentation you may come across down the line.

# Simulation Work

04/09/25 - Today is the firs time that I pulled down Josh's scripts from github and ran them in full, on the cluster. This involved changing some paths in sauger_cross_sim.R and suager_cross_sim.pl. The r script carries out simulating genotype matrices for each generation/group and creating a table of true parents. Then it uses sequoia::GetMaybeRel() to assign parent-offspring pairs and trios. The perl script runs all of those steps in the R script over all the combinations for set of values for nloci and wild_samp, where wild_samp is the percentage of the total simulated population that has been sampled. 

The scripts sauger_cross_check_trios.R and sauger_cross_check_pars.R then check the relationships that were extracted by GetMaybeRel() against the known parents from "true_parents" generated by sauger_cross_sim.R. 

These scripts appear to generate a table that describes how many true and false assignments were made for each combination of nloci and wild_samp.

Both sauger_cross_check_trios.R and sauger_cross_check_pars.R ran successfully, so I'm now curious about adding in missing data.

## Adding missing Data

sauger_cross_sim_sj.pl and sauger_cross_sim_sj.R were created on 04/11/25 to run the simulations on 125 combinations of nloci, wild_samp, and md values.

This was not completed successfully, however, until 04/13/25. Now, in

```{bash}
/project/ysctrout/hatchsauger/SaugerParentage/simulations
```

many text files exist for all combinations. Next step will be to run sauger_cross_check_pars.R and sauger_cross_check_trios.R to understand the number of true and false assignments that were made in each situation. I think I'd like Josh's assistance with this.


This was completed late April 2025 (around and on 4/24/25). See 
 ```{bash}
/project/ysctrout/hatchsauger/SaugerParentage/simulations/sims_with_md
 ```

for visualizations. Keep in mind that the simulations will need to be slightly altered for the correct PTPS scenario. At least for now, results indicate that missing data per site does not have a huge impact. Our ability to detect relationships depends mostly on the number of loci, and the PTPS.

# Bedtools Intersect

When it comes to comparing aln/samse and mem, we've noticed a strange trend that mem adds a TON of reads (many thousands), yet only adds about a thousand sites when aligned to the walleye reference. This is troubling because we're unsure as to where all of those other reads are going, and as to whether or not there is significant overlap in the sites that are being generated by the data when aligned using each algorithm. 

Bedtools is a "swiss army knife" for analyzing genetic data, and has an interesect function that will tell the locations and counts of sites where two files overlap. I created slurm_bedtools_int.sh to investigate. It appears as though .vcf is not a native format for bedtools, and it would easier digest a .bed file. I obtained this line to convert.

```{bash}
cd /project/ysctrout/hatchsauger/sam_sai_svit
bcftools query -f '%CHROM\t%POS0\t%POS\t%ID\n' rehead_variants_rawfiltered_svit_020625.vcf > svit_alnsamse.bed
mv svit_alnsamse.bed ../bedtools/
wc -l svit_mem.bed
# 404644 svit_mem.bed

cd /project/ysctrout/hatchsauger/sam_sai_svit_mem
bcftools query -f '%CHROM\t%POS0\t%POS\t%ID\n' rehead_variants_rawfiltered_svit_mem_031325.vcf > svit_mem.bed
mv svit_mem.bed ../bedtools/
wc -l svit_alnsamse.bed
# 366797 svit_alnsamse.bed
```

Alright. Information is being stored effectively in those .bed files. I'll now try an updated slurm_bedtools_int.sh.

```{bash}
wc -l intersect_output_unfiltered.bed 
# 333134 intersect_output_unfiltered.bed
```

I suppose I can interpret that as: mem is retaining most of the aln/samse snps since the number of snps in the output is pretty close to the number of snps from the aln/samse input, but mem is just grabbing a bunch more snps that aln/samse was not catching.

Per Josh's recommendation, I am running intersect on two filtered (MAF 0.01 and Miss 0.9) .vcfs (for aln/samse mem).
Grabbed filtered .vcfs from sam_sai_svit and sam_sai_svit_mem, copied them to the bedtools directory, and renamed them svit_alnsamse_variants_maf1_miss9.recode.vcf and  svit_mem_variants_maf1_miss9.recode.vcf. 

```{bash}
cd /project/ysctrout/hatchsauger/bedtools
bcftools query -f '%CHROM\t%POS0\t%POS\t%ID\n' svit_mem_variants_maf1_miss9.recode.vcf > svit_mem_maf1_miss9.bed
wc -l svit_mem_maf1_miss9.bed
# 7390 svit_mem_maf1_miss9.bed

cd /project/ysctrout/hatchsauger/bedtools
bcftools query -f '%CHROM\t%POS0\t%POS\t%ID\n' svit_alnsamse_variants_maf1_miss9.recode.vcf > svit_alnsamse_maf1_miss9.bed
wc -l svit_alnsamse_maf1_miss9.bed
# 6488 svit_alnsamse_maf1_miss9.bed
```

Alright. Information is being stored effectively in those .bed files. I'll now try an updated slurm_bedtools_int.sh.

```{bash}
wc -l intersect_output_maf1_miss9.bed 
# 6109 intersect_output_maf1_miss9.bed
```
Intersect shows overlap in ~94% of aln/samse sites. Looks like mem captures most of the same sites as aln/samse, but also picks up a bunch more.

# HipHop Parentage Assignment (+randomForest)

In the few weeks before WDAFS 2025, I was becoming pretty nervous about the state of my results using Sequoia, and sought out to make parentage relationships using the r package hiphop. This is an exclusion method that relies on principals of Mendelian inheritance to give hothiphop (Homozygous Opposite Test, Homozygous Identical Parents, Heterozygous Offspring Precluded) scores to each pair of parents that get assigned to each offspring. 

### There are a series of five r scripts in:
```
/project/ysctrout/hatchsauger/SaugerParentage/r_scripts/hiphop
```

First, I ran hiphop on all of the test F1's and F0 parents and got some pretty good results. I took the parent pairs that were inferred by hiphop and cross checked them against the list of known crosses. When I did that, 85 of the 95 test F1 offspring got assigned to parent crosses that were actually made.

Second, I wanted to see if I could improve that number with more intense filtering of the data. This is the second script. I got up to 92/95 with this "hard" filtering approach. See the filtering approaches in 

```{bash}
/project/ysctrout/hatchsauger/SaugerParentage/slurms/slurm_2softhiphopfilter.sh
```
and 

```{bash}
/project/ysctrout/hatchsauger/SaugerParentage/slurms/slurm_3hardhiphopfilter.sh
```

Third and fourth scripts are exploring the distributions of hothiphop scores when 75%, 50% and 25% of the true F0 parents are sampled. Breaks it down by PTPS, as well as if we have 0, 1, or 2 of the parents sampled for a given individual. Keep in mind that this is only one replication of these operations. We'll need to go back and rereun all of this assignment process with numerous replicates of sampling 75%, 50% and 25% of the true parents.

But, I was in a rush for this conference, and I wanted to see what the scores looked like for all of the F1 generation when they were assigned to ALL of the F0's, both the 2015 and 2016 groups. This is the fifth script.

All of that information was then plotted using
```{bash}
/Users/samjohnson/Documents/Sauger_042225/SaugerParentage/r_scripts/hiphop/sauger_hiphop_plot.R
```

## randomForest

During the week of 060225 I was also interested in running hiphop on the simulated data, though I learned that hiphop requires the parental sexes to be specified in the individuals dataframe. This thwarts our ability to use hiphop on the F2 -> F1 assignments, which has caused me to return to sequoia as the primary candidate for my parentage assignment method.

However, if we could determine sex with reasonable accuracy from the genetic data, as Will did, we'd be in business. I therefore put a yellow perch aligned and filtered .vcf -> genotype matrix into R and cbind()ed it with a df that contained all of the F0's and their sex data. 

```{bash}
/project/ysctrout/hatchsauger/sam_sai_pflav_mem_t2/rf_2/hard_variants_pflav_mem_t2_bial_noindels_q20_mindep3_maxdep75_maf1_miss95.012
```

I then used randomForest to see if these few sites (less than 200) could predict sex with any reasonable accuracy. It cannot. This was done on 07/02/25. (see SaugerParentage/r_scripts/randomForest)

This prompted me to begin to investigate Will's data to see where in the process my data are faltering. Why, when we align to the same reference, do mine perform so poorly by comparision?

# Contaminant Filtering
These developments and continued skepticism about my data have prompted me to attempt contaminant filtering to see if there is a larger problem with the data than anticipated. Contaminant filtering using tapioca was run with the following script, and relies on a few downloaded perl functions and downloaded modules (e.g., File::which). This module can be downloaded using cpanm. First, you must download cpanm to your home directory

```{bash}
cd /home/sjohn208/perl5/bin
curl -L https://cpanmin.us | perl - --local-lib=~/perl5 App::cpanminus
cpanm File::Which
```

and test with

``` {bash}
perl -MFile::Which -e 'print File::Which::which("perl"), "\n";'
```

Then you must set the environment variables in the contaminants script (see line 23 of the slurm script).

``` {bash}
cd /project/ysctrout/hatchsauger/SaugerParentage/slurms/slurm_sauger_contaminants.sh
```

Josh downloaded the illumina, phix, and ecoli files to /project/ysctrout/contaminants and the script goes through all of your raw reads, checks each of them against those databases, and spits out big cleaned .fastq's.

Contaminant filtering was initiated on 7/02/25 at 8:49pm.

Contaminant filtering took out a few of our reads, but allowed me to call more variants than beofre. See these notes from Sauger Meeting Notes F24/S25/Summer25.

"Got contam filtering to finish early Thursday morning, 8/7. Can’t believe it.
Splitting, parsing, concatenating, splitting to individual fastq’s finished 8/8.
Fastp filtering completed 8/8.
Alignment pending EOD 8/8. Fingers crossed.
Amazingly, I called one hundred thousand additional variants than past attempts using the walleye reference and bwa mem (496,823 variants for this round). I suppose that cleaning the dataset more vigorously allowed the algorithm to more effectively search through the remaining reads?"

Every dataset moving forward relied on this filtering. See:

```
cd /project/ysctrout/hatchsauger/sam_sai_contam_fastp_svit_mem
```

## Fastp Filtering

After contaminant filtering was completed, I filtered the raw reads using Fastp (MPR used it and recommended it). This code takes the raw reads and remove any bases with q<20, reads with mean quality less than 20, reads shorter than 50bp, poly g tails, and any remaining sequence adapters.

```
cd /project/ysctrout/hatchsauger/SaugerParentage/slurms/slurm_fastp.sh
```

After running this in early August 2026, all datasets included this filtering. See:

```
cd /project/ysctrout/hatchsauger/sam_sai_contam_fastp_svit_mem
```


# sam_sai_pflav_mem trials and updated location for sam_sai_svit_mem .vcfs
These developments and continued skepticism about my data have ALSO prompted me to try re-aligning and re-filtering the my data to the yellow perch genome using bwa mem. I was concerned about the file basenames having some issues during the first pflav_mem trial (e.g., some files were yellowperch.genome.fna instead of just yellowperch.fna). sam_sai_pflav_mem_t2 was created to address that issue. The alignment scripts (.pl and .sh) were altered and the alignment and variant calling were rerun but with no success.

In contrast, sam_sai_pflav_mem_t3 was created because I saw some pflav .vcfs in sam_sai_svit_mem/Sequoia_Inp/maf30_miss95. First I redid the filtering for svit_mem (with thinning) and placed updated walleye aligned .vcf's in: 

```{bash}
/project/ysctrout/hatchsauger/sam_sai_svit_mem/thin_filtered
```

Then, I created the sam_sai_pflav_mem_t3 to redo the original sam_sai_pflav_mem filtering. That is now updated with thinning filters and should be used from now on. This also prompted me to check on the sam_sai_pflav_mem_t2 files. These t2 files were also freshly filtered on 07/02/25, and were what I used for randomForest on this day (see above).

Counting SNP's became very important for all of this, so to remain consistent, use this command to count the number of SNP's in a .vcf. 

```{bash}
grep -c "^scaff" filename.vcf or *.vcf
```






