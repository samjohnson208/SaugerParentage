cd # SaugerParentage
Created by Sam Johnson, 10/25/24

sjohn208@uwyo.edu

## Project Summary

Below is the code and notes for the Wind River Sauger Parentage Analysis bioinformatics workflow. This began on 10/25/24 with file transfers and establishment of the Sauger Parentage github repository, followed by unzipping files (10/29/24) and parsing (11/01/24).

Problems arose in the alignment phase, where vastly different datasets were produced by aligning to either the yellow perch or walleye reference genome. Information regarding these alignments are deonted by: 
Trial 1: Alignment to the Yellow Perch genome (Perca flavescens) produced very few loci after initial filtering. 
Trial 2: Alignment to the Walleye refernce genome (Sander vitreus)

## Workflow Outline

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

* Alignment (Yellow Perch and Walleye, Trials 1 and 2, respectively)

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

   * Variant Calling BWA MEM, Reheadring

   * SFS Analysis for all combinations of reference genomes and alignment algorithms

* Sequoia

   * Generating input for Sequoia
   
   * Sequoia using F0's & Test F1's (aln/samse and mem)

* Simulation Work

* Bedtools Intersect (aln/samse and mem overlap)

* Entropy

## Demultiplex

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

## Split .fastq
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

## Alignment

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

## BWA MEM Investigation

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

03/14/25 - All four directories aln/mem pflav/svit now have rawfiltered vcfs in them. 

I want to take each of those vcfs and filter them for MAF 0.01 and Miss 0.9. Then, turn those into genotype  matrices, and plot an SFS for each that I can send to them in a 2x2 panel.

First step is to create reheadered vcfs for both mem directories:

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

Now that we have filtered vcfs to MAF 0.01 and Miss 0.90 for all four combinations, let's convert them to genotype matrices and prepare for SFS. I'm also going to convert all of them for Sequoia, since I imagine I'll put them in there at some point. Here's the example for pflav_mem

```{bash}
salloc --account=ysctrout --time=3:00:00
module load arcc/1.0 gcc/14.2.0 vcftools/0.1.17
cd /project/ysctrout/hatchsauger/sam_sai_pflav_mem/
vcftools --vcf variants_maf1_miss9.recode.vcf --012 --out variants_maf1_miss9_mem_pflav
sed "s/-1/-9/g" variants_maf1_miss9_mem_pflav.012 > variants_maf1_miss9_mem_pflav.012_conv
# conv = converted to missing data = "-9"
```

SFS generated in SaugerParentage/r_scripts/sfs_all.R on 03/14/25


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

Sequoia recommends more than 90% genotype call rate, and a MAF of 30% or higher. I have already filtered the rawfiltered vcf,  rehead_variants_rawfiltered_svit_020625.vcf, using: 

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

No relationships found on that data, perhaps because we filtered with maf of 0.03 above instead of 0.3. Ugh. Resume after break and try again. The filtered vcf is already there in sam_sai_svit.

03/24/25 Work: Took the maf30_miss9 vcf and turned it into the genotype matrix. prepared Sequoia_Full.R for a run on that. There are still 1400 sites. May be too many still. We shall see.

Results: Only one parent offspring duo was produced by this run.

03/26/25 Work: Trying maf30_misss9 on svit_mem. Prepared and ran Sequoia_Full.R on that converted matrix. 1540 sites. This should let us know if the number of sites is the problem, or if it's an alignment issue.

Results:

## A Note on Nomenclature

03/14/25 - Changed all of the .bam, .sorted.bam, and .sorted.bam.bai in both mem directories to start with the prefix "mem_" rather than "aln_". The runbwa_mem.pl script was also updated to reflect this change, should we have to run that alignment again. Keep that in mind for any discrepancies in documentation you may come across down the line.

## Simulation Work

04/09/25 - Today is the firs time that I pulled down Josh's scripts from github and ran them in full, on the cluster. This involved changing some paths in sauger_cross_sim.R and suager_cross_sim.pl. The r script carries out simulating genotype matrices for each generation/group and creating a table of true parents. Then it uses sequoia::GetMaybeRel() to assign parent-offspring pairs and trios. The perl script runs all of those steps in the R script over all the combinations for set of values for nloci and wild_samp, where wild_samp is the percentage of the total simulated population that has been sampled. 

The scripts sauger_cross_check_trios.R and sauger_cross_check_pars.R then check the relationships that were extracted by GetMaybeRel() against the known parents from "true_parents" generated by sauger_cross_sim.R. 

These scripts appear to generate a table that describes how many true and false assignments were made for each combination of nloci and wild_samp.

Both sauger_cross_check_trios.R and sauger_cross_check_pars.R ran successfully, so I'm now curious about adding in missing data.

#### Adding missing Data

sauger_cross_sim_sj.pl and sauger_cross_sim_sj.R were created on 04/11/25 to run the simulations on 125 combinations of nloci, wild_samp, and md values.

This was not completed successfully, however, until 04/13/25. Now, in

```{bash}
/project/ysctrout/hatchsauger/SaugerParentage/simulations
```

many text files exist for all combinations. Next step will be to run sauger_cross_check_pars.R and sauger_cross_check_trios.R to understand the number of true and false assignments that were made in each situation. I think I'd like Josh's assistance with this.


## Bedtools Intersect

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

Per Josh's recommendation, I am running intersect on two filtered (MAF 0.01 and Miss 0.9) vcfs (for aln/samse mem).
Grabbed filtered vcfs from sam_sai_svit and sam_sai_svit_mem, copied them to the bedtools directory, and renamed them svit_alnsamse_variants_maf1_miss9.recode.vcf and  svit_mem_variants_maf1_miss9.recode.vcf. 

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




