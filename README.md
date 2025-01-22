# SaugerParentage
Created by Sam Johnson, 10/25/24

sjohn208@uwyo.edu

## Project Summary

Below is the code and notes for the Wind River Sauger Parentage Analysis bioinformatics workflow. This began on 10/25/24 with file transfers and establishment of the Sauger Parentage github repository, followed by unzipping files (10/29/24) and parsing (11/01/24).

## Workflow Outline

* Demultiplex

   * Count Raw Reads

   * Splitting Raw Files

   * Parsing

   * Debugging perl on Medbow

* Split .fastq

* Alignment

* Variant Calling

* Filtering

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

## Alignment

### Yellow Perch Reference
First step is to remove the line endings after 60 characters, using remove60_fasta.pl
```{bash}
## remove 60 line endings
perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/remove60_fasta.pl Perca_flavescens.fasta
```
Second step is to rename the reference names into something a bit more meaningful/reasonable.
```{bash}
## change reference names to something shorter

perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/rename_scaff.pl no60_Perca_flavescens.fasta

mv renamed_no60_Perca_flavescens.fasta.txt yellowperch_genome.fna
```
Third step is to make an index for the reference genome bwa (Burrow Wheeler Align).
```{bash}
sbatch slurm_sauger_bwa_index.sh
```


Run Burrow Wheeler Align on /project/ysctrout/hatchsauger/1Saug/rawfastqs/*.fastq

Generates for you a sam_sai directory that houses the output from this step, aligned reads with assoc. alignment scores. For each individual's split .fastq, you get a .sam and .sai file (standard alignment map, index, respectively). Can less into .sam but not .sai
```{bash}
sbatch slurm_sauger_runbwa.sh
```