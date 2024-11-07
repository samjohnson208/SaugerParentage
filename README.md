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

   * Count mids

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

```{bash}
split -d -l 55000000 --verbose --additional-suffix \.fastq 1SaugEvens.fastq even_fq_
split -d -l 55000000 --verbose --additional-suffix \.fastq 1SaugOdds.fastq odd_fq_
```

### Parsing

Created new scripts to parallelize the parsing process across all split files in evens and odds that start with "even_fq" or "odd_fq". Those scripts are run_parse_evens.pl and run_parse_odds.pl (Located in perl_scripts directory)


```{bash}
perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/run_parse_evens.pl even_fq_*
perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/run_parse_odds.pl odd_fq_*
```

### Count mids

```{bash}
## how many good mids?

grep "Good mids" parsereport_* | cut -f 4 -d " " > good_mids_count.txt
   #for all of these parsed report files, take the fourth column of the line that has "Good mids"
module load arcc/1.0 gcc/12.2.0 r/4.4.0
R

dat <- read.delim("good_mids_count.txt", header = FALSE)
sum(dat[,1])
	## 166 files, 1,394,271,363 good mids

#exit R with 
quit()
```