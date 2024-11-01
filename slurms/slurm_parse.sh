#!/bin/bash

## slurm_parse.sh by SPJ 110124
## PURPOSE: to parse raw .fastq files
## USAGE: sbatch slurm_parse.sh

#SBATCH --job-name=parse
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=ysctrout
#SBATCH --time=1-00:00:00
#SBATCH --mem=1G

module load arcc/1.0 gcc/12.2.0 perl/5.34.1

cd /project/ysctrout/hatchsauger/1Saug/rawreads
perl /project/ysctrout/hatchsauger/SaugerParentage/perl_scripts/parse_barcodes768.pl 1SaugEvens_Demux.csv 1SaugEvens.fastq LH0
