#!/bin/bash

## slurm_sauger_contaminants by JPJ and SPJ 063025
## PURPOSE: to filter raw reads for contaminants
## USAGE: sbatch slurm_sauger_contaminants.sh

#SBATCH --job-name=contams
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=ysctrout
#SBATCH --time=6-23:59:59
#SBATCH --mem=10G
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_contam
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_contam
#SBATCH --mail-type=END

module load arcc/1.0 bowtie2/2.5.4 perl/5.38.0_x86_64

##########################################################################################
## odds library
##########################################################################################

/project/ysctrout/bin/tapioca/src/tap_contam_analysis --db /project/ysctrout/contaminants/illumina_oligos --pct 20 /project/ysctrout/hatchsauger/1Saug/rawreads/1SaugOdds.fastq > 1SaugOdds.readstofilter.ill.txt 

echo "Illumina filtering done for odds"

/project/ysctrout/bin/tapioca/src/tap_contam_analysis --db /project/ysctrout/contaminants/phix174 --pct 80 /project/ysctrout/hatchsauger/1Saug/rawreads/1SaugOdds.fastq > 1SaugOdds.readstofilter.phix.txt 

echo "PhiX filtering done for odds"


/project/ysctrout/bin/tapioca/src/tap_contam_analysis --db /project/ysctrout/contaminants/ecoli-k-12 --pct 80 /project/ysctrout/hatchsauger/1Saug/rawreads/1SaugOdds.fastq > 1SaugOdds.readstofilter.ecoli.txt

echo "ecoli filtering done for odds"


cat /project/ysctrout/hatchsauger/1Saug/rawreads/1SaugOdds.fastq | /project/ysctrout/bin/fqu_cull -r 1SaugOdds.readstofilter.ill.txt 1SaugOdds.readstofilter.phix.txt 1SaugOdds.readstofilter.ecoli.txt > 1SaugOdds.clean.fastq

echo "Clean copy of odds done"


##########################################################################################
## evens library
##########################################################################################

/project/ysctrout/bin/tapioca/src/tap_contam_analysis --db /project/ysctrout/contaminants/illumina_oligos --pct 20 /project/ysctrout/hatchsauger/1Saug/rawreads/1SaugEvens.fastq > 1SaugEvens.readstofilter.ill.txt 

echo "Illumina filtering done for evens"

/project/ysctrout/bin/tapioca/src/tap_contam_analysis --db /project/ysctrout/contaminants/phix174 --pct 80 /project/ysctrout/hatchsauger/1Saug/rawreads/1SaugEvens.fastq > 1SaugEvens.readstofilter.phix.txt 

echo "PhiX filtering done for evens"


/project/ysctrout/bin/tapioca/src/tap_contam_analysis --db /project/ysctrout/contaminants/ecoli-k-12 --pct 80 /project/ysctrout/hatchsauger/1Saug/rawreads/1SaugEvens.fastq > 1SaugEvens.readstofilter.ecoli.txt

echo "ecoli filtering done for evens"


cat /project/ysctrout/hatchsauger/1Saug/rawreads/1SaugEvens.fastq | /project/ysctrout/bin/fqu_cull -r 1SaugEvens.readstofilter.ill.txt 1SaugEvens.readstofilter.phix.txt 1SaugEvens.readstofilter.ecoli.txt > 1SaugEvens.clean.fastq

echo "Clean copy of evens done"




