#!/bin/sh

## slurm_sfs.sh by SPJ on 021925
## PURPOSE: to construct site frequency spectra for .vcf files
## USAGE: sbatch slurm_sfs.sh

#SBATCH --job-name=easysfs
#SBATCH --account=ysctrout
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --mail-type=END
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_easysfs
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_easysfs

## load modules, activate and install conda environment, pull in git repo, and change 
module load arcc/1.0 miniconda3/24.3.0
conda activate easySFS

# navigate to easySFS directory and change script permissions
cd /project/ysctrout/hatchsauger/easySFS
chmod 777 easySFS.py

## choose directory that holds the vcf upon which the script will act
cd /project/ysctrout/hatchsauger/sam_sai_pflav/

## run for that .vcf
/project/ysctrout/hatchsauger/easySFS/easySFS.py -i /project/ysctrout/hatchsauger/sam_sai_pflav/rehead_variants_rawfiltered_pflav_012325.vcf -p /project/ysctrout/hatchsauger/sam_sai_pflav/sauger_ids_pop.txt --proj 10,50,100,500 -o sfs_out_pflav


## choose directory that holds the vcf upon which the script will act
cd /project/ysctrout/hatchsauger/sam_sai_svit/

## run for that .vcf
/project/ysctrout/hatchsauger/easySFS/easySFS.py -i /project/ysctrout/hatchsauger/sam_sai_svit/rehead_variants_rawfiltered_svit_020625.vcf -p /project/ysctrout/hatchsauger/sam_sai_svit/sauger_ids_pop.txt --proj 10,50,100,500 -o sfs_out_svit


