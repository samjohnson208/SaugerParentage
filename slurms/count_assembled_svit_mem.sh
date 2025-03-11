#!/bin/sh

## Run this on medbow inside assem directory (the one that contains bamfiles) to count number of reads assembled and number of raw reads per individuals with samtools
## USAGE: sbatch /project/ysctrout/hatchsauger/SaugerParentage/slurms/count_assembled.sh

#SBATCH --account=ysctrout
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --mail-type=END
#SBATCH -o /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stdout_count_svit_mem
#SBATCH -e /project/ysctrout/hatchsauger/SaugerParentage/slurms/std/stderr_count_svit_mem

module load gcc/14.2.0
module load samtools/1.20
echo "ind raw assembled" > assembled_per_ind_svit_mem.txt

for file in aln_*.sorted.bam
do
indname=`echo $file | sed 's/aln_//g' | sed 's/\.sorted\.bam//g'`
raw=`samtools stats $file | grep "raw total sequences:" | sed 's/SN\t.*:\t//g'`
assembled=`samtools stats $file | grep "reads mapped:" | sed 's/SN\t.*:\t//g'`
echo "$indname $raw $assembled" >> assembled_per_ind_svit_mem.txt
done