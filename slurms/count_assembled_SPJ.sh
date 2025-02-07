#!/bin/sh

## Run this on Teton inside assem directory (the one that contains bamfiles) to count number of reads assembled and number of raw reads per individuals with samtools
## USAGE: sbatch /project/ysctrout/hatchsauger/SaugerParentage/slurms/count_assembled_SPJ.sh

#SBATCH --account=ysctrout
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5G
#SBATCH --mail-type=END

module load swset
module load gcc/7.3.0
module load samtools/1.8
echo "ind raw assembled" > assembled_per_ind.txt

for file in aln_*.sorted.bam
do
  indname=`echo $file | sed 's/aln_//g' | sed 's/\.sorted\.bam//g'`
  
  # Extract raw total sequences
  raw=$(samtools stats $file | grep "^SN\ttotal sequences:" | cut -f 2)
  
  # Extract mapped (assembled) reads
  assembled=$(samtools stats $file | grep "^SN\treads mapped:" | cut -f 2)
  
  # Write results to output file
  echo "$indname $raw $assembled" >> assembled_per_ind.txt
done
