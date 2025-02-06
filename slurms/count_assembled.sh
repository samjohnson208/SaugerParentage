#!/bin/sh

## Run this on Teton inside assem directory (the one that contains bamfiles) to count number of reads assembled and number of raw reads per individuals with samtools

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
raw=`samtools stats $file | grep "raw total sequences:" | sed 's/SN\t.*:\t//g'`
assembled=`samtools stats $file | grep "reads mapped:" | sed 's/SN\t.*:\t//g'`
echo "$indname $raw $assembled" >> assembled_per_ind.txt
done