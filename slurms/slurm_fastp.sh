#!/bin/bash

#SBATCH --job-name=FastP
#SBATCH -A ysctrout
#SBATCH -t 1-12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10   
#SBATCH --mem=50G
#SBATCH --mail-type=END
#SBATCH --array=0-1183

# this script is by MPR. working to edit my own. added to git 070825
# filters split fastq's for: 
# 1. Trimming bases lower than 20 
# 2. Reads shorter than 50 bp 
# 3. Remove poly g tails 
# 4. Remove any sequence adapters (checking against that db she’s got in there). 

# SPJ edits 070825
# USAGE: sbatch slurm_fastp.sh /path/to/split/fastqs /path/to/cleaned/fastqs
# today's usage: sbatch slurm_fastp.sh /project/ysctrout/hatchsauger/1Saug/rawfastqs /project/ysctrout/hatchsauger/1Saug/fastp_cleaned

#SPJ edits 080825
# today's usage: sbatch slurm_fastp.sh /project/ysctrout/hatchsauger/1Saug/contam_cleaned/ind_files /project/ysctrout/hatchsauger/1Saug/contam_cleaned/fastp_cleaned

# Load fastp module
module load fastp/0.23.4

# Define paths
#RAW_DIR="/project/ysctrout/mrodri23/genetic_data/raw_reads/2024_data/1/Prep2/demultiplexinglib1/split_fastqs/"
RAW_DIR=$1
#OUTDIR="/project/ysctrout/mrodri23/genetic_data/raw_reads/2024_data/1/Prep2/demultiplexinglib1/split_fastqs/trim_fastq"
OUTDIR=$2
#REPORTDIR="/project/ysctrout/mrodri23/genetic_data/raw_reads/2024_data/1/Prep2/demultiplexinglib1/split_fastqs/fastp_reports"
#REPORTDIR=$3
ADAPTERS="/project/ysctrout/mrodri23/scripts/rawreads_quality/adapters.fa"


# Get all samples
cd "$RAW_DIR"
#dirs=(*fastq.gz)
dirs=(*.fastq)
#dirs=(*SLEI15C*)
#dirs=(*BTCH20C*)

# Current sample 
sample_dir="${dirs[$SLURM_ARRAY_TASK_ID]}"


# Remove the .fastq.gz. Just needed for the reports
#SAMPLE=${sample_dir%.fastq.gz}
SAMPLE=${sample_dir%.fastq}
#REPORT_JSON="${REPORTDIR}/${SAMPLE}_fastp.json" # commenting otu when not wanting to save reports
#REPORT_HTML="${REPORTDIR}/${SAMPLE}_fastp.html"

# Define the name of output trimmed reads
OUT_SAMPLE="${OUTDIR}/trim_read_${sample_dir}"


# Run fastp

# -Q or --disable_quality_filtering  to disable quality filter per read
#fastp -i "$sample_dir" -o "$OUT_SAMPLE" -h "$REPORT_HTML" -j "$REPORT_JSON" -e 20  --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 20  --length_required 50 --trim_poly_g  -w "$SLURM_CPUS_PER_TASK" --adapter_fasta "$ADAPTERS" 

fastp -i "$sample_dir" -o "$OUT_SAMPLE" -e 20  --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 20  --length_required 50 --trim_poly_g  -w "$SLURM_CPUS_PER_TASK" --adapter_fasta "$ADAPTERS"
