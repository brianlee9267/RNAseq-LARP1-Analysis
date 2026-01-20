#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 24:00:00
#SBATCH --mem=3GB
#SBATCH --chdir=[PATH_TO_PROJECT_DIR]
#SBATCH -o "%x.o%j.log"
#SBATCH --nodes=1
#SBATCH -n 4
#SBATCH --job-name=phs002121.RNAseq
#SBATCH --mail-user=[YOUR_EMAIL]@northwestern.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Create sample sheet for RNAseq pipeline
# Make executable with chmod +x fastq_dir_to_samplesheet.py
# Sample info: AMPK_DK LN samples rerun

./fastq_dir_to_samplesheet.py ./ samplesheet.csv --read1_extension .fastq.gz --read2_extension .fastq.gz --strandedness reverse
