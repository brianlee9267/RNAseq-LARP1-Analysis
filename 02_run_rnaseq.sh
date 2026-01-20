module load singularity
module load graphviz/2.40.1
module load java/jdk11.0.10

# Use screen to run this in the background:
# ctrl+a+d to return to terminal, screen -r to return to screen

# Update paths below before running
INPUT_SHEET="[PATH_TO_SAMPLESHEET]/samplesheet7.csv"
OUT_DIR="[PATH_TO_OUTPUT_DIR]"

nextflow run nf-core/rnaseq \
    --input $INPUT_SHEET \
    --outdir $OUT_DIR \
    --genome hg38 \
    -profile nu_genomics \
    --email [YOUR_EMAIL]@northwestern.edu
