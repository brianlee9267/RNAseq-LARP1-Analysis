# RNA-seq Analysis Pipeline: LARP1KO vs CaS9

This repository contains the code for processing and analyzing RNA-seq data derived from MPN cell lines (HEL/SET-2) to investigate the effects of LARP1KO.

## Workflow Overview

The analysis follows a three-step workflow:

### 1. Sample Sheet Generation
**Script:** `01_generate_samplesheet.sh`
* **Function:** Generates the input CSV required by the Nextflow pipeline.
* **Input:** Directory containing raw `.fastq.gz` files.
* **Output:** `samplesheet.csv`

### 2. Pipeline Execution (Nextflow)
**Script:** `02_run_rnaseq.sh`
* **Function:** Runs the `nf-core/rnaseq` pipeline using Singularity on the compute cluster.
* **Process:** Aligns reads to the `hg38` genome and quantifies gene expression.
* **Output:** BAM files and Salmon gene counts.

### 3. Downstream Analysis (R/DESeq2)
**Script:** `03_deseq2_analysis.R`
* **Function:** Performs differential expression analysis, PCA, and functional enrichment (GO/KEGG/DO).
* **Key Comparison:** LARP1KO vs Cas9 Control.
* **Output:** Volcano plots, Heatmaps, and lists of differentially expressed genes (DEGs).

## Prerequisites

* **Compute Environment:** Linux cluster with SLURM scheduler.
* **Software:**
    * Nextflow
    * Singularity
    * R (v4.0+) with Bioconductor packages (`DESeq2`, `clusterProfiler`, `ComplexHeatmap`).

## Usage Notes

* Ensure paths in the `.sh` scripts are updated to match your directory structure before running.
* The R script requires the `salmon.merged.gene_counts_polysomal_KO.csv` and metadata file `sampleClassificationPolysomalKO.csv` to be present in the working directory.
