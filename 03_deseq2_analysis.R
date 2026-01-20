# ==========================================
# RNA-seq Downstream Analysis: LARP1KO vs Cas9
# Input files required:
# 1. salmon.merged.gene_counts_polysomal_KO.csv
# 2. sampleClassificationPolysomalKO.csv
# ==========================================

# Install necessary packages if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

packages <- c("ComplexHeatmap", "clusterProfiler", "org.Hs.eg.db", 
              "AnnotationDbi", "enrichplot", "DESeq2", "DOSE")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}

library(tidyverse)
library(circlize)
library(pheatmap)

# Load your data
counts <- read.csv('salmon.merged.gene_counts_polysomal_KO.csv', row.names = 1)
colData <- read.csv('sampleClassificationPolysomalKO.csv')

# Filter colData to only include rows for which SRR is in the column names of counts
colData_filtered <- colData[colData$ID %in% colnames(counts), ]
# Set the row names of colData_filtered to the SRR column values
rownames(colData_filtered) <- colData_filtered$ID

# Verify alignment
if (ncol(counts) == nrow(colData_filtered)) {
  print("The datasets are now compatible for DESeq2 analysis.")
} else {
  stop("Mismatch between counts and metadata dimensions.")
}

# Round counts to integers for DESeq2
counts_integer <- round(counts)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_integer,
                              colData = colData_filtered,
                              design = ~ Treatment)

# Run DESeq
dds <- DESeq(dds)

# --- PCA and QC ---
vsd <- vst(dds)
pca_res <- prcomp(t(assay(vsd)))
plotPCA(vsd, intgroup = c("Treatment"))

# --- Differential Expression Analysis ---
# Results extraction
res <- results(dds)

# Filter DEGs (padj < 0.05)
degs <- rownames(res[!is.na(res$padj) & res$padj < 0.05, ])

# --- Heatmap Visualization ---
norm_data <- assay(vsd)
annotation_for_heatmap <- data.frame(Treatment = colData_filtered$Treatment)
rownames(annotation_for_heatmap) <- rownames(colData_filtered)

# Subset to DEGs
norm_data_degs <- norm_data[degs, ]
final_annotation <- annotation_for_heatmap[colnames(norm_data_degs), , drop = FALSE]

if(length(degs) > 0) {
  pheatmap::pheatmap(norm_data_degs,
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     scale = "row",
                     show_rownames = FALSE,
                     show_colnames = TRUE,
                     color = colorRampPalette(c("blue", "white", "red"))(255),
                     annotation_col = final_annotation)
} else {
  print("No differentially expressed genes found.")
}

# --- Functional Enrichment (GO/KEGG) ---
# Contrast: LARP1KO vs Cas9
res_CTRL_vs_MF <- results(dds, contrast=c("Treatment", "LARP1KO", "Cas9"))
degs_CTRL_vs_MF <- rownames(res_CTRL_vs_MF[which(res_CTRL_vs_MF$padj < 0.05), ])

# Map Symbols to Entrez IDs
entrezIds <- mapIds(org.Hs.eg.db,
                    keys = degs_CTRL_vs_MF,
                    column = "ENTREZID",
                    keytype = "SYMBOL",
                    multiVals = "first")
entrezIds <- na.omit(entrezIds)

# GO Enrichment
ego <- enrichGO(gene          = entrezIds,
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05,
                readable      = TRUE)

dotplot(ego, showCategory=20)
ego_sim <- pairwise_termsim(ego)
emapplot(ego_sim, showCategory=20)

# KEGG Enrichment
kegg_res <- enrichKEGG(gene = entrezIds,
                       organism = 'hsa',
                       pAdjustMethod = "BH",
                       qvalueCutoff = 0.05)
barplot(kegg_res)

# Disease Ontology
do_res <- enrichDO(gene = entrezIds,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05)
barplot(do_res)

# --- Export Results ---
res_Treatment_df <- as.data.frame(res_CTRL_vs_MF)
res_Treatment_df <- res_Treatment_df[!is.na(res_Treatment_df$padj) & res_Treatment_df$padj < 1, ]

final_Treatment_df <- data.frame(
  Gene = rownames(res_Treatment_df),
  log2FC = res_Treatment_df$log2FoldChange,
  padj = res_Treatment_df$padj
)

write.csv(final_Treatment_df, file = "degs_polysomal_LARP1KO_vs_Cas9_genes.csv", row.names = FALSE)
