# CRISPER_DEG-50GENES
## Differential Gene Expression Analysis of TCGA Data
This project performs differential gene expression analysis on RNA-seq data from The Cancer Genome Atlas (TCGA), comparing tumor samples to normal samples.

### Overview
Data preprocessing and filtering using edgeR
Normalization of gene expression data
Differential expression analysis
Visualization of results including MD plots and volcano plots
Generation of heatmaps for top differentially expressed genes

### Key Steps
Created a DGEList object from TCGA count data
Filtered low-expression genes using filterByExpr()
Normalized the data with calcNormFactors()
Estimated dispersion and fit a quasi-likelihood negative binomial model
Performed differential expression testing with glmQLFTest()
Identified significantly up- and down-regulated genes
Created visualizations including MD plots, volcano plots, and heatmaps

### Results
Identified 4,351 up-regulated and 2,291 down-regulated genes in tumor samples
Generated a heatmap of the top 50 differentially expressed genes
Saved full results to "differential_expression_results.csv"

#### Files
script.R: R script containing the full analysis pipeline
differential_expression_results.csv: CSV file with all differentially expressed genes
heatmap_top50_DEGs.png: Heatmap image of top 50 differentially expressed genes

### Dependencies
R (version used: [Your R version])
Bioconductor packages: edgeR, limma
Other R packages: pheatmap
