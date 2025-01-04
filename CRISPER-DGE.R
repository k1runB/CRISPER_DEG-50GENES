
install.packages("BiocManager")

# install GEOquery
BiocManager::install("GEOquery")
library(GEOquery)

# Load the dataset
gse <- getGEO("GSE120861", GSEMatrix = TRUE)
data <- exprs(gse[[1]])

# View the data
head(data)
# Filter out low-count genes
filtered_data <- data[rowSums(data > 10) > 5, ]

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("locfit")
BiocManager::install("Rcpp")

library(edgeR)

# Create DGEList object
dge <- DGEList(counts = filtered_data)

# Normalize using TMM (Trimmed Mean of M-values)
dge <- calcNormFactors(dge)

dim(filtered_data)
summary(rowSums(data))
hist(log2(rowSums(data) + 1))
head(data)
str(data)
dim(data)

library(GEOquery)

# Load the dataset
gse <- getGEO("GSE120861", GSEMatrix = TRUE)
data <- exprs(gse[[1]])

# Check data characteristics
dim(data)
summary(data)
head(data)

# Verify if it's count data or normalized data
hist(log2(data[,1] + 1), main="Distribution of expression values")

# If it's normalized data, you might not need to filter as aggressively
# Adjust filtering if needed
filtered_data <- data[rowSums(data > 0) > 5, ]

# Create DGEList
library(edgeR)
dge <- DGEList(counts = filtered_data)

# Check the result
dim(dge$counts)
head(dge$counts)

library(GEOquery)

# Reload the dataset
gse <- getGEO("GSE120861", GSEMatrix = TRUE)

# Check what's in the gse object
length(gse)
class(gse[[1]])

# Try to access the expression data differently
if (is(gse[[1]], "ExpressionSet")) {
  data <- exprs(gse[[1]])
} else {
  print("Data is not in ExpressionSet format")
}

# Check the data
dim(data)
class(data)
head(data[,1:5])  # View first 5 columns of data

# If data is still not correct, try accessing it as a data frame
phenoData <- pData(gse[[1]])
head(phenoData)

# Check if the actual data is stored elsewhere in the gse object
names(gse[[1]])

# Convert to matrix if necessary
tpm_matrix <- as.matrix(tpm_data)

# Log-transform the data (common for TPM values)
log_tpm <- log2(tpm_matrix + 1)
library(edgeR)


# Assuming you already have your dge object with counts
BiocManager::install("edgeR", update = TRUE, ask = FALSE)

# Extract the necessary components
counts <- dge$counts
samples <- dge$samples
genes <- data.frame(genes = rownames(counts))

# Create a new DGEList object
new_dge <- DGEList(counts = counts, samples = samples, genes = genes)

# Verify the structure
str(new_dge)

# Try filtering again
keep <- filterByExpr(new_dge, group = new_dge$samples$group)
new_dge_filtered <- new_dge[keep, , keep.lib.sizes=FALSE]

# Check the dimensions of the filtered data
dim(new_dge_filtered)

# Assuming you have already created new_dge and performed filtering

# Normalize the data
new_dge_filtered <- calcNormFactors(new_dge_filtered)

# Set up the design matrix
design <- model.matrix(~new_dge_filtered$samples$group)

# Estimate dispersions
new_dge_filtered <- estimateDisp(new_dge_filtered, design)

# Fit the model
fit <- glmQLFit(new_dge_filtered, design)

# Perform differential expression analysis
qlf <- glmQLFTest(fit, coef=2)

# View top differentially expressed genes
topTags(qlf)

# Get all results
results <- topTags(qlf, n=Inf)

# Write results to a file
write.csv(results, file="differential_expression_results.csv")

# Create an MD plot
plotMD(qlf)

# Create a volcano plot
with(results$table, plot(logFC, -log10(FDR), pch=20, main="Volcano plot", xlim=c(-3,3)))
with(subset(results$table, FDR<0.05 & abs(logFC)>1), points(logFC, -log10(FDR), pch=20, col="red"))

# Identify significantly differentially expressed genes
sig_genes <- decideTestsDGE(qlf)
summary(sig_genes)

install.packages("pheatmap")

# Install and load required packages if not already done
if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")
library(pheatmap)

# After your existing code, add:

# Extract normalized log-CPM values
logCPM <- cpm(new_dge_filtered, log=TRUE)

# Get the top 50 differentially expressed genes
top_genes <- rownames(topTags(qlf, n=50))

# Subset the logCPM matrix for these genes
top_gene_matrix <- logCPM[top_genes,]

# Create a dataframe with sample information
sample_info <- data.frame(Group = new_dge_filtered$samples$group)
rownames(sample_info) <- colnames(top_gene_matrix)

# Generate heatmap
pheatmap(top_gene_matrix,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = sample_info,
         main = "Top 50 Differentially Expressed Genes")

# If you want to save the heatmap to a file:
png("heatmap_DEGs.png", width = 1200, height = 1000, res = 150)
pheatmap(top_gene_matrix,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = sample_info,
         main = "Top 50 Differentially Expressed Genes")
dev.off()

if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")
library(pheatmap)

