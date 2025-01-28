#################################################################
##  Script: RNA_transcript_analysis.R
##  Step 1: Load and merge ENCODE transcript quantifications
##
##  Author: Daniel Batyrev
#################################################################

# Load required libraries
library(data.table)  # Required for fread
library(rtracklayer) # Required for reading GTF
library(DESeq2)      # Required for differential expression analysis
library(ggplot2)     # Required for plotting

# Set working directory
base_dir <- "D:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/RNA-seq"
setwd(base_dir)

# Create necessary folders if they don't exist
dirs <- c("raw_data2", "processed_data", "results", "plots")
for (d in dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Define file paths
metadata_file <- file.path(base_dir, "raw_data2", "metadata.tsv")
expression_output <- file.path(base_dir, "processed_data", "transcript_matrix.csv")

#################################################################
##  Step 1: Process Metadata
#################################################################

cat("Processing metadata...\n")

# Read metadata
metadata <- fread(metadata_file, sep = "\t", header = TRUE)

# Extract relevant columns
metadata_clean <- metadata[, .(
  File_accession = `File accession`,
  Experiment_accession = `Experiment accession`,
  Biosample = `Biosample term name`,
  Assay = Assay,
  Output_type = `Output type`,
  File_assembly = `File assembly`,
  File_status = `File Status`,
  File_format = `File format`,
  File_analysis = `File analysis title`
)]

# Save cleaned metadata
metadata_clean_file <- file.path(base_dir, "processed_data", "metadata_clean.csv")
fwrite(metadata_clean, metadata_clean_file)

cat("Metadata processing complete. Saved to:", metadata_clean_file, "\n")

#################################################################
##  Step 2: Merge Transcript Quantifications
#################################################################

cat("Merging transcript quantification files...\n")

# List all transcript quantification files (EXCLUDING metadata.tsv)
rna_files <- list.files(file.path(base_dir, "raw_data2"), pattern = "*.tsv", full.names = TRUE)
rna_files <- rna_files[!grepl("metadata.tsv", rna_files)]  # Exclude metadata file


read_expression_data <- function(file) {
  df <- fread(file, sep = "\t", header = TRUE)
  
  # Rename transcript_id(s) column to remove special characters
  setnames(df, old = "transcript_id(s)", new = "transcript_id")
  
  # Select only necessary columns
  df <- df[, .(gene_id, transcript_id, expected_count)]
  
  # Extract sample name
  sample_name <- tools::file_path_sans_ext(basename(file))
  setnames(df, "expected_count", sample_name)
  
  return(df)
}


# Merge all files into a single transcript expression matrix
expression_data <- Reduce(function(x, y) merge(x, y, by=c("gene_id", "transcript_id"), all=TRUE),
                          lapply(rna_files, read_expression_data))


# Save merged transcript expression matrix
fwrite(expression_data, expression_output)
cat("Transcript expression matrix saved to:", expression_output, "\n")

#################################################################
##  Step 3: Clean and Annotate Transcript Data
#################################################################

cat("Cleaning and annotating transcript expression data...\n")

# Load transcript expression data
expression_data <- fread(expression_output, sep = ",", header = TRUE)

# Convert all expression columns to numeric
expression_data[, (3:ncol(expression_data)) := lapply(.SD, as.numeric), .SDcols = 3:ncol(expression_data)]

# Remove rows where all expression values are zero
expression_data <- expression_data[rowSums(expression_data[, -c(1,2), with = FALSE]) > 0]

library(rtracklayer)

# Load Gencode GTF file
gencode <- import("gencode.v24.primary_assembly.annotation.gtf")
gencode_df <- as.data.frame(gencode)

# Keep only relevant columns for transcript-level annotation
gencode_df <- gencode_df[, c("transcript_id", "gene_id", "gene_name", "gene_type")]
# Remove exact duplicates from Gencode annotations
gencode_df <- unique(gencode_df)
gencode_df <- gencode_df[gencode_df$gene_type == "protein_coding", ]


# Convert to data.table if necessary
if (!is.data.table(expression_data)) setDT(expression_data)
if (!is.data.table(gencode_df)) setDT(gencode_df)

# Remove version suffix from transcript_id
expression_data[, transcript_id := sub("\\..*", "", transcript_id)]
gencode_df[, transcript_id := sub("\\..*", "", transcript_id)]

# Merge expression data with Gencode annotations using transcript_id
expression_data_annotated <- merge(expression_data, gencode_df, by = "transcript_id", all.x = TRUE)

# Remove rows where gene_name is missing (NA)
expression_data_annotated <- expression_data_annotated[!is.na(gene_name)]

# Save final annotated transcript matrix
annotated_output <- file.path(base_dir, "processed_data", "transcript_matrix_annotated.csv")
fwrite(expression_data_annotated, annotated_output)

cat("Annotated transcript expression matrix saved to:", annotated_output, "\n")

#################################################################
##  Step 4: Differential Expression Analysis Across All 4 Biosamples
#################################################################

cat("Running edgeR analysis for multi-sample comparison...\n")

library(edgeR)

# Load annotated transcript matrix
annotated_output <- file.path(base_dir, "processed_data", "transcript_matrix_annotated.csv")
expression_data_annotated <- fread(annotated_output, sep = ",", header = TRUE)

# Remove metadata columns (keep only transcript_id + expression counts)
expression_counts <- expression_data_annotated[, c("transcript_id", grep("^ENCFF", colnames(expression_data_annotated), value = TRUE)), with = FALSE]

# Convert to data frame first (to ensure row names are preserved)
expression_counts <- as.data.frame(expression_counts)

# Ensure transcript_id is set as rownames before converting
rownames(expression_counts) <- expression_counts$transcript_id

# Remove transcript_id column (no need in matrix)
expression_counts$transcript_id <- NULL

# Convert to matrix
expression_counts <- as.matrix(expression_counts)

# Verify row names
head(rownames(expression_counts))  # Should now return transcript IDs

# Define biological groups
group <- factor(metadata_clean$Biosample)  # Defines biosample category

# Create a design matrix without replicate_id (to avoid collinearity)
design <- model.matrix(~ 0 + group)  # No + replicate_id
colnames(design) <- gsub("group", "", colnames(design))  # Clean column names

y <- DGEList(counts = expression_counts, group = group)

# Normalization (TMM for library size differences)
y <- calcNormFactors(y)

# Estimate dispersions
y <- estimateDisp(y, design, robust=TRUE)

# Fit GLM with Quasi-Likelihood
fit <- glmQLFit(y, design, robust=TRUE)

# ANOVA-like test to find **any** differences across all 4 biosamples
qlf <- glmQLFTest(fit, coef = 2:ncol(design))  # Ensure intercept is not tested

# Extract results
res_edgeR <- topTags(qlf, n = Inf, adjust.method = "BH")

# Convert results to a data table
res_edgeR_dt <- as.data.table(res_edgeR$table)
res_edgeR_dt$transcript_id <- rownames(res_edgeR)

cat("Global test for differential expression across all biosamples complete. Results saved.\n")


#################################################################
##  Step 5: Store Differential Expression Results for Visualization
#################################################################

cat("Saving all differential expression results...\n")

# Extract results
res_edgeR <- topTags(qlf, n = Inf, adjust.method = "BH")

# Convert results to a data table
res_edgeR_dt <- as.data.table(res_edgeR$table)
res_edgeR_dt$transcript_id <- rownames(res_edgeR)

# Save all results (including non-significant genes)
fwrite(res_edgeR_dt, file.path(base_dir, "processed_data", "edgeR_all_genes.csv"))

# Adjust p-values for multiple testing (Benjamini-Hochberg)
res_edgeR_dt[, padj := p.adjust(PValue, method = "BH")]

# Filter significant genes (padj < 0.05)
significant_genes <- res_edgeR_dt[padj < 0.05]

# Save only significant genes
fwrite(significant_genes, file.path(base_dir, "results", "edgeR_significant_genes.csv"))

cat("All differential expression results saved for visualization.\n")
cat("Significant differentially expressed genes saved separately.\n")

#################################################################
##  Step 6: LogCPM Barplots for Specific Transcripts (Ordered by Biosample & Experiment)
#################################################################

cat("Creating logCPM barplots for selected transcripts...\n")

library(ggplot2)
library(data.table)

# 1) Compute logCPM from normalized counts
logCPM <- cpm(y, log=TRUE)  # Convert to log counts per million (logCPM)

# 2) Convert logCPM matrix to data.table for visualization
logCPM_df <- as.data.table(logCPM, keep.rownames = "transcript_id")  # Keep transcript_id for labeling

# 3) Reshape to long format
logCPM_long <- melt(
  logCPM_df, 
  id.vars = "transcript_id", 
  variable.name = "File_accession", 
  value.name = "logCPM"
)

# 4) Merge with metadata (to get Biosample, etc.)
logCPM_long <- merge(
  logCPM_long, 
  metadata_clean[, .(File_accession, Biosample, Experiment_accession)], 
  by = "File_accession",
  all.x = TRUE
)

# 5) Define plot colors and ordering
biosample_colors <- c(
  "A549"    = "#E76BF3", 
  "GM12878" = "#A3A500", 
  "HepG2"   = "#F8766D", 
  "K562"    = "#00BFC4"
)

logCPM_long[, Biosample := factor(Biosample, 
                                  levels = c("A549", "GM12878", "HepG2", "K562"))]
setorder(logCPM_long, Biosample, Experiment_accession)

# 6) Merge edgeR results with transcript annotations
#    We assume 'res_edgeR_dt' has at least: transcript_id, FDR, PValue, etc.
#    and 'expression_data_annotated' has transcript_id -> gene_name
res_edgeR_dt_annot <- merge(
  res_edgeR_dt, 
  unique(expression_data_annotated[, .(transcript_id, gene_name)]),
  by = "transcript_id",
  all.x = TRUE
)

# 7) List of genes of interest (epigenetic regulators, etc.)
genes_of_interest <- c(
  "DNMT1", "DNMT3A", "DNMT3B",
  "TET1", "TET2", "TET3",
  "EZH2", "SUZ12", "EED",
  "EHMT2", "SUV39H1", "SUV39H2",
  "KMT2A", "KMT2B",
  "PRMT1", "PRMT5",
  "KDM1A", "KDM4A", "KDM5B", "KDM6A",
  "MECP2", "MBD1", "MBD2", "MBD3", "MBD4",
  "EP300", "CREBBP", "KAT2B", "KAT2A",
  "HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC5",
  "HDAC6", "HDAC7", "HDAC8", "HDAC9", "HDAC10", "HDAC11",
  "SIRT1", "SIRT2", "SIRT3", "SIRT4", "SIRT5", "SIRT6", "SIRT7"
)

# 8) Identify matching transcripts in our annotated expression dataset
genes_in_data <- expression_data_annotated[gene_name %in% genes_of_interest]

# Collect the transcript IDs you want to plot
selected_transcripts <- unique(genes_in_data$transcript_id)

# Ensure the output directory exists
plot_dir <- file.path(base_dir, "plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# 9) Generate barplots for each transcript
for (tx in selected_transcripts) {
  
  # (a) Subset the expression data for this transcript
  transcript_data <- logCPM_long[transcript_id == tx]
  if (nrow(transcript_data) == 0) {
    next  # Skip if no expression data
  }
  
  # (b) Find the gene name
  gene_name <- unique(genes_in_data[transcript_id == tx]$gene_name)
  if (length(gene_name) == 0 || is.na(gene_name)) {
    gene_name <- "UnknownGene"
  }
  
  # (c) Look up FDR from the merged results
  fdr_row <- res_edgeR_dt_annot[transcript_id == tx]
  fdr_val <- NA
  if (nrow(fdr_row) > 0 && !is.na(fdr_row$FDR[1])) {
    fdr_val <- fdr_row$FDR[1]
  }
  
  # Format the FDR for display (e.g., scientific notation)
  fdr_str <- if (!is.na(fdr_val)) {
    formatC(fdr_val, format = "e", digits = 2)
  } else {
    "NA"
  }
  
  # (d) Order x-axis by Biosample and Experiment
  transcript_data[, File_accession := factor(
    File_accession,
    levels = unique(transcript_data$File_accession[order(
      Biosample, Experiment_accession
    )])
  )]
  
  # (e) Create the plot title: include gene name + transcript ID + FDR
  plot_title <- paste0(
    "Expression of ", gene_name, " (", tx, ")\n", 
    "FDR = ", fdr_str
  )
  
  # (f) Make the barplot
  p <- ggplot(transcript_data, aes(x = File_accession, y = logCPM, fill = Biosample)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = biosample_colors) +
    theme_minimal() +
    labs(
      title = plot_title,
      x = "Experiment (File Accession)",
      y = "log Counts Per Million (logCPM)"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # (g) Define the output filename (include gene name & transcript ID)
  plot_file_name <- paste0("logCPM_", gene_name, "_", tx, ".png")
  plot_path <- file.path(plot_dir, plot_file_name)
  
  # (h) Save the plot
  ggsave(plot_path, plot = p, width = 10, height = 6)
  
  cat("Saved:", plot_path, "\n")
}

cat("logCPM barplots for selected transcripts saved in 'plots/' folder.\n")

#################################################################
##  Step 7: PCA Analysis of RNA-seq Data
#################################################################

cat("Performing PCA analysis on transcript expression data...\n")

# Load required libraries
library(ggplot2)
library(data.table)
library(ggrepel)

# Load normalized logCPM values
logCPM_matrix <-  cpm(y, log=TRUE)  # Convert to log counts per million (logCPM)

# Convert to numeric matrix
logCPM_matrix <- as.matrix(logCPM_matrix)

# Perform PCA
pca <- prcomp(t(logCPM_matrix), scale. = TRUE)

# Extract PCA results
pca_data <- as.data.table(pca$x)
pca_data$Sample <- colnames(logCPM_matrix)  # Add sample names

# Merge with metadata to get biosample information
pca_data <- merge(pca_data, metadata_clean, by.x = "Sample", by.y = "File_accession")

# Define percentage variance explained by PC1 and PC2
pca_var <- summary(pca)$importance[2, ]
PC1_var <- round(pca_var[1] * 100, 2)
PC2_var <- round(pca_var[2] * 100, 2)

# PCA plot
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Biosample, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3) +
  theme_minimal() +
  labs(title = "PCA of RNA-Seq Data",
       x = paste0("PC1: ", PC1_var, "% variance"),
       y = paste0("PC2: ", PC2_var, "% variance")) +
  scale_color_manual(values = c("A549" = "#E76BF3", "GM12878" = "#A3A500", 
                                "HepG2" = "#F8766D", "K562" = "#00BFC4"))

# Save PCA plot
pca_plot_path <- file.path(base_dir, "plots", "PCA_plot.png")
ggsave(pca_plot_path, plot = p, width = 8, height = 6)

cat("PCA analysis complete. Plot saved at:", pca_plot_path, "\n")

#################################################################
## Perform t-SNE Analysis
#################################################################
library(Rtsne)  # Required for t-SNE

cat("Performing t-SNE analysis on transcript expression data...\n")

# Run t-SNE
set.seed(123)  # For reproducibility
tsne_results <- Rtsne(t(logCPM_matrix), perplexity = 4, check_duplicates = FALSE)

# Convert t-SNE output to data.table
tsne_data <- as.data.table(tsne_results$Y)
colnames(tsne_data) <- c("tSNE1", "tSNE2")
tsne_data$Sample <- colnames(logCPM_matrix)

# Merge with metadata
tsne_data <- merge(tsne_data, metadata_clean, by.x = "Sample", by.y = "File_accession")

# t-SNE plot
tsne_plot <- ggplot(tsne_data, aes(x = tSNE1, y = tSNE2, color = Biosample, label = Sample)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3) +
  theme_minimal() +
  labs(title = "t-SNE of RNA-Seq Data",
       x = "tSNE1",
       y = "tSNE2") +
  scale_color_manual(values = c("A549" = "#E76BF3", "GM12878" = "#A3A500", 
                                "HepG2" = "#F8766D", "K562" = "#00BFC4"))

# Save t-SNE plot
tsne_plot_path <- file.path(base_dir, "plots", "tSNE_plot.png")
ggsave(tsne_plot_path, plot = tsne_plot, width = 8, height = 6)

cat("t-SNE analysis complete. Plot saved at:", tsne_plot_path, "\n")

#################################################################
##  Step 8: Volcano Plots for Sample Pairs
#################################################################

cat("Generating Volcano Plots for all sample pairs...\n")

# Load DE results
deg_results <- fread(file.path(base_dir, "processed_data", "edgeR_all_genes.csv"))

# Define sample pairs for comparison
sample_pairs <- combn(unique(metadata_clean$Biosample), 2, simplify = FALSE)

for (pair in sample_pairs) {
  sample1 <- pair[1]
  sample2 <- pair[2]
  
  # Extract relevant comparison columns
  deg_subset <- deg_results[, .(transcript_id, logFC = get(paste0("logFC.", sample2)), PValue)]
  
  # Define significance thresholds
  logFC_threshold <- 1.0  # Adjust if needed
  pval_threshold <- 0.05
  
  # Create significance column
  deg_subset[, Significance := ifelse(PValue < pval_threshold & abs(logFC) > logFC_threshold,
                                      ifelse(logFC > 0, "Upregulated", "Downregulated"),
                                      "Not Significant")]
  
  # Volcano plot
  volcano_plot <- ggplot(deg_subset, aes(x = logFC, y = -log10(PValue), color = Significance)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
    theme_minimal() +
    labs(title = paste("Volcano Plot:", sample1, "vs", sample2),
         x = "Log2 Fold Change",
         y = "-log10 P-Value") +
    theme(legend.title = element_blank())
  
  # Save the plot
  volcano_plot_path <- file.path(base_dir, "plots", paste0("volcano_plot_", sample1, "_vs_", sample2, ".png"))
  ggsave(volcano_plot_path, plot = volcano_plot, width = 8, height = 6)
  
  cat("Saved Volcano Plot for:", sample1, "vs", sample2, "\n")
}

cat("Volcano Plots generation complete.\n")