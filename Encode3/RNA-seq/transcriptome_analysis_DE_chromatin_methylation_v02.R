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
base_dir <- getwd()
#"D:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/RNA-seq"
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
##  Step 4a: Summarize technical replicates by Experiment_accession
##           and run global ANOVA-like DE test
#################################################################

cat("Merging technical replicates by Experiment_accession, then running edgeR.\n")

library(data.table)
library(edgeR)

# -------------------------------
# (A) Load Data & Metadata
# -------------------------------

# 1) Annotated expression matrix (or counts), e.g. 'expression_counts'
#    - rows = genes/transcripts
#    - columns = File_accession
#    - we'll assume you already loaded it into 'expression_counts' object
#    - also assume 'metadata_clean' is available

# Example (if you were reading from CSV):
# expression_counts <- fread(file.path(base_dir, "processed_data", "transcript_matrix_annotated.csv"))

# 2) Ensure your 'expression_counts' is a matrix/data.table with rownames = gene_id
#    and colnames = File_accession
#    If it's not, adapt accordingly.

# -------------------------------
# (B) Melt, Merge, Summarize
# -------------------------------

# Convert to data.table
dt_counts <- as.data.table(expression_data_annotated)

# Keep only relevant columns (File_accession from metadata_clean)
relevant_columns <- c("transcript_id", metadata_clean$File_accession)
dt_counts <- dt_counts[, ..relevant_columns]

# Set transcript_id as row names
setkey(dt_counts, transcript_id)

# Print confirmation
cat("Filtered expression matrix saved in dt_counts\n")
print(head(dt_counts))


# Melt to long format: one row per (gene_id, File_accession)
long_counts <- melt(
  dt_counts,
  id.vars = "transcript_id",
  variable.name = "File_accession",
  value.name = "count"
)

# Merge with metadata to find each row's Experiment_accession
long_counts <- merge(
  long_counts,
  metadata_clean[, .(File_accession, Experiment_accession, Biosample)],
  by = "File_accession"
)

# Sum counts by (gene_id, Experiment_accession)
summed_counts <- long_counts[
  , .(count = sum(count)), 
  by = .(transcript_id, Experiment_accession)
]

# Cast back to wide format: one column per Experiment_accession
merged_counts <- dcast(
  summed_counts,
  transcript_id ~ Experiment_accession,
  value.var = "count",
  fill = 0
)

# Convert to matrix (row.names=gene_id)
transcript_id <- merged_counts$transcript_id
merged_counts$transcript_id <- NULL
expression_merged <- as.matrix(merged_counts)
rownames(expression_merged) <- transcript_id

# Create a condensed metadata for each Experiment_accession
exp_meta <- unique(metadata_clean[, .(Experiment_accession, Biosample)])
# Ensure the order matches columns of 'expression_merged'
setkey(exp_meta, Experiment_accession)
col_order <- data.table(Experiment_accession = colnames(expression_merged))
setkey(col_order, Experiment_accession)
exp_meta <- exp_meta[col_order]

cat("Merged technical replicates -> now have", ncol(expression_merged), "columns.\n")
cat("Unique Experiment_accession per biosample:\n")
print(table(exp_meta$Biosample))

# -------------------------------
# (C) edgeR Setup
# -------------------------------

# 1) Define group factor from the condensed metadata
group <- factor(exp_meta$Biosample) 
cat("Levels of 'group':\n")
print(levels(group))

# 2) Build an intercept-based design
design <- model.matrix(~ group)
colnames(design)
cat("Design matrix:\n")
print(design)

# 3) Create DGEList and normalize
y <- DGEList(counts = expression_merged, group = group)
y <- calcNormFactors(y)

# 4) Estimate dispersion and fit the model
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)

# 5) Omnibus (ANOVA-like) test: test all non-intercept coefficients
#    i.e., do any biosamples differ from the baseline?
qlf <- glmQLFTest(fit, coef = 2:ncol(design))

# 6) Extract and save results
res_edgeR <- topTags(qlf, n = Inf, adjust.method = "BH")
res_edgeR_dt <- as.data.table(res_edgeR$table)
res_edgeR_dt[, transcript_id := rownames(res_edgeR$table)]

# Write to file
out_global <- file.path(base_dir, "processed_data", "edgeR_omnibus_all_genes.csv")
fwrite(res_edgeR_dt, out_global)

cat("Global ANOVA-like test results saved to:", out_global, "\n")

# 7) Optional: filter significant genes
res_edgeR_dt[, padj := p.adjust(PValue, method = "BH")]
sig_genes <- res_edgeR_dt[padj < 0.05]
cat("Number of significant transcripts (FDR < 0.05):", nrow(sig_genes), "\n")

sig_file <- file.path(base_dir, "results", "edgeR_omnibus_significant.csv")
fwrite(sig_genes, sig_file)
cat("Significant transcripts saved to:", sig_file, "\n")

cat("Step 4a complete: technical replicates merged, omnibus test done.\n\n")


#################################################################
##  Step 4b: Pairwise Contrasts for Biosamples + Volcano Plots
#################################################################

cat("Performing pairwise comparisons among biosamples with merged replicates...\n")

library(ggplot2)

# (A) Define thresholds
fdr_threshold <- 0.05
logfc_threshold <- 1.0  # e.g., 2-fold

# (B) Identify unique biosamples from 'group'
all_groups <- levels(group)

# (C) Generate all pairwise combos
all_pairs <- combn(all_groups, 2, simplify = FALSE)

# (D) Create output dirs
pairwise_dir <- file.path(base_dir, "results", "pairwise_merged")
volcano_dir  <- file.path(base_dir, "plots", "volcano_pairwise_merged")
if (!dir.exists(pairwise_dir)) dir.create(pairwise_dir, recursive = TRUE)
if (!dir.exists(volcano_dir)) dir.create(volcano_dir, recursive = TRUE)

# (E) Helper to build contrast vector for group1 - group2
makeContrastVector <- function(g1, g2, design_colnames) {
  cv <- rep(0, length(design_colnames))
  names(cv) <- design_colnames
  
  # e.g., colnames might be: (Intercept), groupGM12878, groupHepG2, groupK562
  # If A549 is baseline, "groupA549" won't appear because it's intercept
  g1_col <- paste0("group", g1)
  g2_col <- paste0("group", g2)
  
  has_g1 <- g1_col %in% design_colnames
  has_g2 <- g2_col %in% design_colnames
  
  if (has_g1 && !has_g2) {
    # g2 is baseline => test g1 vs. baseline
    cv[g1_col] <- 1
  } else if (!has_g1 && has_g2) {
    # g1 is baseline => test g2 vs. baseline
    # but this might invert the direction, depending on what you want 
    # e.g., if you want "g1 - g2", you'd put -1 in g2_col
    # We'll keep it symmetrical for example
    cv[g2_col] <- 1
  } else if (has_g1 && has_g2) {
    # both non-baseline => group1 - group2
    cv[g1_col] <- 1
    cv[g2_col] <- -1
  } else {
    warning("No recognized columns for groups: ", g1, " vs ", g2)
  }
  return(cv)
}

# (F) Loop over pairs
for (p in all_pairs) {
  g1 <- p[1]
  g2 <- p[2]
  
  contrast_label <- paste0(g1, "_vs_", g2)
  
  # Build contrast vector
  contrast_vec <- makeContrastVector(g1, g2, colnames(design))
  
  # Run test
  qlf_contrast <- glmQLFTest(fit, contrast = contrast_vec)
  
  # Extract results
  top_res <- topTags(qlf_contrast, n = Inf, adjust.method = "BH")
  pairwise_dt <- as.data.table(top_res$table)
  pairwise_dt[, transcript_id := rownames(top_res$table)]
  
  # Save table
  out_file <- file.path(pairwise_dir, paste0("edgeR_pairwise_", contrast_label, ".csv"))
  fwrite(pairwise_dt, out_file)
  
  # Classify significance
  pairwise_dt[, abs_logFC := abs(logFC)]
  pairwise_dt[, Significance := "Not Significant"]
  pairwise_dt[FDR < fdr_threshold & logFC >  logfc_threshold, Significance := "Upregulated"]
  pairwise_dt[FDR < fdr_threshold & logFC < -logfc_threshold, Significance := "Downregulated"]
  
  # Count
  num_sig <- nrow(pairwise_dt[Significance != "Not Significant"])
  cat("Pairwise comparison:", contrast_label, 
      "- DE Genes (FDR<", fdr_threshold, "& |logFC|>", logfc_threshold, "):", 
      num_sig, "/", nrow(pairwise_dt), "\n")
  
  # Volcano plot
  volcano_plot <- ggplot(pairwise_dt, aes(x = logFC, y = -log10(PValue), color = Significance)) +
    geom_point(alpha = 0.8, size = 1.5) +
    theme_minimal() +
    scale_color_manual(values = c("Upregulated" = "red",
                                  "Downregulated" = "blue",
                                  "Not Significant" = "gray70")) +
    labs(
      title = paste("Volcano Plot:", contrast_label),
      subtitle = paste(num_sig, "DE genes at FDR<", fdr_threshold),
      x = "Log2 Fold Change",
      y = "-log10(PValue)"
    ) +
    theme(legend.title = element_blank())
  
  plot_file <- file.path(volcano_dir, paste0("volcano_", contrast_label, ".png"))
  ggsave(plot_file, volcano_plot, width = 8, height = 6)
  cat("Saved volcano plot to:", plot_file, "\n\n")
}

cat("Step 4b complete: Pairwise contrasts and Volcano plots done.\n")

#################################################################
##  Step 5: Store Differential Expression Results for Visualization
#################################################################

cat("Storing ANOVA-like differential expression results...\n")

# 'qlf' was created in Step 4a after the omnibus test with merged replicates
# For clarity, we extract results again (though you may already have done so):
res_edgeR <- topTags(qlf, n = Inf, adjust.method = "BH")

# Convert to data.table
res_edgeR_dt <- as.data.table(res_edgeR$table)
res_edgeR_dt[, transcript_id := rownames(res_edgeR$table)]

# Save all results (including non-significant)
out_all <- file.path(base_dir, "processed_data", "edgeR_all_genes.csv")
fwrite(res_edgeR_dt, out_all)

cat("All DE results saved to:", out_all, "\n")

# Add an FDR column if not already present
if (!"FDR" %in% colnames(res_edgeR_dt)) {
  # 'topTags' typically has a column named FDR, but in case it doesn't, you can do:
  res_edgeR_dt[, FDR := p.adjust(PValue, method = "BH")]
}

# Filter significant genes (FDR < 0.05)
sig_genes <- res_edgeR_dt[FDR < 0.05]
cat("Number of significant transcripts (FDR < 0.05):", nrow(sig_genes), "\n")

# Save only significant genes
out_sig <- file.path(base_dir, "results", "edgeR_significant_genes.csv")
fwrite(sig_genes, out_sig)

cat("Significant DE genes saved to:", out_sig, "\n")
cat("Step 5 complete.\n\n")


#################################################################
##  Step 6: LogCPM Barplots for Specific Transcripts
#################################################################

cat("Creating logCPM barplots for selected transcripts...\n")

library(ggplot2)
library(data.table)

# 1) Compute logCPM from merged, normalized counts
#    'y' is the DGEList with merged replicates
logCPM <- cpm(y, log = TRUE)  # one column per Experiment_accession

# Convert to data.table for visualization
logCPM_dt <- as.data.table(logCPM, keep.rownames = "transcript_id")
# 'transcript_id' is now the first column, other columns are the experiment IDs

# 2) Melt to long format
logCPM_long <- melt(
  logCPM_dt,
  id.vars = "transcript_id",
  variable.name = "Experiment_accession",
  value.name = "logCPM"
)

# 3) Merge with the condensed metadata (exp_meta)
#    exp_meta has columns: Experiment_accession, Biosample
#    So each row is one biological replicate
logCPM_long <- merge(
  logCPM_long,
  exp_meta, 
  by = "Experiment_accession"
)

# 4) Order the factor levels if desired
biosample_colors <- c(
  "A549"    = "#E76BF3",
  "GM12878" = "#A3A500",
  "HepG2"   = "#F8766D",
  "K562"    = "#00BFC4"
)

logCPM_long[, Biosample := factor(Biosample, 
                                  levels = c("A549", "GM12878", "HepG2", "K562"))]

# 5) Merge edgeR results with transcript annotations
#    Suppose 'res_edgeR_dt' has transcript_id, FDR, etc.
#    Suppose you have 'expression_data_annotated' or similar, 
#    which has (transcript_id, gene_name).
res_edgeR_dt_annot <- merge(
  res_edgeR_dt,
  unique(expression_data_annotated[, .(transcript_id, gene_name)]),
  by = "transcript_id",
  all.x = TRUE
)

# 6) Choose transcripts of interest
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

# Identify matching transcripts
genes_in_data <- expression_data_annotated[gene_name %in% genes_of_interest]
selected_transcripts <- unique(genes_in_data$transcript_id)

# 7) Output directory
plot_dir <- file.path(base_dir, "plots", "barplots_merged")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# 8) Generate barplots
for (tx in selected_transcripts) {
  
  # Subset expression
  sub_dt <- logCPM_long[transcript_id == tx]
  if (nrow(sub_dt) == 0) next
  
  # Find gene name
  gene_name <- unique(genes_in_data[transcript_id == tx]$gene_name)
  if (length(gene_name) == 0 || is.na(gene_name)) gene_name <- "UnknownGene"
  
  # Look up FDR from DE results
  fdr_val <- res_edgeR_dt_annot[transcript_id == tx, FDR][1]  # might be NA if not in table
  # Conditional formatting: scientific notation for small values, regular for larger ones
  fdr_str <- if (!is.na(fdr_val)) {
    if (fdr_val < 0.01) {
      formatC(fdr_val, format = "e", digits = 2)  # Scientific notation for very small values
    } else {
      formatC(fdr_val, format = "f", digits = 3)  # Regular float format
    }
  } else "NA"
  
  # Create the plot title
  plot_title <- paste0("Expression of ", gene_name, " (", tx, ")\nFDR = ", fdr_str)
  
  # Reorder x-axis by Biosample
  sub_dt[, Experiment_accession := factor(Experiment_accession,
                                          levels = unique(sub_dt$Experiment_accession[order(Biosample)]))]
  
  p <- ggplot(sub_dt, aes(x = Experiment_accession, y = logCPM, fill = Biosample)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = biosample_colors) +
    theme_minimal() +
    labs(
      title = plot_title,
      x = "Experiment Accession",
      y = "log Counts Per Million (logCPM)"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save
  plot_file <- file.path(plot_dir, paste0("logCPM_", gene_name, "_", tx, ".png"))
  ggsave(plot_file, plot = p, width = 10, height = 6)
  
  cat("Saved barplot for", tx, "->", plot_file, "\n")
}

cat("Barplots for selected transcripts saved in:", plot_dir, "\n")


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
pca_data$Experiment_accession<- colnames(logCPM_matrix)  # Add sample names

# Merge with metadata to get biosample information
pca_data <- merge(pca_data, metadata_clean, by.x = "Experiment_accession", by.y = "Experiment_accession")

# Define percentage variance explained by PC1 and PC2
pca_var <- summary(pca)$importance[2, ]
PC1_var <- round(pca_var[1] * 100, 2)
PC2_var <- round(pca_var[2] * 100, 2)

# PCA plot
p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Biosample, label = Experiment_accession)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3) +
  theme_minimal() +
  labs(title = "PCA of RNA-Seq Data",
       x = paste0("PC1: ", PC1_var, "% variance"),
       y = paste0("PC2: ", PC2_var, "% variance")) +
  scale_color_manual(values = c("A549" = "#E76BF3", "GM12878" = "#A3A500", 
                                "HepG2" = "#F8766D", "K562" = "#00BFC4"))

# Save PCA plot
pca_plot_path <- file.path(base_dir, "plots", "PCA_plot_Experiment_accession.png")
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
tsne_data$Experiment_accession <- colnames(logCPM_matrix)

# Merge with metadata
tsne_data <- merge(tsne_data, metadata_clean, by.x = "Experiment_accession", by.y = "Experiment_accession")

# t-SNE plot
tsne_plot <- ggplot(tsne_data, aes(x = tSNE1, y = tSNE2, color = Biosample, label = Experiment_accession)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3) +
  theme_minimal() +
  labs(title = "t-SNE of RNA-Seq Data",
       x = "tSNE1",
       y = "tSNE2") +
  scale_color_manual(values = c("A549" = "#E76BF3", "GM12878" = "#A3A500", 
                                "HepG2" = "#F8766D", "K562" = "#00BFC4"))

# Save t-SNE plot
tsne_plot_path <- file.path(base_dir, "plots", "tSNE_plot_Experiment_accession.png")
ggsave(tsne_plot_path, plot = tsne_plot, width = 8, height = 6)

cat("t-SNE analysis complete. Plot saved at:", tsne_plot_path, "\n")
