#################################################################
##  
##  RNA-seq Expression Processing
##
##  Input:  
##   - Metadata file for sample labeling
##   - RNA-seq gene quantification files
##
##  Biosample Mapping:
##   - E123 = K562
##   - E118 = HepG2
##
##  Author: Daniel Batyrev 777634015
#################################################################

# Clear R working environment
rm(list = ls())

# Set working directory
this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this.dir)

# Function to detach all non-base packages
detachAllPackages <- function() {
  basic.packages <-
    c("package:stats", "package:graphics", "package:grDevices",
      "package:utils", "package:datasets", "package:methods", "package:base")
  
  package.list <- search()[grepl("package:", search())]
  package.list <- setdiff(package.list, basic.packages)
  
  if (length(package.list) > 0)
    for (package in package.list) detach(package, character.only = TRUE)
}
detachAllPackages()

################################################ LIBRARIES #################################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

################################################ INPUT ###################################################

# Define paths
base_folder <- "D:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/expression"
metadata_file <- file.path(base_folder, "metadata.tsv")
data_folder <- file.path(base_folder, "rna-seq")

# Suggested colors for the biosamples
group.colors <- c(
  HepG2    = "#F8766D",   # Warm reddish color
  K562     = "#00BFC4",   # Cool cyan color
  GM12878  = "#A3A500",   # Yellow-green color
  A549     = "#E76BF3"    # Vibrant purple
)

################################################ LOAD METADATA ###################################################
metadata <- fread(metadata_file, sep = "\t", header = TRUE, quote = "")

# Extract relevant columns for file ID to biosample mapping
metadata_clean <- metadata %>%
  select(`File accession`, `Biosample term name`) %>%
  rename(File_ID = `File accession`, Biosample = `Biosample term name`)

################################################ LOAD EXPRESSION FILES ###################################################

# List all `.tsv` expression files in `rna-seq/`
file_list <- list.files(data_folder, pattern = "*.tsv", full.names = TRUE)

# Function to read and extract relevant columns
read_expression_data <- function(file) {
  df <- fread(file, sep = "\t")
  
  # Select only gene_id and TPM (Transcripts Per Million)
  df <- df %>% select(gene_id, TPM)
  
  # Extract file name (ENCFF ID) from full path
  file_id <- tools::file_path_sans_ext(basename(file))
  
  # Match file ID with biosample name from metadata
  sample_name <- metadata_clean %>%
    filter(File_ID == file_id) %>%
    pull(Biosample)
  
  # If biosample name is found, rename column
  if (length(sample_name) > 0) {
    colnames(df)[2] <- sample_name
  } else {
    colnames(df)[2] <- file_id  # Fallback: use file ID if no match
  }
  
  return(df)
}

# Read and merge all expression files
expression_data <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE),
                          lapply(file_list, read_expression_data))

# Save summary table
summary_file <- file.path(base_folder, "expression_summary.csv")
fwrite(expression_data, summary_file)

################################################ PLOT GENE EXPRESSION ###################################################

# Create output directory for plots
output_folder <- file.path(base_folder, "plots")
dir.create(output_folder, showWarnings = FALSE)

# Function to generate plots for each gene
plot_gene_expression <- function(gene) {
  gene_data <- expression_data %>%
    filter(gene_id == gene) %>%
    pivot_longer(-gene_id, names_to = "Sample", values_to = "TPM")
  
  p <- ggplot(gene_data, aes(x = Sample, y = TPM, fill = Sample)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = group.colors) +
    theme_minimal() +
    labs(title = paste("Expression of Gene", gene), x = "Sample", y = "TPM") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(filename = file.path(output_folder, paste0(gene, ".png")),
         plot = p, width = 8, height = 6, limitsize = FALSE)
}

# Generate plots for all genes
unique_genes <- unique(expression_data$gene_id)
lapply(unique_genes, plot_gene_expression)

cat("Processing complete!\n- Summary saved as 'expression_summary.csv'\n- Plots saved in 'plots' folder.\n")
