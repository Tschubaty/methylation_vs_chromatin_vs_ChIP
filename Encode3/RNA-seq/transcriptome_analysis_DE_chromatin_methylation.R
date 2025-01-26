#################################################################
##  Script: RNA_transcript_analysis.R
##  Step 1: Load and merge ENCODE transcript quantifications
##
##  Author: Daniel Batyrev
#################################################################

# Load required libraries
library(data.table)
library(DESeq2)
library(ggplot2)

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

# List all transcript quantification files
rna_files <- list.files(file.path(base_dir, "raw_data2"), pattern = "*.tsv", full.names = TRUE)

# Function to extract "expected_count" for transcripts
read_expression_data <- function(file) {
  df <- fread(file, sep = "\t", header = TRUE)[, .(gene_id, transcript_id = `transcript_id(s)`, expected_count)]
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
