#################################################################
##  ENCODE Metadata Input Creator
##
##  Input:  meta_summery.RDS from Encode3/ChIP-seq/2.filter_samples
##  Output: encode_metadata_input.csv
##  Purpose: Prepare correct metadata input for ENCODE antibody scraper
##
##  Author: Daniel Batyrev
##  Date: 2026
#################################################################

rm(list = ls())

cluster <- FALSE

if (cluster) {
  base_dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/Encode3"
  this.dir <- file.path(base_dir, "Antibody_Info")
} else {
  this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  base_dir <- dirname(this.dir)
}

setwd(this.dir)

library(dplyr)
library(readr)

cat("===================================================\n")
cat("ENCODE Metadata Input Creator\n")
cat("===================================================\n\n")

# Correct metadata source
meta_file <- file.path(
  base_dir,
  "ChIP-seq",
  "2.filter_samples",
  "meta_summery.RDS"
)

if (!file.exists(meta_file)) {
  stop("meta_summery.RDS not found at: ", meta_file)
}

cat("Loading metadata from:\n")
cat(meta_file, "\n\n")

meta_summery <- readRDS(meta_file)

cat("Rows loaded:", nrow(meta_summery), "\n")
cat("Unique files:", n_distinct(meta_summery$File.accession), "\n")
cat("Unique experiments:", n_distinct(meta_summery$Experiment.accession), "\n\n")

# Convert to Python scraper input format
encode_metadata_input <- meta_summery %>%
  transmute(
    File_Accession = File.accession,
    Experiment_Accession = Experiment.accession,
    Target = Experiment.target,
    Biosample = Biosample.term.name,
    Biological_Replicate = Biological.replicate.s.,
    Technical_Replicate = Technical.replicate.s.,
       Experiment_Date_Released = Experiment.date.released,
  ) %>%
  distinct(File_Accession, .keep_all = TRUE)

cat("After deduplication by File_Accession:\n")
cat("Rows:", nrow(encode_metadata_input), "\n")
cat("Unique files:", n_distinct(encode_metadata_input$File_Accession), "\n")
cat("Unique experiments:", n_distinct(encode_metadata_input$Experiment_Accession), "\n\n")

# Save CSV for Python scraper
output_file <- file.path(this.dir, "encode_metadata_input.csv")

write_csv(encode_metadata_input, output_file)

cat("Saved Python input CSV to:\n")
cat(output_file, "\n\n")

cat("First 10 rows:\n")
print(head(encode_metadata_input, 10))

cat("\n===================================================\n")
cat("COMPLETE! Ready for Python scraper.\n")
cat("===================================================\n")