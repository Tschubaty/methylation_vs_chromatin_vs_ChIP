#################################################################
##  ENCODE Metadata Aggregator
##
##  Input: metadata.tsv files from all ChIP-seq sample directories
##  Output: encode_metadata_input.csv (combined and cleaned)
##  Purpose: Prepare data for ENCODE antibody scraper
##
##  Author: Your Name
##  Date: 2024
#################################################################

# Clear environment
#Clear R working environment
rm(list = ls())
cluster <- FALSE
if (cluster) {
  this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/Encode3/"
  picuture_file_extension <- "pdf"
} else{
  this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  picuture_file_extension <- "png"
}

setwd(this.dir)

# Set working directory to Encode3/ChIP-seq/1.download_from_Encode3
# Adjust this path to your actual location
input_dir <- "C:/Users/Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/ChIP-seq/1.download_from_Encode3"

# Load libraries
library(dplyr)
library(data.table)
library(readr)

cat("===================================================\n")
cat("ENCODE Metadata Aggregator\n")
cat("===================================================\n\n")

# Define sample directories to process
sample_dirs <- c(
  "A549july2023",
  "K562july2023", 
  "HepG2july2023", 
  "GM12878july2023"
)

# Initialize list to store all metadata
all_metadata <- list()

# Loop through each directory
for (sample_dir in sample_dirs) {
  metadata_file <- file.path(input_dir, sample_dir, "metadata.tsv")
  print(metadata_file)
  cat(sprintf("Processing: %s\n", sample_dir))
  
  # Check if file exists
  if (!file.exists(metadata_file)) {
    cat(sprintf("  WARNING: File not found - %s\n", metadata_file))
    next
  }
  
  tryCatch({
    # Read metadata file
    cat("  Reading metadata.tsv...")
    meta_raw <- read.table(
      metadata_file, 
      header = TRUE, 
      sep = "\t", 
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = ""
    )
    cat(sprintf(" OK (%d rows)\n", nrow(meta_raw)))
    
    # Select and rename relevant columns
    # Note: Column names in ENCODE metadata use dots instead of spaces
    meta_filtered <- meta_raw %>%
      filter(File.Status == "released") %>%
      select(
        File.accession,
        Experiment.accession,
        Experiment.target,
        Biosample.term.name,
        Biological.replicate.s.,
        Technical.replicate.s.,
        Lab,
        Experiment.date.released,
        File.analysis.status
      ) %>%
      mutate(
        Sample = sample_dir,
        Sample_Name = sub("july2023", "", sample_dir)
      )
    
    # Store in list
    all_metadata[[sample_dir]] <- meta_filtered
    
    cat(sprintf("  Extracted %d released files\n", nrow(meta_filtered)))
    
  }, error = function(e) {
    cat(sprintf("  ERROR reading file: %s\n", e$message))
  })
}

# Combine all metadata into single dataframe
cat("\n---------------------------------------------------\n")
cat("Combining all samples...\n")

combined_metadata <- bind_rows(all_metadata)

cat(sprintf("Total records combined: %d\n", nrow(combined_metadata)))
cat(sprintf("Unique experiments: %d\n", length(unique(combined_metadata$Experiment.accession))))
cat(sprintf("Unique files: %d\n", length(unique(combined_metadata$File.accession))))
cat(sprintf("Samples: %s\n", paste(unique(combined_metadata$Sample_Name), collapse = ", ")))

# Remove duplicates (keep first occurrence of each experiment accession)
combined_metadata_dedup <- combined_metadata %>%
  distinct(Experiment.accession, .keep_all = TRUE)

cat(sprintf("After removing experiment duplicates: %d records\n", nrow(combined_metadata_dedup)))

# Rename columns for clarity in output
combined_metadata_final <- combined_metadata_dedup %>%
  rename(
    File_Accession = File.accession,
    Experiment_Accession = Experiment.accession,
    Target = Experiment.target,
    Biosample = Biosample.term.name,
    Biological_Replicate = Biological.replicate.s.,
    Technical_Replicate = Technical.replicate.s.,
    Lab = Lab,
    Experiment_Date_Released = Experiment.date.released,
    File_Analysis_Status = File.analysis.status,
    Sample_Directory = Sample,
    Sample_Name = Sample_Name
  ) %>%
  select(
    File_Accession,
    Experiment_Accession,
    Target,
    Biosample,
    Biological_Replicate,
    Technical_Replicate,
    Sample_Name,
    Sample_Directory,
    Lab,
    Experiment_Date_Released,
    File_Analysis_Status
  )

# Display summary statistics
cat("\n---------------------------------------------------\n")
cat("Summary Statistics:\n")
cat("---------------------------------------------------\n")

cat("\nTargets found:\n")
targets_summary <- combined_metadata_final %>%
  group_by(Target) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  arrange(desc(Count))
print(targets_summary)

cat("\nSamples breakdown:\n")
sample_summary <- combined_metadata_final %>%
  group_by(Sample_Name) %>%
  summarise(
    Files = n(), 
    Experiments = n_distinct(Experiment_Accession), 
    .groups = 'drop'
  )
print(sample_summary)

# Save to CSV
output_file <- file.path(this.dir, "encode_metadata_input.csv")
cat(sprintf("\n---------------------------------------------------\n"))
cat(sprintf("Saving to: %s\n", output_file))

write.csv(
  combined_metadata_final,
  output_file,
  row.names = FALSE,
  quote = TRUE
)

cat(sprintf("Successfully saved %d records to encode_metadata_input.csv\n", nrow(combined_metadata_final)))

# Display first few rows
cat("\nFirst 10 rows of output:\n")
print(head(combined_metadata_final, 10))

cat("\n===================================================\n")
cat("COMPLETE! Ready for Python scraper.\n")
cat("===================================================\n")
