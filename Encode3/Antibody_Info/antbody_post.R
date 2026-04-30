#################################################################
##  ENCODE Metadata Aggregator
##
##  Input:  final_summary.rds from analysis_chromatin_state/
##  Output: antibody/QC comparison tables and plots
##  Purpose: merge chromatin-state QC summaries with ENCODE antibody metadata
##
##  Author: Daniel Batyrev
##  Date: 2026
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

# Load libraries
library(dplyr)
library(data.table)
library(readr)
library(tidyr)
library(ggplot2)



# Define project folders
setwd(this.dir)
base_dir <- dirname(this.dir)
analysis_dir <- file.path(base_dir, "analysis_chromatin_state")
antibody_dir <- this.dir
output_dir_plots <- file.path(antibody_dir, "plots")

dir.create(output_dir_plots, recursive = TRUE, showWarnings = FALSE)


# Load final chromatin summary
final_summary_file <- file.path(analysis_dir, "final_summary.rds")

if (!file.exists(final_summary_file)) {
  stop("final_summary.rds not found at: ", final_summary_file)
}


final_summary <- readRDS(final_summary_file)

cat("Loaded final_summary from:\n", final_summary_file, "\n")
cat("Rows:", nrow(final_summary), "\n")
cat("Columns:", paste(colnames(final_summary), collapse = ", "), "\n")


final_summary <- final_summary %>%
  tidyr::separate(
    protein_hit,
    into = c("Protein", "Experiment_Accession", "Motif"),
    sep = "_",
    remove = FALSE,
    extra = "merge"
  )

# load antibody mapping (output of your python script)
antibody_df <- read_csv("ENCODE_antibody_complete_mapping.csv")


final_summary_merged <- final_summary %>%
  left_join(
    antibody_df,
    by = c(
      "Experiment_Accession" = "File_Accession",
      "biosample" = "Biosample"
    )
  )

targets_multi_antibody <- antibody_df %>%
  filter(!is.na(Antibody_Accession), Antibody_Accession != "N/A") %>%
  distinct(Target, Experiment_Accession, File_Accession, Biosample, Antibody_Accession, Antibody_Title) %>%
  group_by(Target) %>%
  filter(n_distinct(Antibody_Accession) > 1) %>%
  arrange(Target, Antibody_Accession, Experiment_Accession, File_Accession) %>%
  ungroup()

n_targets <- length(unique(targets_multi_antibody$Target))
n_antibodies <- length(unique(targets_multi_antibody$Antibody_Accession))

cat("\n===== Antibody Diversity Summary =====\n")
cat(sprintf("Targets with >1 antibody : %d\n", n_targets))
cat(sprintf("Unique antibodies        : %d\n", n_antibodies))
cat("=====================================\n\n")