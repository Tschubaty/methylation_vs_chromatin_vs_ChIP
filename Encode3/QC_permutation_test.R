#################################################################
##  Chip-seq vs Methylation vs Chromatin Analysis
##
##  Description: Processes ChIP-seq, methylation, and chromatin data
##               to generate summary statistics and ratios.
##
##  Input: Methylation and protein peak data (HG38)
##         Located in Encode3/meme/fimo_single_experiments_on_known_motif_and_peaks_add_methylation_v3/
##
##  Output: Summary statistics saved in:
##          Encode3/permutation_test/
##
##  Version: 24.12.2024
##  Author: Daniel Batyrev (777634015)
#################################################################

# Clear R working environment
rm(list = ls())
cluster <- FALSE

# Set working directory
if (cluster) {
  this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/Encode3/"
  picuture_file_extension <- "pdf"
} else {
  this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  picuture_file_extension <- "svg"
}
setwd(this.dir)

# Function to detach all loaded packages except base ones
detachAllPackages <- function() {
  basic.packages <- c(
    "package:stats",
    "package:graphics",
    "package:grDevices",
    "package:utils",
    "package:datasets",
    "package:methods",
    "package:base"
  )
  package.list <- search()[ifelse(unlist(gregexpr("package:", search())) == 1, TRUE, FALSE)]
  package.list <- setdiff(package.list, basic.packages)
  if (length(package.list) > 0) {
    for (package in package.list)
      detach(package, character.only = TRUE)
  }
}

detachAllPackages()

#################################### Libraries ###################################
# Load necessary libraries
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rlang)

#################################### Constants ###################################
# Start time for the script
start_script <- Sys.time()

# Define input and output directories
input_folder <- file.path(
  this.dir,
  "permutation_test_2"
)
output_folder <- file.path(this.dir, "QC_permutation_test")

# Create output directory if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

group.colors <- c(
  HepG2    = "#F8766D",
  # A warm reddish color
  K562     = "#00BFC4",
  # A cool cyan color
  GM12878  = "#A3A500",
  # A yellow-green color, good contrast with the others
  A549     = "#E76BF3"   # A vibrant purple color
)

# Define chromatin state colors
chromatin_state_colors_short <- c(
  "1" = "#FF0000",
  "2" = "#FF4500",
  "3" = "#FF9900",
  "4" = "#FFCC00",
  "5" = "#00CC00",
  "6" = "#006400",
  "7" = "#FFD700",
  "8" = "#FFD700",
  "9" = "#FFFF00",
  "10" = "#FFDD00",
  "11" = "#FFEA73",
  "12" = "#9370DB",
  "13" = "#C0C0C0",
  "14" = "#FF4500",
  "15" = "#FFDD00",
  "16" = "#808080",
  "17" = "#A9A9A9",
  "18" = "#000000"
)

input_rds <- "stratified_test_100k_26.10.2025.rds"

stratified_test <- readRDS(file.path(input_folder,input_rds))

str(stratified_test)
# ============================================================
# QC 1: How many proteins have enough CpGs in both groups
#       in at least 1, 2, 3, 4, or 5 chromatin states?
# ============================================================

min_cpg <- 50

# Total number of unique proteins in the full result table,
# independent of whether they pass the CpG threshold.
total_proteins <- n_distinct(stratified_test$protein)

# First, mark each row/test as valid only if both groups have enough CpGs.
# Then collapse to protein + chromatin state level.
# This avoids counting the same state multiple times if one protein has
# several motifs or biosample comparisons in the same chromatin state.
protein_valid_states <- stratified_test %>%
  mutate(
    chromatin_state = str_remove(data_set, "state_"),
    enough_both = !is.na(n_1) & !is.na(n_2) & n_1 >= min_cpg & n_2 >= min_cpg
  ) %>%
  filter(enough_both) %>%
  distinct(protein, chromatin_state) %>%
  count(protein, name = "n_valid_states")

# Count how many proteins have enough CpGs in both groups
# in at least 1, 2, 3, 4, or 5 chromatin states.
qc_proteins_by_valid_states <- tibble(
  min_valid_states = 1:5
) %>%
  mutate(
    total_proteins = total_proteins,
    n_proteins = sapply(
      min_valid_states,
      function(x) sum(protein_valid_states$n_valid_states >= x)
    ),
    percent_proteins = round(100 * n_proteins / total_proteins, 2)
  )

print(qc_proteins_by_valid_states)

# Plot the result.
# Interpretation:
# If the numbers drop very fast from 1 to 5, then most proteins only have
# enough CpGs in very few chromatin states.
p_qc_proteins_by_valid_states <- ggplot(
  qc_proteins_by_valid_states,
  aes(x = factor(min_valid_states), y = n_proteins)
) +
  geom_col() +
  geom_text(
    aes(label = paste0(n_proteins, " / ", total_proteins)),
    vjust = -0.3,
    size = 4
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = paste0(
      "Proteins with enough CpGs in both groups, min CpGs = ",
      min_cpg,
      " | total proteins = ",
      total_proteins
    ),
    x = "Minimum number of valid chromatin states",
    y = "Number of proteins"
  )

print(p_qc_proteins_by_valid_states)



# # ============================================================
# # QC 2: Threshold curve
# # How many proteins have enough CpGs in both groups
# # in at least one chromatin state?
# # ============================================================
# 
# # We test many possible minimum CpG thresholds.
# # For each threshold x, a protein is counted if it has at least one
# # chromatin state where both n_1 and n_2 are >= x.
# min_cpg_grid <- seq(1, 300, by = 1)
# 
# qc_threshold_curve <- tibble(
#   min_cpg = min_cpg_grid
# ) %>%
#   rowwise() %>%
#   mutate(
#     n_proteins = stratified_test %>%
#       mutate(
#         chromatin_state = str_remove(data_set, "state_"),
#         enough_both = !is.na(n_1) & !is.na(n_2) & n_1 >= min_cpg & n_2 >= min_cpg
#       ) %>%
#       filter(enough_both) %>%
#       distinct(protein, chromatin_state) %>%
#       distinct(protein) %>%
#       nrow()
#   ) %>%
#   ungroup()
# 
# print(qc_threshold_curve)
# 
# # Plot:
# # The best cutoff is usually not where the curve is still falling sharply.
# # A useful cutoff is often around the point where the curve becomes flatter.
# p_qc_threshold_curve <- ggplot(
#   qc_threshold_curve,
#   aes(x = min_cpg, y = n_proteins)
# ) +
#   geom_line(linewidth = 1) +
#   theme_bw(base_size = 14) +
#   labs(
#     title = "QC: proteins retained across minimum CpG thresholds",
#     subtitle = "Protein counted if both n_1 and n_2 pass the threshold in at least one chromatin state",
#     x = "Minimum CpGs required in both groups",
#     y = "Number of retained proteins"
#   )

print(p_qc_threshold_curve)
