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
  "permutation_test_v4"
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

input_rds <- "stratified_test_all_pairs_1e+05perm.rds"

stratified_test_df <- readRDS(file.path(input_folder,input_rds))

str(stratified_test_df)
# ============================================================
# QC 1: How many proteins have enough CpGs in both groups
#       in at least 1, 2, 3, 4, or 5 chromatin states?
# ============================================================

library(dplyr)
library(stringr)

alpha <- 0.05

stratified_test_qc <- stratified_test_df %>%
  mutate(
    # A row is considered tested only if the statistic and p-values exist.
    tested = test_status == "tested" &
      !is.na(observed_S_statistic) &
      !is.na(p_value_S_bigger) &
      !is.na(p_value_S_smaller),
    
    # Direction of the observed S statistic.
    #
    # Interpretation:
    #   positive = binding associated with relatively higher methylation
    #   negative = binding associated with relatively lower methylation
    S_direction = case_when(
      !tested ~ NA_character_,
      observed_S_statistic > 0 ~ "positive",
      observed_S_statistic < 0 ~ "negative",
      observed_S_statistic == 0 ~ "zero"
    ),
    
    # One-sided significance in the positive direction.
    #
    # p_value_S_bigger is small when observed S is unusually large
    # compared with the permutation distribution.
    sig_positive = tested &
      observed_S_statistic > 0 &
      p_value_S_bigger < alpha,
    
    # One-sided significance in the negative direction.
    #
    # p_value_S_smaller is small when observed S is unusually small
    # compared with the permutation distribution.
    sig_negative = tested &
      observed_S_statistic < 0 &
      p_value_S_smaller < alpha,
    
    # Compact label for every row.
    sig_direction = case_when(
      !tested ~ "not_tested",
      sig_positive ~ "significant_positive",
      sig_negative ~ "significant_negative",
      S_direction == "positive" ~ "non_sig_positive",
      S_direction == "negative" ~ "non_sig_negative",
      S_direction == "zero" ~ "zero",
      TRUE ~ "check"
    )
  )

overall_direction_counts <- stratified_test_qc %>%
  count(sig_direction, name = "n_rows") %>%
  mutate(percent = round(100 * n_rows / sum(n_rows), 2))

print(overall_direction_counts)

tested_direction_counts <- stratified_test_qc %>%
  filter(tested) %>%
  count(sig_direction, name = "n_tested_rows") %>%
  mutate(percent = round(100 * n_tested_rows / sum(n_tested_rows), 2))

print(tested_direction_counts)

direction_by_state <- stratified_test_qc %>%
  filter(tested) %>%
  count(chromatin_state, sig_direction, name = "n") %>%
  group_by(chromatin_state) %>%
  mutate(percent = round(100 * n / sum(n), 2)) %>%
  ungroup()

print(direction_by_state)

p_direction_by_state <- ggplot(
  direction_by_state,
  aes(x = factor(chromatin_state), y = n, fill = sig_direction)
) +
  geom_col(position = "stack") +
  theme_bw(base_size = 14) +
  labs(
    title = paste0("Permutation-test direction by chromatin state, alpha = ", alpha),
    x = "Chromatin state",
    y = "Number of tested rows",
    fill = "Direction"
  )

print(p_direction_by_state)


consistency_by_motif_state <- stratified_test_qc %>%
  group_by(protein, motif, chromatin_state) %>%
  summarise(
    # Number of rows/pairs in this protein-motif-state group
    n_total_pairs = n(),
    n_tested_pairs = sum(tested),
    
    # Significant directions
    n_sig_positive = sum(sig_positive, na.rm = TRUE),
    n_sig_negative = sum(sig_negative, na.rm = TRUE),
    n_sig_total = n_sig_positive + n_sig_negative,
    
    # Non-significant observed directions
    n_non_sig_positive = sum(tested & !sig_positive & !sig_negative &
                               S_direction == "positive", na.rm = TRUE),
    n_non_sig_negative = sum(tested & !sig_positive & !sig_negative &
                               S_direction == "negative", na.rm = TRUE),
    
    # All observed directions, significant or not
    n_observed_positive = sum(tested & S_direction == "positive", na.rm = TRUE),
    n_observed_negative = sum(tested & S_direction == "negative", na.rm = TRUE),
    
    # Median group sizes for QC
    median_n_sample1 = median(n_sample1, na.rm = TRUE),
    median_n_sample2 = median(n_sample2, na.rm = TRUE),
    median_n_both = median(n_both, na.rm = TRUE),
    
    min_n_sample1 = min(n_sample1, na.rm = TRUE),
    min_n_sample2 = min(n_sample2, na.rm = TRUE),
    min_n_both = min(n_both, na.rm = TRUE),
    
    # Median absolute effect size
    median_abs_S = median(abs(observed_S_statistic), na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    # Significant-direction consistency:
    # This answers your main question.
    significant_consistency = case_when(
      n_tested_pairs == 0 ~ "no_tested_pairs",
      n_sig_total == 0 ~ "no_significant_pairs",
      n_sig_positive > 0 & n_sig_negative == 0 ~ "consistent_significant_positive",
      n_sig_negative > 0 & n_sig_positive == 0 ~ "consistent_significant_negative",
      n_sig_positive > 0 & n_sig_negative > 0 ~ "conflicting_significant_directions",
      TRUE ~ "check"
    ),
    
    # Observed-direction consistency:
    # This is stricter because it also counts non-significant direction flips.
    observed_direction_consistency = case_when(
      n_tested_pairs == 0 ~ "no_tested_pairs",
      n_observed_positive > 0 & n_observed_negative == 0 ~ "all_observed_positive",
      n_observed_negative > 0 & n_observed_positive == 0 ~ "all_observed_negative",
      n_observed_positive > 0 & n_observed_negative > 0 ~ "mixed_observed_directions",
      TRUE ~ "check"
    )
  )

print(consistency_by_motif_state)


conflicting_significant_cases <- consistency_by_motif_state %>%
  filter(significant_consistency == "conflicting_significant_directions") %>%
  arrange(desc(n_sig_total), protein, motif, chromatin_state)

print(conflicting_significant_cases)

conflicting_rows <- stratified_test_qc %>%
  semi_join(
    conflicting_significant_cases,
    by = c("protein", "motif", "chromatin_state")
  ) %>%
  arrange(
    protein,
    motif,
    chromatin_state,
    sig_direction,
    p_value_two_sided
  ) %>%
  select(
    protein,
    motif,
    chromatin_state,
    biosample1,
    experiment_id1,
    biosample2,
    experiment_id2,
    n_sample1,
    n_sample2,
    n_both,
    observed_S_statistic,
    p_value_S_bigger,
    p_value_S_smaller,
    p_value_two_sided,
    sig_direction,
    test_status
  )

print(conflicting_rows)


mixed_observed_cases <- consistency_by_motif_state %>%
  filter(
    significant_consistency %in% c(
      "consistent_significant_positive",
      "consistent_significant_negative"
    ),
    observed_direction_consistency == "mixed_observed_directions"
  ) %>%
  arrange(desc(n_sig_total), protein, motif, chromatin_state)

print(mixed_observed_cases)

consistency_summary <- consistency_by_motif_state %>%
  count(significant_consistency, observed_direction_consistency, name = "n_motif_states") %>%
  arrange(desc(n_motif_states))

print(consistency_summary)

consistency_by_state <- consistency_by_motif_state %>%
  count(chromatin_state, significant_consistency, name = "n_motif_states") %>%
  group_by(chromatin_state) %>%
  mutate(percent = round(100 * n_motif_states / sum(n_motif_states), 2)) %>%
  ungroup()

print(consistency_by_state)

p_consistency_by_state <- ggplot(
  consistency_by_state,
  aes(x = factor(chromatin_state), y = n_motif_states, fill = significant_consistency)
) +
  geom_col() +
  theme_bw(base_size = 14) +
  labs(
    title = "Consistency of significant direction by chromatin state",
    x = "Chromatin state",
    y = "Number of protein-motif-state groups",
    fill = "Consistency"
  )

print(p_consistency_by_state)

#stratified_test_df

#pairwise_jaccard_values_df <- read.csv(file.path(this.dir,"jaccard_analysis/Jaccard_grouped_analysis/ALL_pairwise_jaccard_values_with_groups.csv"))

all_jaccard_pairs_with_qc <- readRDS(file = file.path(this.dir,
  "QC_motif_CpG_antibody_jaccard",
  "all_jaccard_pairs_with_qc.rds"
))



library(dplyr)

# ============================================================
# Extract experiment-level no-motif QC values
# from pairwise Jaccard/QC table
# ============================================================

experiment_no_motif_long <- bind_rows(
  
  # Side 1 of pair
  all_jaccard_pairs_with_qc %>%
    transmute(
      protein,
      motif,
      biosample = biosample_1,
      experiment_id = experiment_id_1,
      antibody_label = antibody_label_1,
      No_motif_in_peak_Percentage = No_motif_in_peak_Percentage_1,
      No_CG_in_motif_Percentage = No_CG_in_motif_Percentage_1,
      No_State_assignment_for_CG_Percentage = No_State_assignment_for_CG_Percentage_1,
      Usable_1_18_Count = Usable_1_18_Count_1,
      Usable_1_18_Percentage = Usable_1_18_Percentage_1
    ),
  
  # Side 2 of pair
  all_jaccard_pairs_with_qc %>%
    transmute(
      protein,
      motif,
      biosample = biosample_2,
      experiment_id = experiment_id_2,
      antibody_label = antibody_label_2,
      No_motif_in_peak_Percentage = No_motif_in_peak_Percentage_2,
      No_CG_in_motif_Percentage = No_CG_in_motif_Percentage_2,
      No_State_assignment_for_CG_Percentage = No_State_assignment_for_CG_Percentage_2,
      Usable_1_18_Count = Usable_1_18_Count_2,
      Usable_1_18_Percentage = Usable_1_18_Percentage_2
    )
)

# Remove repeated appearances of the same experiment caused by pairwise comparisons.
experiment_no_motif_unique <- experiment_no_motif_long %>%
  distinct(
    protein,
    motif,
    biosample,
    experiment_id,
    antibody_label,
    No_motif_in_peak_Percentage,
    No_CG_in_motif_Percentage,
    No_State_assignment_for_CG_Percentage,
    Usable_1_18_Count,
    Usable_1_18_Percentage
  )

print(experiment_no_motif_unique)

experiment_no_motif_conflicts <- experiment_no_motif_unique %>%
  group_by(protein, motif, biosample, experiment_id) %>%
  summarise(
    n_distinct_no_motif_values = n_distinct(No_motif_in_peak_Percentage),
    values = paste(sort(unique(No_motif_in_peak_Percentage)), collapse = "; "),
    .groups = "drop"
  ) %>%
  filter(n_distinct_no_motif_values > 1)

print(experiment_no_motif_conflicts)

experiment_no_motif_clean <- experiment_no_motif_unique %>%
  group_by(protein, motif, biosample, experiment_id, antibody_label) %>%
  summarise(
    No_motif_in_peak_Percentage = median(No_motif_in_peak_Percentage, na.rm = TRUE),
    No_CG_in_motif_Percentage = median(No_CG_in_motif_Percentage, na.rm = TRUE),
    No_State_assignment_for_CG_Percentage = median(No_State_assignment_for_CG_Percentage, na.rm = TRUE),
    Usable_1_18_Count = median(Usable_1_18_Count, na.rm = TRUE),
    Usable_1_18_Percentage = median(Usable_1_18_Percentage, na.rm = TRUE),
    .groups = "drop"
  )


best_motif_per_protein <- experiment_no_motif_unique %>%
  group_by(protein, motif) %>%
  summarise(
    n_experiments = n_distinct(experiment_id),
    n_biosamples = n_distinct(biosample),
    median_no_motif = median(No_motif_in_peak_Percentage, na.rm = TRUE),
    mean_no_motif = mean(No_motif_in_peak_Percentage, na.rm = TRUE),
    max_no_motif = max(No_motif_in_peak_Percentage, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(protein) %>%
  arrange(
    median_no_motif,
    mean_no_motif,
    max_no_motif,
    desc(n_biosamples),
    desc(n_experiments),
    motif,
    .by_group = TRUE
  ) %>%
  slice(1) %>%
  ungroup()

print(best_motif_per_protein)

# Save selected best motif per protein
saveRDS(
  best_motif_per_protein,
  file = file.path(output_folder, "best_motif_per_protein_by_median_no_motif.rds")
)

write.csv(
  best_motif_per_protein,
  file = file.path(output_folder, "best_motif_per_protein_by_median_no_motif.csv"),
  row.names = FALSE
)

if (require(openxlsx)) {
  openxlsx::write.xlsx(
    best_motif_per_protein,
    file = file.path(output_folder, "best_motif_per_protein_by_median_no_motif.xlsx"),
    overwrite = TRUE
  )
}

## ============================================================
# Join permutation-test results with Jaccard/QC table
# ONLY for selected best motif per protein
# ============================================================
#
# Required object:
#   best_motif_per_protein
#
# Must contain:
#   protein
#   motif
#
# This keeps only one motif per protein before joining.

library(dplyr)

# -----------------------------
# 0. Keep only selected good motifs
# -----------------------------

selected_good_motifs <- best_motif_per_protein %>%
  select(protein, motif) %>%
  distinct()

cat("Selected protein-motif pairs:", nrow(selected_good_motifs), "\n")

# Filter permutation-test table to selected motifs only
stratified_test_df_good_motifs <- stratified_test_df %>%
  semi_join(
    selected_good_motifs,
    by = c("protein", "motif")
  )

# Filter Jaccard/QC table to selected motifs only
all_jaccard_pairs_with_qc_good_motifs <- all_jaccard_pairs_with_qc %>%
  semi_join(
    selected_good_motifs,
    by = c("protein", "motif")
  )

cat("Rows in stratified_test_df before motif filter:", nrow(stratified_test_df), "\n")
cat("Rows in stratified_test_df after motif filter: ", nrow(stratified_test_df_good_motifs), "\n")

cat("Rows in all_jaccard_pairs_with_qc before motif filter:", nrow(all_jaccard_pairs_with_qc), "\n")
cat("Rows in all_jaccard_pairs_with_qc after motif filter: ", nrow(all_jaccard_pairs_with_qc_good_motifs), "\n")


# -----------------------------
# 1. Keep only essential Jaccard/QC columns
# -----------------------------

jaccard_forward <- all_jaccard_pairs_with_qc_good_motifs %>%
  transmute(
    protein,
    motif,
    
    biosample1 = biosample_1,
    experiment_id1 = experiment_id_1,
    antibody_label1 = antibody_label_1,
    
    biosample2 = biosample_2,
    experiment_id2 = experiment_id_2,
    antibody_label2 = antibody_label_2,
    
    jaccard_index,
    same_biosample,
    same_antibody,
    combined_group,
    
    usable_count1 = Usable_1_18_Count_1,
    usable_percentage1 = Usable_1_18_Percentage_1,
    no_motif_percentage1 = No_motif_in_peak_Percentage_1,
    no_cg_percentage1 = No_CG_in_motif_Percentage_1,
    no_state_assignment_percentage1 = No_State_assignment_for_CG_Percentage_1,
    
    usable_count2 = Usable_1_18_Count_2,
    usable_percentage2 = Usable_1_18_Percentage_2,
    no_motif_percentage2 = No_motif_in_peak_Percentage_2,
    no_cg_percentage2 = No_CG_in_motif_Percentage_2,
    no_state_assignment_percentage2 = No_State_assignment_for_CG_Percentage_2,
    
    pair_low_usable_count,
    pair_high_no_motif,
    pair_low_usable_percentage,
    pair_high_no_cg,
    pair_QC_reason
  )

# Reversed version:
# side 1 and side 2 are swapped so that the join aligns correctly
# even if the pair appears in the opposite order.
jaccard_reversed <- all_jaccard_pairs_with_qc_good_motifs %>%
  transmute(
    protein,
    motif,
    
    biosample1 = biosample_2,
    experiment_id1 = experiment_id_2,
    antibody_label1 = antibody_label_2,
    
    biosample2 = biosample_1,
    experiment_id2 = experiment_id_1,
    antibody_label2 = antibody_label_1,
    
    jaccard_index,
    same_biosample,
    same_antibody,
    combined_group,
    
    usable_count1 = Usable_1_18_Count_2,
    usable_percentage1 = Usable_1_18_Percentage_2,
    no_motif_percentage1 = No_motif_in_peak_Percentage_2,
    no_cg_percentage1 = No_CG_in_motif_Percentage_2,
    no_state_assignment_percentage1 = No_State_assignment_for_CG_Percentage_2,
    
    usable_count2 = Usable_1_18_Count_1,
    usable_percentage2 = Usable_1_18_Percentage_1,
    no_motif_percentage2 = No_motif_in_peak_Percentage_1,
    no_cg_percentage2 = No_CG_in_motif_Percentage_1,
    no_state_assignment_percentage2 = No_State_assignment_for_CG_Percentage_1,
    
    pair_low_usable_count,
    pair_high_no_motif,
    pair_low_usable_percentage,
    pair_high_no_cg,
    pair_QC_reason
  )

# Combine forward and reversed versions
jaccard_join_table <- bind_rows(
  jaccard_forward,
  jaccard_reversed
) %>%
  distinct(
    protein,
    motif,
    biosample1,
    experiment_id1,
    biosample2,
    experiment_id2,
    .keep_all = TRUE
  )


# ----------------------------------------------------------
# 2. Keep only essential columns from permutation results
# ----------------------------------------------------------

stratified_test_compact <- stratified_test_df_good_motifs %>%
  select(
    protein,
    motif,
    chromatin_state,
    
    biosample1,
    experiment_id1,
    
    biosample2,
    experiment_id2,
    
    n_sample1,
    n_sample2,
    n_both,
    
    observed_S_statistic,
    p_value_S_bigger,
    p_value_S_smaller,
    p_value_two_sided,
    
    test_status
  )


# -----------------------------
# 3. Join permutation + Jaccard/QC table
# -----------------------------

stratified_test_with_jaccard_qc <- stratified_test_compact %>%
  left_join(
    jaccard_join_table,
    by = c(
      "protein",
      "motif",
      "biosample1",
      "experiment_id1",
      "biosample2",
      "experiment_id2"
    )
  )


# -----------------------------
# 4. Check join quality
# -----------------------------

cat("Rows before join:", nrow(stratified_test_compact), "\n")
cat("Rows after join: ", nrow(stratified_test_with_jaccard_qc), "\n")

cat(
  "Rows without Jaccard/QC match:",
  sum(is.na(stratified_test_with_jaccard_qc$jaccard_index)),
  "\n"
)

unmatched_jaccard_rows <- stratified_test_with_jaccard_qc %>%
  filter(is.na(jaccard_index)) %>%
  distinct(
    protein,
    motif,
    biosample1,
    experiment_id1,
    biosample2,
    experiment_id2
  )

print(unmatched_jaccard_rows)

# shoudl retrurn empy 

cat("Unique proteins:\n")
print(
  stratified_test_with_jaccard_qc %>%
    summarise(n_proteins = n_distinct(protein))
)

cat("\nUnique protein-motif pairs:\n")
print(
  stratified_test_with_jaccard_qc %>%
    distinct(protein, motif) %>%
    summarise(n_protein_motif_pairs = n())
)

min_n_per_group <- 50

stratified_test_with_jaccard_qc_n50 <- stratified_test_with_jaccard_qc %>%
  filter(
    n_sample1 >= min_n_per_group,
    n_sample2 >= min_n_per_group
  )

cat("Rows before n50 filter:", nrow(stratified_test_with_jaccard_qc), "\n")
cat("Rows after n50 filter: ", nrow(stratified_test_with_jaccard_qc_n50), "\n")

cat("\nUnique protein-motif pairs:\n")
print(
  stratified_test_with_jaccard_qc_n50 %>%
    distinct(protein, motif) %>%
    summarise(n_protein_motif_pairs = n())
)

# -----------------------------
# 5. Save joined table
# -----------------------------

saveRDS(
  stratified_test_with_jaccard_qc,
  file = file.path(output_folder, "stratified_test_with_jaccard_qc_best_motif_only.rds")
)

write.csv(
  stratified_test_with_jaccard_qc,
  file = file.path(output_folder, "stratified_test_with_jaccard_qc_best_motif_only.csv"),
  row.names = FALSE
)

if (require(openxlsx)) {
  openxlsx::write.xlsx(
    stratified_test_with_jaccard_qc,
    file = file.path(output_folder, "stratified_test_with_jaccard_qc_best_motif_only.xlsx"),
    overwrite = TRUE
  )
}

##################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

alpha <- 0.05

qc_output_folder <- file.path(output_folder, "QC_significance_vs_pair_QC")

if (!dir.exists(qc_output_folder)) {
  dir.create(qc_output_folder, recursive = TRUE)
}


library(dplyr)
library(ggplot2)

alpha <- 0.05

# ============================================================
# Direction and consistency plots for cleaned n50 table
# ============================================================

stratified_test_qc_n50 <- stratified_test_with_jaccard_qc_n50 %>%
  mutate(
    tested = test_status == "tested" &
      !is.na(observed_S_statistic) &
      !is.na(p_value_S_bigger) &
      !is.na(p_value_S_smaller),
    
    S_direction = case_when(
      !tested ~ NA_character_,
      observed_S_statistic > 0 ~ "positive",
      observed_S_statistic < 0 ~ "negative",
      observed_S_statistic == 0 ~ "zero",
      TRUE ~ NA_character_
    ),
    
    sig_positive = tested &
      observed_S_statistic > 0 &
      p_value_S_bigger < alpha,
    
    sig_negative = tested &
      observed_S_statistic < 0 &
      p_value_S_smaller < alpha,
    
    sig_direction = case_when(
      !tested ~ "not_tested",
      sig_positive ~ "significant_positive",
      sig_negative ~ "significant_negative",
      S_direction == "positive" ~ "non_sig_positive",
      S_direction == "negative" ~ "non_sig_negative",
      S_direction == "zero" ~ "zero",
      TRUE ~ "check"
    )
  )

direction_by_state_n50 <- stratified_test_qc_n50 %>%
  filter(tested) %>%
  count(chromatin_state, sig_direction, name = "n") %>%
  group_by(chromatin_state) %>%
  mutate(percent = round(100 * n / sum(n), 2)) %>%
  ungroup()

print(direction_by_state_n50)

p_direction_by_state_n50 <- ggplot(
  direction_by_state_n50,
  aes(x = factor(chromatin_state), y = n, fill = sig_direction)
) +
  geom_col(position = "stack") +
  theme_bw(base_size = 14) +
  labs(
    title = paste0("Permutation-test direction by chromatin state, n ≥ 50, alpha = ", alpha),
    x = "Chromatin state",
    y = "Number of tested rows",
    fill = "Direction"
  )

print(p_direction_by_state_n50)

consistency_by_motif_state_n50 <- stratified_test_qc_n50 %>%
  group_by(protein, motif, chromatin_state) %>%
  summarise(
    n_total_pairs = n(),
    n_tested_pairs = sum(tested),
    
    n_sig_positive = sum(sig_positive, na.rm = TRUE),
    n_sig_negative = sum(sig_negative, na.rm = TRUE),
    n_sig_total = n_sig_positive + n_sig_negative,
    
    n_observed_positive = sum(tested & S_direction == "positive", na.rm = TRUE),
    n_observed_negative = sum(tested & S_direction == "negative", na.rm = TRUE),
    
    median_n_sample1 = median(n_sample1, na.rm = TRUE),
    median_n_sample2 = median(n_sample2, na.rm = TRUE),
    median_n_both = median(n_both, na.rm = TRUE),
    median_abs_S = median(abs(observed_S_statistic), na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    significant_consistency = case_when(
      n_tested_pairs == 0 ~ "no_tested_pairs",
      n_sig_total == 0 ~ "no_significant_pairs",
      n_sig_positive > 0 & n_sig_negative == 0 ~ "consistent_significant_positive",
      n_sig_negative > 0 & n_sig_positive == 0 ~ "consistent_significant_negative",
      n_sig_positive > 0 & n_sig_negative > 0 ~ "conflicting_significant_directions",
      TRUE ~ "check"
    ),
    
    observed_direction_consistency = case_when(
      n_tested_pairs == 0 ~ "no_tested_pairs",
      n_observed_positive > 0 & n_observed_negative == 0 ~ "all_observed_positive",
      n_observed_negative > 0 & n_observed_positive == 0 ~ "all_observed_negative",
      n_observed_positive > 0 & n_observed_negative > 0 ~ "mixed_observed_directions",
      TRUE ~ "check"
    )
  )

print(consistency_by_motif_state_n50)


consistency_by_state_n50 <- consistency_by_motif_state_n50 %>%
  count(chromatin_state, significant_consistency, name = "n_motif_states") %>%
  group_by(chromatin_state) %>%
  mutate(percent = round(100 * n_motif_states / sum(n_motif_states), 2)) %>%
  ungroup()

print(consistency_by_state_n50)

p_consistency_by_state_n50 <- ggplot(
  consistency_by_state_n50,
  aes(x = factor(chromatin_state), y = n_motif_states, fill = significant_consistency)
) +
  geom_col() +
  theme_bw(base_size = 14) +
  labs(
    title = "Consistency of significant direction by chromatin state, n ≥ 50",
    x = "Chromatin state",
    y = "Number of protein-motif-state groups",
    fill = "Consistency"
  )

print(p_consistency_by_state_n50)

conflicting_significant_cases_n50 <- consistency_by_motif_state_n50 %>%
  filter(significant_consistency == "conflicting_significant_directions") %>%
  arrange(desc(n_sig_total), protein, motif, chromatin_state)

print(conflicting_significant_cases_n50)

conflicting_rows_n50 <- stratified_test_qc_n50 %>%
  semi_join(
    conflicting_significant_cases_n50,
    by = c("protein", "motif", "chromatin_state")
  ) %>%
  arrange(
    protein,
    motif,
    chromatin_state,
    sig_direction,
    p_value_two_sided
  )

print(conflicting_rows_n50)


conflicting_rows_n50 <- stratified_test_qc_n50 %>%
  semi_join(
    consistency_by_motif_state_n50 %>%
      filter(significant_consistency == "conflicting_significant_directions"),
    by = c("protein", "motif", "chromatin_state")
  )

print(conflicting_rows_n50)

length(unique(c(stratified_test_qc_n50$experiment_id1,stratified_test_qc_n50$experiment_id2)))
length(unique(c(stratified_test_qc_n50$protein,stratified_test_qc_n50$protein)))



############# median jaccard 

library(dplyr)


# 1. Remove chromatin-state duplication
unique_pairs_n50 <- stratified_test_qc_n50 %>%
  distinct(
    protein,
    motif,
    biosample1,
    experiment_id1,
    antibody_label1,
    biosample2,
    experiment_id2,
    antibody_label2,
    jaccard_index,
    .keep_all = TRUE
  )

# 2. Convert pair table into experiment-level table
experiment_jaccard_n50 <- bind_rows(
  unique_pairs_n50 %>%
    transmute(
      protein,
      motif,
      biosample = biosample1,
      experiment_id = experiment_id1,
      antibody_label = antibody_label1,
      other_biosample = biosample2,
      other_experiment_id = experiment_id2,
      jaccard_index
    ),
  
  unique_pairs_n50 %>%
    transmute(
      protein,
      motif,
      biosample = biosample2,
      experiment_id = experiment_id2,
      antibody_label = antibody_label2,
      other_biosample = biosample1,
      other_experiment_id = experiment_id1,
      jaccard_index
    )
) %>%
  distinct()

experiment_median_jaccard_n50 <- experiment_jaccard_n50 %>%
  group_by(protein, motif, biosample, experiment_id, antibody_label) %>%
  summarise(
    n_jaccard_comparisons = n(),
    median_jaccard = median(jaccard_index, na.rm = TRUE),
    mean_jaccard = mean(jaccard_index, na.rm = TRUE),
    min_jaccard = min(jaccard_index, na.rm = TRUE),
    max_jaccard = max(jaccard_index, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(protein, biosample, desc(median_jaccard))

print(experiment_median_jaccard_n50)

best_experiment_by_median_jaccard_n50 <- experiment_median_jaccard_n50 %>%
  group_by(protein, motif, biosample) %>%
  arrange(
    desc(median_jaccard),
    desc(mean_jaccard),
    desc(n_jaccard_comparisons),
    experiment_id,
    .by_group = TRUE
  ) %>%
  slice(1) %>%
  ungroup()

print(best_experiment_by_median_jaccard_n50)

library(dplyr)

# Best experiment table: one row per protein + motif + biosample
best_exp_keys <- best_experiment_by_median_jaccard_n50 %>%
  select(protein, motif, biosample, experiment_id) %>%
  distinct()

# Keep only rows where experiment_id1 is best for biosample1
# AND experiment_id2 is best for biosample2
stratified_test_qc_n50_best_experiments <- stratified_test_qc_n50 %>%
  semi_join(
    best_exp_keys %>%
      rename(
        biosample1 = biosample,
        experiment_id1 = experiment_id
      ),
    by = c("protein", "motif", "biosample1", "experiment_id1")
  ) %>%
  semi_join(
    best_exp_keys %>%
      rename(
        biosample2 = biosample,
        experiment_id2 = experiment_id
      ),
    by = c("protein", "motif", "biosample2", "experiment_id2")
  )

cat("Rows before best-experiment filter:", nrow(stratified_test_qc_n50), "\n")
cat("Rows after best-experiment filter: ", nrow(stratified_test_qc_n50_best_experiments), "\n")

summary_best_exp_n50 <- stratified_test_qc_n50_best_experiments %>%
  summarise(
    n_rows = n(),
    n_proteins = n_distinct(protein),
    n_motifs = n_distinct(motif),
    n_protein_motif_pairs = n_distinct(paste(protein, motif)),
    n_biosamples = n_distinct(c(biosample1, biosample2)),
    n_experiment_ids = n_distinct(c(experiment_id1, experiment_id2)),
    n_experiment_pairs = n_distinct(
      paste(protein, motif, biosample1, experiment_id1, biosample2, experiment_id2)
    ),
    n_chromatin_states = n_distinct(chromatin_state)
  )

print(summary_best_exp_n50)


chromatin_state_counts_best_exp <- stratified_test_qc_n50_best_experiments %>%
  count(chromatin_state, name = "n_rows") %>%
  tidyr::complete(
    chromatin_state = 1:18,
    fill = list(n_rows = 0)
  ) %>%
  mutate(
    percent = round(100 * n_rows / sum(n_rows), 2)
  ) %>%
  arrange(chromatin_state)

p_chromatin_state_counts_best_exp <- ggplot(
  chromatin_state_counts_best_exp,
  aes(x = factor(chromatin_state), y = n_rows, fill = factor(chromatin_state))
) +
  geom_col() +
  geom_text(
    aes(label = paste0(n_rows, "\n", percent, "%")),
    vjust = -0.25,
    size = 3.5
  ) +
  scale_fill_manual(
    values = chromatin_state_colors_short,
    drop = FALSE
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Chromatin states represented after best-experiment filtering",
    x = "Chromatin state",
    y = "Number of rows",
    fill = "Chromatin state"
  )

print(p_chromatin_state_counts_best_exp)


library(dplyr)

alpha_fdr <- 0.05

# ============================================================
# Add FDR-corrected two-sided p-values
# and redefine significance direction by FDR
# ============================================================
# stratified_test_qc_n50_best_experiments <- stratified_test_qc_n50_best_experiments_FDR
stratified_test_qc_n50_best_experiments <- stratified_test_qc_n50_best_experiments %>%
  mutate(
    # Keep old direction label for comparison
    sig_direction_nominal_old = sig_direction,

    tested = test_status == "tested" &
      !is.na(observed_S_statistic) &
      !is.na(p_value_two_sided)
  )

# Add BH/FDR correction only among actually tested rows
stratified_test_qc_n50_best_experiments$q_value_two_sided <- NA_real_

stratified_test_qc_n50_best_experiments$q_value_two_sided[
  stratified_test_qc_n50_best_experiments$tested
] <- p.adjust(
  stratified_test_qc_n50_best_experiments$p_value_two_sided[
    stratified_test_qc_n50_best_experiments$tested
  ],
  method = "BH"
)

# Redefine significance using FDR-corrected two-sided p-value
stratified_test_qc_n50_best_experiments <- stratified_test_qc_n50_best_experiments %>%
  mutate(
    sig_positive = tested &
      observed_S_statistic > 0 &
      q_value_two_sided < alpha_fdr,

    sig_negative = tested &
      observed_S_statistic < 0 &
      q_value_two_sided < alpha_fdr,

    sig_any = sig_positive | sig_negative,

    sig_direction = case_when(
      !tested ~ "not_tested",
      sig_positive ~ "significant_positive",
      sig_negative ~ "significant_negative",
      tested & observed_S_statistic > 0 ~ "non_sig_positive",
      tested & observed_S_statistic < 0 ~ "non_sig_negative",
      tested & observed_S_statistic == 0 ~ "zero",
      TRUE ~ "check"
    )
  )

stratified_test_qc_n50_best_experiments %>%
  count(sig_direction, name = "n_rows") %>%
  mutate(percent = round(100 * n_rows / sum(n_rows), 2)) %>%
  print()


saveRDS(
  stratified_test_qc_n50_best_experiments,
  file = file.path(output_folder, "stratified_test_qc_n50_best_experiments_FDR.rds")
)

write.csv(
  stratified_test_qc_n50_best_experiments,
  file = file.path(output_folder, "stratified_test_qc_n50_best_experiments_FDR.csv"),
  row.names = FALSE
)



# Recalculate consistency after keeping only best experiments
consistency_best_exp_n50 <- stratified_test_qc_n50_best_experiments %>%
  group_by(protein, motif, chromatin_state) %>%
  summarise(
    n_total_pairs = n(),
    n_tested_pairs = sum(tested, na.rm = TRUE),

    n_sig_positive = sum(sig_positive, na.rm = TRUE),
    n_sig_negative = sum(sig_negative, na.rm = TRUE),
    n_sig_total = n_sig_positive + n_sig_negative,

    .groups = "drop"
  ) %>%
  mutate(
    significant_consistency = case_when(
      n_tested_pairs == 0 ~ "no_tested_pairs",
      n_sig_total == 0 ~ "no_significant_pairs",
      n_sig_positive > 0 & n_sig_negative == 0 ~ "consistent_significant_positive",
      n_sig_negative > 0 & n_sig_positive == 0 ~ "consistent_significant_negative",
      n_sig_positive > 0 & n_sig_negative > 0 ~ "conflicting_significant_directions",
      TRUE ~ "check"
    )
  )

# Summary
consistency_best_exp_n50 %>%
  count(significant_consistency, name = "n_motif_states") %>%
  arrange(desc(n_motif_states))


conflicting_best_exp_n50 <- consistency_best_exp_n50 %>%
  filter(significant_consistency == "conflicting_significant_directions") %>%
  arrange(desc(n_sig_total), protein, chromatin_state)

print(conflicting_best_exp_n50)



# cheklc results
protein_sig_summary_FDR <- stratified_test_qc_n50_best_experiments %>%
  group_by(protein) %>%
  summarise(
    n_rows = n(),
    n_tested = sum(tested, na.rm = TRUE),

    n_sig_positive = sum(sig_direction == "significant_positive", na.rm = TRUE),
    n_sig_negative = sum(sig_direction == "significant_negative", na.rm = TRUE),
    n_sig_total = n_sig_positive + n_sig_negative,

    n_non_sig_positive = sum(sig_direction == "non_sig_positive", na.rm = TRUE),
    n_non_sig_negative = sum(sig_direction == "non_sig_negative", na.rm = TRUE),

    percent_sig_positive = round(100 * n_sig_positive / n_tested, 2),
    percent_sig_negative = round(100 * n_sig_negative / n_tested, 2),
    percent_sig_total = round(100 * n_sig_total / n_tested, 2),

    dominant_sig_direction = case_when(
      n_sig_positive > 0 & n_sig_negative == 0 ~ "positive_only",
      n_sig_negative > 0 & n_sig_positive == 0 ~ "negative_only",
      n_sig_positive > 0 & n_sig_negative > 0 ~ "both_directions",
      n_sig_total == 0 ~ "no_significant",
      TRUE ~ "check"
    ),

    .groups = "drop"
  ) %>%
  arrange(desc(n_sig_total), desc(abs(n_sig_positive - n_sig_negative)), protein)

print(protein_sig_summary_FDR)



contradictory_proteins_FDR <- protein_sig_summary_FDR %>%
  filter(dominant_sig_direction == "both_directions") %>%
  arrange(desc(n_sig_total), protein)

print(contradictory_proteins_FDR)

contradictory_rows_FDR <- stratified_test_qc_n50_best_experiments %>%
  semi_join(
    contradictory_proteins_FDR %>% select(protein),
    by = "protein"
  ) %>%
  filter(sig_direction %in% c("significant_positive", "significant_negative")) %>%
  arrange(protein, chromatin_state, sig_direction, q_value_two_sided) %>%
  select(
    protein,
    motif,
    chromatin_state,
    biosample1,
    experiment_id1,
    biosample2,
    experiment_id2,
    n_sample1,
    n_sample2,
    observed_S_statistic,
    p_value_two_sided,
    q_value_two_sided,
    sig_direction
  )

print(contradictory_rows_FDR)




r_table <- stratified_test_qc_n50_best_experiments %>%
  count(sig_direction, name = "n_rows") %>%
  mutate(percent = round(100 * n_rows / sum(n_rows), 2))


n_sig_negative <- r_table$n_rows[3]
n_sig_positive <- r_table$n_rows[4]
alpha_fdr <- 0.05

expected_false_negative <- alpha_fdr * n_sig_negative
expected_false_positive <- alpha_fdr * n_sig_positive
expected_false_total <- alpha_fdr * (n_sig_negative + n_sig_positive)

expected_false_negative
expected_false_positive
expected_false_total

saveRDS(object =  stratified_test_qc_n50_best_experiments,
  file = file.path(output_folder, "stratified_test_qc_n50_best_experiments_FDR.rds")
)

write.csv(
  stratified_test_qc_n50_best_experiments,
  file = file.path(output_folder, "stratified_test_qc_n50_best_experiments_FDR.csv"),
  row.names = FALSE
)

if (require(openxlsx)) {
  openxlsx::write.xlsx(
    stratified_test_qc_n50_best_experiments,
    file = file.path(output_folder, "stratified_test_qc_n50_best_experiments_FDR.xlsx"),
    overwrite = TRUE
  )
}


library(dplyr)

experiment_counts <- dplyr::bind_rows(
  stratified_test_qc_n50_best_experiments %>%
    transmute(protein, ExperimentID = experiment_id1, Biosample = biosample1),
  stratified_test_qc_n50_best_experiments %>%
    transmute(protein, ExperimentID = experiment_id2, Biosample = biosample2)
) %>%
  distinct() %>%
  summarise(
    n_proteins = n_distinct(protein),
    n_experiments = n_distinct(ExperimentID),
    n_biosamples = n_distinct(Biosample),
    n_protein_experiment_pairs = n_distinct(paste(protein, ExperimentID, sep = "_"))
  )

experiment_counts


library(dplyr)

state1_results <- stratified_test_qc_n50_best_experiments %>%
  filter(chromatin_state == 1, tested == TRUE)

state1_protein_summary <- state1_results %>%
  group_by(protein) %>%
  summarise(
    n_tests = n(),
    n_sig_positive = sum(sig_positive, na.rm = TRUE),
    n_sig_negative = sum(sig_negative, na.rm = TRUE),
    n_sig_any = sum(sig_any, na.rm = TRUE),
    has_sig_positive = any(sig_positive, na.rm = TRUE),
    has_sig_negative = any(sig_negative, na.rm = TRUE),
    has_sig_any = any(sig_any, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    protein_result = case_when(
      has_sig_positive & has_sig_negative ~ "significant_both_directions",
      has_sig_positive ~ "significant_higher_methylation",
      has_sig_negative ~ "significant_lower_methylation",
      TRUE ~ "not_significant"
    )
  )

state1_protein_counts <- state1_protein_summary %>%
  count(protein_result)

state1_protein_counts


stratified_test_qc_n50_best_experiments[stratified_test_qc_n50_best_experiments$sig_direction == "significant_positive",]
unique(stratified_test_qc_n50_best_experiments$protein[stratified_test_qc_n50_best_experiments$sig_direction == "significant_positive"])
# CEBPB  significant_positive
# outlier: MAFF MAFF.H12CORE.0.PSM.A              17      HepG2    ENCFF643LNZ       K562    ENCFF939WDM       326        76          0.128511607           0.00250      0.0061970684 significant_positive
# REST significant_negative  outlier:  REST.H12CORE.0.P.B              18       A549    ENCFF871AVP    GM12878    ENCFF348LKE       147       150          0.118382048           0.00606      0.0137661493 significant_positive
# ZBTB40 significant_positive outlier: ZBT40.H12CORE.0.P.B               2      HepG2    ENCFF703AUF       K562    ENCFF745YDH        79        88         -0.077510311           0.02306      0.0440921106 significant_negative
# ZNF121 significant_negative outlier : ZNF121  ZN121.H12CORE.0.P.B              18      HepG2    ENCFF817UOJ       K562    ENCFF650DWZ       556       955          0.092266522           0.00006      0.0002094495 significant_positive