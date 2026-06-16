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
# Join permutation-test results with Jaccard/QC table
# ============================================================
#
# Important:
# stratified_test_df uses:
#   biosample1, experiment_id1
#   biosample2, experiment_id2
#
# all_jaccard_pairs_with_qc uses:
#   biosample_1, experiment_id_1
#   biosample_2, experiment_id_2
#
# The pair orientation may not always be identical.
# Therefore we create:
#   1. forward Jaccard table
#   2. reversed Jaccard table, where side 1 and side 2 are swapped
#
# This ensures that QC columns ending in _1 always correspond to
# stratified_test_df$biosample1 / experiment_id1,
# and QC columns ending in _2 always correspond to
# stratified_test_df$biosample2 / experiment_id2.

# -----------------------------
# 1. Keep only essential columns
# -----------------------------

jaccard_forward <- all_jaccard_pairs_with_qc %>%
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
# side 1 and side 2 are swapped so that the join can still align correctly
# if the pair appears in the opposite order.
jaccard_reversed <- all_jaccard_pairs_with_qc %>%
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

# Combine forward and reversed versions.
# distinct() avoids duplicate matches if a pair is symmetrical.
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
# 2. Keep only essential columns from the permutation results
# ----------------------------------------------------------

stratified_test_compact <- stratified_test_df %>%
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
# 3. Join permutation + QC table
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

# Show unmatched rows if any exist
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

# -----------------------------
# 5. Optional: save joined table
# -----------------------------

saveRDS(
  stratified_test_with_jaccard_qc,
  file = file.path(output_folder, "stratified_test_with_jaccard_qc.rds")
)

if (require(openxlsx)) {
  openxlsx::write.xlsx(
    stratified_test_with_jaccard_qc,
    file = file.path(output_folder, "stratified_test_with_jaccard_qc.xlsx"),
    overwrite = TRUE
  )
}


library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

alpha <- 0.05

qc_output_folder <- file.path(output_folder, "QC_significance_vs_pair_QC")

if (!dir.exists(qc_output_folder)) {
  dir.create(qc_output_folder, recursive = TRUE)
}


# ============================================================
# Add significance direction and pair-QC labels
# ============================================================

stratified_qc_direction <- stratified_test_with_jaccard_qc %>%
  mutate(
    # Row was actually tested only if S statistic and p-values exist
    tested = test_status == "tested" &
      !is.na(observed_S_statistic) &
      !is.na(p_value_S_bigger) &
      !is.na(p_value_S_smaller),
    
    # Direction of observed S statistic
    S_direction = case_when(
      !tested ~ NA_character_,
      observed_S_statistic > 0 ~ "positive",
      observed_S_statistic < 0 ~ "negative",
      observed_S_statistic == 0 ~ "zero",
      TRUE ~ NA_character_
    ),
    
    # One-sided significance in positive direction
    sig_positive = tested &
      observed_S_statistic > 0 &
      p_value_S_bigger < alpha,
    
    # One-sided significance in negative direction
    sig_negative = tested &
      observed_S_statistic < 0 &
      p_value_S_smaller < alpha,
    
    sig_any = sig_positive | sig_negative,
    
    # Compact row-level category
    row_sig_category = case_when(
      !tested ~ "not_tested",
      sig_positive ~ "significant_positive",
      sig_negative ~ "significant_negative",
      tested & !sig_any & S_direction == "positive" ~ "non_sig_positive",
      tested & !sig_any & S_direction == "negative" ~ "non_sig_negative",
      tested & !sig_any & S_direction == "zero" ~ "zero",
      TRUE ~ "check"
    ),
    
    # Pair QC flag.
    # Assumes pair_QC_reason has "Pass" for clean pairs.
    pair_QC_reason_chr = as.character(pair_QC_reason),
    
    pair_QC_flag = case_when(
      is.na(pair_QC_reason_chr) ~ NA,
      pair_QC_reason_chr == "Pass" ~ FALSE,
      TRUE ~ TRUE
    ),
    
    pair_QC_status = case_when(
      is.na(pair_QC_flag) ~ "No Jaccard/QC match",
      pair_QC_flag ~ "QC flagged comparison",
      !pair_QC_flag ~ "QC pass comparison"
    )
  )

# ============================================================
# Row-level summary:
# Among tested rows, how often are results non-significant?
# ============================================================

row_qc_summary <- stratified_qc_direction %>%
  filter(tested) %>%
  group_by(pair_QC_status) %>%
  summarise(
    n_tested_rows = n(),
    n_significant = sum(sig_any, na.rm = TRUE),
    n_non_significant = sum(!sig_any, na.rm = TRUE),
    percent_non_significant = round(100 * n_non_significant / n_tested_rows, 2),
    percent_significant = round(100 * n_significant / n_tested_rows, 2),
    .groups = "drop"
  )

print(row_qc_summary)

p_non_sig_by_pair_qc <- ggplot(
  row_qc_summary,
  aes(x = pair_QC_status, y = percent_non_significant)
) +
  geom_col() +
  geom_text(
    aes(label = paste0(percent_non_significant, "%\n", n_non_significant, "/", n_tested_rows)),
    vjust = -0.25,
    size = 4
  ) +
  theme_bw(base_size = 14) +
  labs(
    title = "Non-significant tests by pair QC status",
    subtitle = paste0("alpha = ", alpha, "; only tested rows included"),
    x = "Pair QC status",
    y = "Non-significant tests (%)"
  )

print(p_non_sig_by_pair_qc)

ggsave(
  filename = file.path(qc_output_folder, "non_significant_by_pair_QC_status.svg"),
  plot = p_non_sig_by_pair_qc,
  width = 8,
  height = 5,
  limitsize = FALSE
)

# ============================================================
# Row-level significance category by exact QC reason
# ============================================================

row_sig_by_qc_reason <- stratified_qc_direction %>%
  filter(tested) %>%
  mutate(
    pair_QC_reason_chr = ifelse(
      is.na(pair_QC_reason_chr),
      "No Jaccard/QC match",
      pair_QC_reason_chr
    )
  ) %>%
  count(pair_QC_reason_chr, row_sig_category, name = "n") %>%
  group_by(pair_QC_reason_chr) %>%
  mutate(
    total = sum(n),
    percent = round(100 * n / total, 2)
  ) %>%
  ungroup()

print(row_sig_by_qc_reason)

p_sig_category_by_qc_reason <- ggplot(
  row_sig_by_qc_reason,
  aes(x = pair_QC_reason_chr, y = percent, fill = row_sig_category)
) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = "Significance category by pair QC reason",
    subtitle = paste0("alpha = ", alpha, "; only tested rows included"),
    x = "Pair QC reason",
    y = "Rows (%)",
    fill = "Test category"
  )

print(p_sig_category_by_qc_reason)

ggsave(
  filename = file.path(qc_output_folder, "significance_category_by_pair_QC_reason.svg"),
  plot = p_sig_category_by_qc_reason,
  width = 10,
  height = 6,
  limitsize = FALSE
)

# ============================================================
# Consistency per protein + motif + chromatin state
# ============================================================

motif_state_consistency_qc <- stratified_qc_direction %>%
  group_by(protein, motif, chromatin_state) %>%
  summarise(
    n_total_rows = n(),
    n_tested_pairs = sum(tested, na.rm = TRUE),
    
    n_sig_positive = sum(sig_positive, na.rm = TRUE),
    n_sig_negative = sum(sig_negative, na.rm = TRUE),
    n_sig_total = n_sig_positive + n_sig_negative,
    
    n_non_sig = sum(tested & !sig_any, na.rm = TRUE),
    
    n_observed_positive = sum(tested & S_direction == "positive", na.rm = TRUE),
    n_observed_negative = sum(tested & S_direction == "negative", na.rm = TRUE),
    
    # QC information across all pair comparisons belonging to this motif-state
    n_pairs_with_qc_match = sum(!is.na(pair_QC_flag)),
    n_qc_flagged_pairs = sum(pair_QC_flag, na.rm = TRUE),
    n_qc_pass_pairs = sum(pair_QC_flag == FALSE, na.rm = TRUE),
    
    any_qc_flagged_pair = any(pair_QC_flag == TRUE, na.rm = TRUE),
    all_pairs_qc_pass = all(pair_QC_flag == FALSE, na.rm = TRUE),
    
    median_jaccard = median(jaccard_index, na.rm = TRUE),
    min_jaccard = min(jaccard_index, na.rm = TRUE),
    
    median_n_sample1 = median(n_sample1, na.rm = TRUE),
    median_n_sample2 = median(n_sample2, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  mutate(
    any_qc_flagged_pair = ifelse(n_pairs_with_qc_match == 0, NA, any_qc_flagged_pair),
    all_pairs_qc_pass = ifelse(n_pairs_with_qc_match == 0, NA, all_pairs_qc_pass),
    
    pair_QC_group = case_when(
      is.na(any_qc_flagged_pair) ~ "No Jaccard/QC match",
      any_qc_flagged_pair ~ "At least one QC-flagged comparison",
      !any_qc_flagged_pair ~ "All comparisons QC pass"
    ),
    
    significant_consistency = case_when(
      n_tested_pairs == 0 ~ "no_tested_pairs",
      n_sig_total == 0 ~ "no_significant_pairs",
      n_sig_positive > 0 & n_sig_negative == 0 ~ "consistent_significant_positive",
      n_sig_negative > 0 & n_sig_positive == 0 ~ "consistent_significant_negative",
      n_sig_positive > 0 & n_sig_negative > 0 ~ "conflicting_significant_directions",
      TRUE ~ "check"
    ),
    
    problem_group = significant_consistency %in% c(
      "no_significant_pairs",
      "conflicting_significant_directions"
    )
  )

print(motif_state_consistency_qc)

# ============================================================
# Motif-state level:
# Are problematic consistency outcomes enriched in QC-flagged groups?
# ============================================================

motif_state_problem_by_qc <- motif_state_consistency_qc %>%
  filter(n_tested_pairs > 0) %>%
  count(pair_QC_group, significant_consistency, name = "n_motif_states") %>%
  group_by(pair_QC_group) %>%
  mutate(
    total = sum(n_motif_states),
    percent = round(100 * n_motif_states / total, 2)
  ) %>%
  ungroup()

print(motif_state_problem_by_qc)

p_motif_state_problem_by_qc <- ggplot(
  motif_state_problem_by_qc,
  aes(x = pair_QC_group, y = percent, fill = significant_consistency)
) +
  geom_col() +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = "Motif-state consistency outcome by pair QC status",
    subtitle = "Grouped by protein + motif + chromatin state; only groups with tested pairs included",
    x = "Pair QC group",
    y = "Motif-state groups (%)",
    fill = "Consistency outcome"
  )

print(p_motif_state_problem_by_qc)

ggsave(
  filename = file.path(qc_output_folder, "motif_state_consistency_by_pair_QC_status.svg"),
  plot = p_motif_state_problem_by_qc,
  width = 10,
  height = 6,
  limitsize = FALSE
)


library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

alpha <- 0.05

# ============================================================
# Classify row-level significance direction
# ============================================================

stratified_eval <- stratified_test_with_jaccard_qc %>%
  mutate(
    tested = test_status == "tested" &
      !is.na(observed_S_statistic) &
      !is.na(p_value_S_bigger) &
      !is.na(p_value_S_smaller),
    
    sig_positive = tested &
      observed_S_statistic > 0 &
      p_value_S_bigger < alpha,
    
    sig_negative = tested &
      observed_S_statistic < 0 &
      p_value_S_smaller < alpha,
    
    sig_any = sig_positive | sig_negative,
    
    sig_direction = case_when(
      sig_positive ~ "positive",
      sig_negative ~ "negative",
      tested & !sig_any ~ "not_significant",
      TRUE ~ "not_tested"
    ),
    
    pair_QC_reason_chr = as.character(pair_QC_reason),
    
    pair_QC_flag = case_when(
      is.na(pair_QC_reason_chr) ~ NA,
      pair_QC_reason_chr == "Pass" ~ FALSE,
      TRUE ~ TRUE
    )
  )

# ============================================================
# Conflict is defined at protein + motif + chromatin_state level
# ============================================================

conflict_by_motif_state <- stratified_eval %>%
  group_by(protein, motif, chromatin_state) %>%
  summarise(
    n_tested_pairs = sum(tested, na.rm = TRUE),
    n_sig_positive = sum(sig_positive, na.rm = TRUE),
    n_sig_negative = sum(sig_negative, na.rm = TRUE),
    n_sig_total = n_sig_positive + n_sig_negative,
    
    conflicting_significant_direction =
      n_sig_positive > 0 & n_sig_negative > 0,
    
    .groups = "drop"
  )

stratified_eval_conflict <- stratified_eval %>%
  left_join(
    conflict_by_motif_state,
    by = c("protein", "motif", "chromatin_state")
  ) %>%
  mutate(
    conflict_group = case_when(
      conflicting_significant_direction ~ "In conflicting motif-state",
      TRUE ~ "Not conflicting"
    )
  )

# ============================================================
# Pair-state level comparison of QC flags
# ============================================================

qc_flag_comparison_pair_state <- stratified_eval_conflict %>%
  filter(tested) %>%
  select(
    conflict_group,
    pair_low_usable_count,
    pair_high_no_motif,
    pair_low_usable_percentage,
    pair_high_no_cg,
    pair_QC_flag
  ) %>%
  pivot_longer(
    cols = c(
      pair_low_usable_count,
      pair_high_no_motif,
      pair_low_usable_percentage,
      pair_high_no_cg,
      pair_QC_flag
    ),
    names_to = "qc_metric",
    values_to = "qc_flag"
  ) %>%
  filter(!is.na(qc_flag)) %>%
  group_by(conflict_group, qc_metric) %>%
  summarise(
    n_rows = n(),
    n_flagged = sum(qc_flag, na.rm = TRUE),
    percent_flagged = round(100 * n_flagged / n_rows, 2),
    .groups = "drop"
  )

print(qc_flag_comparison_pair_state)

p_qc_flags_conflict_pair_state <- ggplot(
  qc_flag_comparison_pair_state,
  aes(x = qc_metric, y = percent_flagged, fill = conflict_group)
) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_text(
    aes(label = paste0(percent_flagged, "%")),
    position = position_dodge(width = 0.8),
    vjust = -0.25,
    size = 3.5
  ) +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = "QC flags in conflicting vs non-conflicting pair-state rows",
    subtitle = paste0("Only tested rows included; alpha = ", alpha),
    x = "QC metric",
    y = "Flagged rows (%)",
    fill = "Conflict status"
  )

print(p_qc_flags_conflict_pair_state)

# ============================================================
# Unique pair-level version
# One row per protein + motif + experiment pair
# ============================================================

pair_conflict_status <- stratified_eval_conflict %>%
  group_by(
    protein,
    motif,
    biosample1,
    experiment_id1,
    biosample2,
    experiment_id2
  ) %>%
  summarise(
    n_tested_states = sum(tested, na.rm = TRUE),
    n_conflicting_states = sum(conflicting_significant_direction & tested, na.rm = TRUE),
    
    ever_in_conflict = n_conflicting_states > 0,
    
    conflict_group = ifelse(
      ever_in_conflict,
      "Pair ever in conflicting motif-state",
      "Pair never in conflicting motif-state"
    ),
    
    jaccard_index = first(jaccard_index),
    pair_low_usable_count = first(pair_low_usable_count),
    pair_high_no_motif = first(pair_high_no_motif),
    pair_low_usable_percentage = first(pair_low_usable_percentage),
    pair_high_no_cg = first(pair_high_no_cg),
    pair_QC_flag = first(pair_QC_flag),
    pair_QC_reason_chr = first(pair_QC_reason_chr),
    
    usable_count1 = first(usable_count1),
    usable_percentage1 = first(usable_percentage1),
    no_motif_percentage1 = first(no_motif_percentage1),
    no_cg_percentage1 = first(no_cg_percentage1),
    
    usable_count2 = first(usable_count2),
    usable_percentage2 = first(usable_percentage2),
    no_motif_percentage2 = first(no_motif_percentage2),
    no_cg_percentage2 = first(no_cg_percentage2),
    
    .groups = "drop"
  ) %>%
  filter(n_tested_states > 0)


qc_flag_comparison_unique_pair <- pair_conflict_status %>%
  select(
    conflict_group,
    pair_low_usable_count,
    pair_high_no_motif,
    pair_low_usable_percentage,
    pair_high_no_cg,
    pair_QC_flag
  ) %>%
  pivot_longer(
    cols = c(
      pair_low_usable_count,
      pair_high_no_motif,
      pair_low_usable_percentage,
      pair_high_no_cg,
      pair_QC_flag
    ),
    names_to = "qc_metric",
    values_to = "qc_flag"
  ) %>%
  filter(!is.na(qc_flag)) %>%
  group_by(conflict_group, qc_metric) %>%
  summarise(
    n_pairs = n(),
    n_flagged = sum(qc_flag, na.rm = TRUE),
    percent_flagged = round(100 * n_flagged / n_pairs, 2),
    .groups = "drop"
  )

print(qc_flag_comparison_unique_pair)

p_qc_flags_conflict_unique_pair <- ggplot(
  qc_flag_comparison_unique_pair,
  aes(x = qc_metric, y = percent_flagged, fill = conflict_group)
) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_text(
    aes(label = paste0(percent_flagged, "%")),
    position = position_dodge(width = 0.8),
    vjust = -0.25,
    size = 3.5
  ) +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = "QC flags in pairs ever involved in conflict vs never involved",
    subtitle = "Unique experiment-pair level",
    x = "QC metric",
    y = "Flagged pairs (%)",
    fill = "Conflict status"
  )

print(p_qc_flags_conflict_unique_pair)

#7. Plot exact QC reasons

#This shows whether conflicts are enriched in specific QC failure reas
#''


qc_reason_comparison <- pair_conflict_status %>%
  mutate(
    pair_QC_reason_chr = ifelse(
      is.na(pair_QC_reason_chr),
      "No QC match",
      pair_QC_reason_chr
    )
  ) %>%
  count(conflict_group, pair_QC_reason_chr, name = "n_pairs") %>%
  group_by(conflict_group) %>%
  mutate(
    total_pairs = sum(n_pairs),
    percent = round(100 * n_pairs / total_pairs, 2)
  ) %>%
  ungroup()

print(qc_reason_comparison)

p_qc_reason_conflict <- ggplot(
  qc_reason_comparison,
  aes(x = pair_QC_reason_chr, y = percent, fill = conflict_group)
) +
  geom_col(position = position_dodge(width = 0.8)) +
  coord_flip() +
  theme_bw(base_size = 14) +
  labs(
    title = "QC reasons in conflicting vs non-conflicting pairs",
    subtitle = "Unique experiment-pair level",
    x = "Pair QC reason",
    y = "Pairs (%)",
    fill = "Conflict status"
  )

print(p_qc_reason_conflict)

numeric_qc_comparison <- pair_conflict_status %>%
  mutate(
    low_usable_count = pmin(usable_count1, usable_count2, na.rm = TRUE),
    low_usable_percentage = pmin(usable_percentage1, usable_percentage2, na.rm = TRUE),
    high_no_motif_percentage = pmax(no_motif_percentage1, no_motif_percentage2, na.rm = TRUE),
    high_no_cg_percentage = pmax(no_cg_percentage1, no_cg_percentage2, na.rm = TRUE)
  ) %>%
  select(
    conflict_group,
    jaccard_index,
    low_usable_count,
    low_usable_percentage,
    high_no_motif_percentage,
    high_no_cg_percentage
  ) %>%
  pivot_longer(
    cols = c(
      jaccard_index,
      low_usable_count,
      low_usable_percentage,
      high_no_motif_percentage,
      high_no_cg_percentage
    ),
    names_to = "qc_metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

p_numeric_qc_conflict <- ggplot(
  numeric_qc_comparison,
  aes(x = conflict_group, y = value, fill = conflict_group)
) +
  geom_boxplot(outlier.alpha = 0.25) +
  facet_wrap(~ qc_metric, scales = "free_y") +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Numeric QC metrics in conflicting vs non-conflicting pairs",
    subtitle = "Unique experiment-pair level",
    x = "Conflict status",
    y = "QC metric value"
  )

print(p_numeric_qc_conflict)

###########################################################

protein_motif_pairs <- stratified_test_with_jaccard_qc %>%
  distinct(protein, motif) %>%
  arrange(protein, motif)


cat("Unique proteins:", n_distinct(stratified_test_with_jaccard_qc$protein), "\n")
cat("Unique protein-motif pairs:", nrow(protein_motif_pairs), "\n")