#################################################################
##  QC_motif_CpG_antibody_jaccard_analysis_v1.R
##
##  Purpose:
##
##  This script performs quality-control analysis for the
##  methylation-vs-ChIP project.
##
##  It combines:
##    1. Pairwise Jaccard overlap between ChIP-seq peak sets
##       to evaluate experiment comparability by biosample and antibody.
##
##    2. Motif/CpG/chromatin-state QC information from final_summary2
##       to evaluate whether each ChIP-seq experiment is usable for
##       motif-level CpG methylation sensitivity analysis.
##
##  Main QC questions:
##    - How many peaks have no motif in peak?
##    - How many motif instances have no CpG in motif?
##    - How many motif-CpG sites lack chromatin-state assignment?
##    - How many usable motif-CpG sites remain in chromatin states 1-18?
##    - Do experiments cluster or differ by biosample, antibody, or both?
##    - Which experiments should be excluded or flagged before downstream
##      methylation sensitivity testing?
##
##  Input:
##    Everything is in hg38.
##
##    1. Pairwise Jaccard table:
##       Encode3/jaccard_analysis/Jaccard_grouped_analysis/
##       ALL_pairwise_jaccard_values_with_groups.csv
##
##    2. Chromatin-state / motif-CpG QC summary:
##       Encode3/analysis_chromatin_state/final_summary2.rds
##
##  Output:
##    Encode3/QC_motif_CpG_antibody_jaccard/
##
##  Version:
##    v1
##
##  Date:
##    2026-05-26
##
##  Author:
##    Daniel Batyrev
#################################################################

################################ Clear environment ##############################

rm(list = ls())

################################ Settings #######################################

#Clear R working environment

cluster <- FALSE
if (cluster) {
  this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/Encode3/"
  picuture_file_extension <- "pdf"
} else{
  this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  picuture_file_extension <- "svg"
}

setwd(this.dir)

################################ Libraries ######################################

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)


# Define custom colors
group.colors <- c(
  HepG2    = "#F8766D",  # Warm reddish
  K562     = "#00BFC4",  # Cool cyan
  GM12878  = "#A3A500",  # Yellow-green
  A549     = "#E76BF3"   # Vibrant purple
)

# Define chromatin state names and colors
chromatin_state_colors <- c(
  "#FF0000",  # 1_TssA       -> Red
  "#FF4500",  # 2_TssFlnk    -> Orange-Red
  "#FF9900",  # 3_TssFlnkU   -> Orange
  "#FFCC00",  # 4_TssFlnkD   -> Yellow-Orange
  "#00CC00",  # 5_Tx         -> Green
  "#006400",  # 6_TxWk       -> Dark Green
  "#FFD700",  # 7_EnhG1      -> Gold
  "#FFD700",  # 8_EnhG2      -> Gold (same as EnhG1)
  "#FFFF00",  # 9_EnhA1      -> Yellow
  "#FFDD00",  # 10_EnhA2     -> Yellow-Orange
  "#FFEA73",  # 11_EnhWk     -> Light Yellow
  "#9370DB",  # 12_ZNF_Rpts  -> Purple
  "#C0C0C0",  # 13_Het       -> Light Gray
  "#FF4500",  # 14_TssBiv    -> Orange-Red (similar to TssFlnk)
  "#FFDD00",  # 15_EnhBiv    -> Yellow-Orange (similar to EnhA2)
  "#808080",  # 16_ReprPC    -> Dark Gray
  "#A9A9A9",  # 17_ReprPCWk  -> Light Gray
  "#000000",  # 18_Quies     -> Black
  "#C2185B",  # No motif in peak -> Pink
  "#E91E63",  # No CG in motif -> Dodger Blue
  "#F48FB1",   # No State assignment for CG -> Lime Green,
  "#1E90FF"  # "Usable_1_18"
  
)

chromatin_state_labels <- c(
  "1_TssA", "2_TssFlnk", "3_TssFlnkU", "4_TssFlnkD",
  "5_Tx", "6_TxWk", "7_EnhG1", "8_EnhG2",
  "9_EnhA1", "10_EnhA2", "11_EnhWk", "12_ZNF_Rpts",
  "13_Het", "14_TssBiv", "15_EnhBiv", "16_ReprPC",
  "17_ReprPCWk", "18_Quies", 
  "No motif in peak", "No CG in motif", "No State assignment for CG","Usable_1_18"
)


################################ Input files ####################################

jaccard_pairs_file <- file.path(
  this.dir,
  "jaccard_analysis",
  "Jaccard_grouped_analysis",
  "ALL_pairwise_jaccard_values_with_groups.csv"
)

final_summary2_file <- file.path(
  this.dir,
  "analysis_chromatin_state",
  "final_summary2.rds"
)

################################ Output folder ##################################

qc_output_dir <- file.path(
  this.dir,
  "QC_motif_CpG_antibody_jaccard"
)

if (!dir.exists(qc_output_dir)) {
  dir.create(qc_output_dir, recursive = TRUE)
}

################################ Check input files ###############################

if (!file.exists(jaccard_pairs_file)) {
  stop(paste("Jaccard pairwise file not found:", jaccard_pairs_file))
}

if (!file.exists(final_summary2_file)) {
  stop(paste("final_summary2 RDS file not found:", final_summary2_file))
}

################################ Load data ######################################

print("Loading pairwise Jaccard table...")

all_jaccard_pairs <- readr::read_csv(
  jaccard_pairs_file,
  show_col_types = FALSE,
  col_types = readr::cols(
    .default = readr::col_character(),
    jaccard_index = readr::col_double(),
    same_biosample = readr::col_logical(),
    same_antibody = readr::col_logical(),
    antibody_known_1 = readr::col_logical(),
    antibody_known_2 = readr::col_logical(),
    both_antibodies_known = readr::col_logical(),
    same_antibody_strict = readr::col_logical(),
    different_antibody_strict = readr::col_logical()
  )
)

print("Loading final_summary2...")

final_summary2 <- readRDS(final_summary2_file)

################################ Basic checks ###################################

print("Structure of all_jaccard_pairs:")
print(str(all_jaccard_pairs))

print("Structure of final_summary2:")
print(str(final_summary2))

print("Chromatin_State levels in final_summary2:")
print(levels(final_summary2$Chromatin_State))

################################ Clean variables ################################
qc_percentage_states <- c(
  "No motif in peak",
  "No CG in motif",
  "No State assignment for CG"
)

usable_state <- "Usable_1_18"

final_summary2_qc_wide <- final_summary2 %>%
  dplyr::filter(
    Chromatin_State %in% c(qc_percentage_states, usable_state)
  ) %>%
  dplyr::select(
    Protein,
    ExperimentID,
    Biosample,
    Chromatin_State,
    Count,
    Percentage
  ) %>%
  tidyr::pivot_wider(
    names_from = Chromatin_State,
    values_from = c(Count, Percentage),
    names_glue = "{Chromatin_State}_{.value}"
  ) %>%
  dplyr::select(
    Protein,
    ExperimentID,
    Biosample,
    `No motif in peak_Percentage`,
    `No CG in motif_Percentage`,
    `No State assignment for CG_Percentage`,
    Usable_1_18_Count,
    Usable_1_18_Percentage
  ) %>%
  dplyr::rename_with(~ stringr::str_replace_all(.x, " ", "_")) %>%
  dplyr::rename_with(~ stringr::str_replace_all(.x, "-", "_"))

##########
################################ Compact Jaccard pair table #####################

all_jaccard_pairs_compact <- all_jaccard_pairs %>%
  dplyr::select(
    protein,
    motif,
    jaccard_index,
    
    biosample_1,
    experiment_id_1,
    antibody_label_1,
    
    biosample_2,
    experiment_id_2,
    antibody_label_2,
    
    same_biosample,
    same_antibody,

    biosample_group,
    antibody_group,
    combined_group,
    biosample_pair,
  )

################################ Add motif/CpG QC to Jaccard pairs ##############
################################ Add motif/CpG QC to Jaccard pairs ##############

# Make experiment-level QC table for joining
experiment_qc_lookup <- final_summary2_qc_wide %>%
  dplyr::mutate(
    Protein = as.character(Protein),
    ExperimentID = as.character(ExperimentID),
    Biosample = as.character(Biosample)
  ) %>%
  dplyr::select(
    Protein,
    ExperimentID,
    Biosample,
    No_motif_in_peak_Percentage,
    No_CG_in_motif_Percentage,
    No_State_assignment_for_CG_Percentage,
    Usable_1_18_Count,
    Usable_1_18_Percentage
  ) %>%
  dplyr::distinct(
    Protein,
    ExperimentID,
    Biosample,
    .keep_all = TRUE
  )

# Add QC info for both experiments in each Jaccard pair
all_jaccard_pairs_with_qc <- all_jaccard_pairs_compact %>%
  
  # Add QC info for experiment 1
  dplyr::left_join(
    experiment_qc_lookup,
    by = c(
      "protein" = "Protein",
      "experiment_id_1" = "ExperimentID",
      "biosample_1" = "Biosample"
    )
  ) %>%
  dplyr::rename(
    No_motif_in_peak_Percentage_1 = No_motif_in_peak_Percentage,
    No_CG_in_motif_Percentage_1 = No_CG_in_motif_Percentage,
    No_State_assignment_for_CG_Percentage_1 = No_State_assignment_for_CG_Percentage,
    Usable_1_18_Count_1 = Usable_1_18_Count,
    Usable_1_18_Percentage_1 = Usable_1_18_Percentage
  ) %>%
  
  # Add QC info for experiment 2
  dplyr::left_join(
    experiment_qc_lookup,
    by = c(
      "protein" = "Protein",
      "experiment_id_2" = "ExperimentID",
      "biosample_2" = "Biosample"
    )
  ) %>%
  dplyr::rename(
    No_motif_in_peak_Percentage_2 = No_motif_in_peak_Percentage,
    No_CG_in_motif_Percentage_2 = No_CG_in_motif_Percentage,
    No_State_assignment_for_CG_Percentage_2 = No_State_assignment_for_CG_Percentage,
    Usable_1_18_Count_2 = Usable_1_18_Count,
    Usable_1_18_Percentage_2 = Usable_1_18_Percentage
  )

####

library(viridis)

p_usable_vs_usable_jaccard <- ggplot(
  all_jaccard_pairs_with_qc,
  aes(
    x = Usable_1_18_Percentage_1,
    y = Usable_1_18_Percentage_2,
    color = jaccard_index
  )
) +
  geom_point(alpha = 0.75, size = 2.5) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title = "Relationship between usable motif-CpG percentage and Jaccard overlap",
    subtitle = "Each point is one pairwise ChIP-seq experiment comparison",
    x = "Usable motif-CpG percentage: experiment 1",
    y = "Usable motif-CpG percentage: experiment 2",
    color = "Jaccard index"
  ) +
  theme_minimal(base_size = 16)

print(p_usable_vs_usable_jaccard)

ggsave(
  filename = file.path(qc_output_dir, "scatter_usable_percentage_exp1_vs_exp2_colored_by_jaccard.png"),
  plot = p_usable_vs_usable_jaccard,
  width = 10,
  height = 8,
  dpi = 180,
  limitsize = FALSE
)

################################ QC histograms: motif/CpG usability #############


library(patchwork)

# Make binned usable-count variable:
# 0-500, 500-1000, ..., 9500-10000, >10000
final_summary2_qc_wide_plot <- final_summary2_qc_wide %>%
  dplyr::mutate(
    Usable_1_18_Count_bin = cut(
      Usable_1_18_Count,
      breaks = c(seq(0, 10000, by = 500), Inf),
      include.lowest = TRUE,
      right = FALSE,
      labels = c(
        paste0(seq(0, 9500, by = 500), "-", seq(500, 10000, by = 500)),
        ">10000"
      )
    )
  )

# 1. Histogram: No motif in peak percentage
p_no_motif <- ggplot(
  final_summary2_qc_wide_plot,
  aes(x = No_motif_in_peak_Percentage)
) +
  geom_histogram(
    bins = 40,
    color = "black",
    fill = "#CC79A7"
  ) +
  labs(
    title = "No motif in peak",
    x = "Percentage",
    y = "Number of experiments"
  ) +
  theme_minimal(base_size = 16)

# 2. Histogram: No CG in motif percentage
p_no_cg <- ggplot(
  final_summary2_qc_wide_plot,
  aes(x = No_CG_in_motif_Percentage)
) +
  geom_histogram(
    bins = 40,
    color = "black",
    fill = "#E69F00"
  ) +
  labs(
    title = "No CpG in motif",
    x = "Percentage",
    y = "Number of experiments"
  ) +
  theme_minimal(base_size = 16)

# 3. Histogram: Usable percentage
p_usable_percentage <- ggplot(
  final_summary2_qc_wide_plot,
  aes(x = Usable_1_18_Percentage)
) +
  geom_histogram(
    bins = 40,
    color = "black",
    fill = "#009E73"
  ) +
  labs(
    title = "Usable motif-CpG sites",
    x = "Percentage",
    y = "Number of experiments"
  ) +
  theme_minimal(base_size = 16)

# 4. Barplot: Usable count, capped at >10000
p_usable_count <- ggplot(
  final_summary2_qc_wide_plot,
  aes(x = Usable_1_18_Count_bin)
) +
  geom_bar(
    color = "black",
    fill = "#0072B2"
  ) +
  labs(
    title = "Number of usable motif-CpG sites",
    subtitle = "Experiments with >10,000 usable sites are grouped into one bin",
    x = "Usable_1_18 count bin",
    y = "Number of experiments"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Combine into 2x2 grid
p_qc_histograms_2x2 <- (
  p_no_motif + p_no_cg
) / (
  p_usable_percentage + p_usable_count
) +
  plot_annotation(
    title = "Motif/CpG usability QC across ChIP-seq experiments"
  )

# Preview in RStudio
print(p_qc_histograms_2x2)

# Save
ggsave(
  filename = file.path(qc_output_dir, "QC_histograms_motif_CpG_usability_2x2.png"),
  plot = p_qc_histograms_2x2,
  width = 16,
  height = 11,
  dpi = 180,
  limitsize = FALSE

  
)

################################ QC histograms with threshold highlighting ########

library(patchwork)

# Thresholds for initial QC flagging
cutoff_no_motif <- 50
cutoff_no_cg <- 75
cutoff_usable_percentage <- 10
cutoff_usable_count <- 500

################################ Dynamic QC labels and colors ###################

qc_label_low_usable_count <- paste0(
  "Low usable count (<", cutoff_usable_count, ")"
)

qc_label_high_no_motif <- paste0(
  "High no motif in peak (>", cutoff_no_motif, "%)"
)

qc_label_low_usable_percentage <- paste0(
  "Low usable percentage (<", cutoff_usable_percentage, "%)"
)

qc_label_high_no_cg <- paste0(
  "High no CpG in motif (>", cutoff_no_cg, "%)"
)

qc_label_ok <- "OK"

qc_reason_levels <- c(
  qc_label_low_usable_count,
  qc_label_high_no_motif,
  qc_label_low_usable_percentage,
  qc_label_high_no_cg,
  qc_label_ok
)

qc_reason_colors <- c(
  "#004C73",
  "#8B1C62",
  "#00664F",
  "#9A6700",
  "#80CDC1"
)

names(qc_reason_colors) <- qc_reason_levels

# Make binned usable-count variable:
# 0-500, 500-1000, ..., 9500-10000, >10000
final_summary2_qc_wide_plot <- final_summary2_qc_wide %>%
  dplyr::mutate(
    Usable_1_18_Count_bin = cut(
      Usable_1_18_Count,
      breaks = c(seq(0, 10000, by = 500), Inf),
      include.lowest = TRUE,
      right = FALSE,
      labels = c(
        paste0(seq(0, 9500, by = 500), "-", seq(500, 10000, by = 500)),
        ">10000"
      )
    ),
    Usable_1_18_Count_bin_problematic = Usable_1_18_Count < cutoff_usable_count
  )

# 1. Histogram: No motif in peak percentage
p_no_motif <- ggplot(
  final_summary2_qc_wide_plot,
  aes(
    x = No_motif_in_peak_Percentage,
    fill = after_stat(x > cutoff_no_motif)
  )
) +
  geom_histogram(
    bins = 40,
    color = "black"
  ) +
  geom_vline(
    xintercept = cutoff_no_motif,
    color = "#8B1C62",
    linewidth = 0.6
  ) +
  scale_fill_manual(
    values = c(
      "FALSE" = "#E8B7D4",
      "TRUE" = "#8B1C62"
    ),
    labels = c(
      "OK",
      "Problematic"
    ),
    name = "QC flag"
  ) +
  labs(
    title = "No motif in peak",
    subtitle = paste0("Problematic: > ", cutoff_no_motif, "%"),
    x = "Percentage",
    y = "Number of experiments"
  ) +
  theme_minimal(base_size = 16)

# 2. Histogram: No CpG in motif percentage
p_no_cg <- ggplot(
  final_summary2_qc_wide_plot,
  aes(
    x = No_CG_in_motif_Percentage,
    fill = after_stat(x > cutoff_no_cg)
  )
) +
  geom_histogram(
    bins = 40,
    color = "black"
  ) +
  geom_vline(
    xintercept = cutoff_no_cg,
    color = "#9A6700",
    linewidth = 0.6
  ) +
  scale_fill_manual(
    values = c(
      "FALSE" = "#F5C66A",
      "TRUE" = "#9A6700"
    ),
    labels = c(
      "OK",
      "Problematic"
    ),
    name = "QC flag"
  ) +
  labs(
    title = "No CpG in motif",
    subtitle = paste0("Warning: > ", cutoff_no_cg, "%"),
    x = "Percentage",
    y = "Number of experiments"
  ) +
  theme_minimal(base_size = 16)

# 3. Histogram: Usable percentage
p_usable_percentage <- ggplot(
  final_summary2_qc_wide_plot,
  aes(
    x = Usable_1_18_Percentage,
    fill = after_stat(x < cutoff_usable_percentage)
  )
) +
  geom_histogram(
    bins = 40,
    color = "black"
  ) +
  geom_vline(
    xintercept = cutoff_usable_percentage,
    color = "#00664F",
    linewidth = 0.6
  ) +
  scale_fill_manual(
    values = c(
      "FALSE" = "#80CDC1",
      "TRUE" = "#00664F"
    ),
    labels = c(
      "OK",
      "Problematic"
    ),
    name = "QC flag"
  ) +
  labs(
    title = "Usable motif-CpG sites",
    subtitle = paste0("Problematic: < ", cutoff_usable_percentage, "%"),
    x = "Percentage",
    y = "Number of experiments"
  ) +
  theme_minimal(base_size = 16)

# 4. Barplot: Usable count, capped at >10000
p_usable_count <- ggplot(
  final_summary2_qc_wide_plot,
  aes(
    x = Usable_1_18_Count_bin,
    fill = Usable_1_18_Count_bin_problematic
  )
) +
  geom_bar(
    color = "black"
  ) +
  geom_vline(
    xintercept = 1.5,
    color = "#004C73",
    linewidth = 0.6
  ) +
  scale_fill_manual(
    values = c(
      "FALSE" = "#8ECAE6",
      "TRUE" = "#004C73"
    ),
    labels = c(
      "OK",
      "Problematic"
    ),
    name = "QC flag"
  ) +
  labs(
    title = "Number of usable motif-CpG sites",
    subtitle = paste0(
      "Problematic: < ", cutoff_usable_count,
      "; experiments with >10,000 usable sites grouped into one bin"
    ),
    x = "Usable_1_18 count bin",
    y = "Number of experiments"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Combine into 2x2 grid
p_qc_histograms_2x2 <- (
  p_no_motif + p_no_cg
) / (
  p_usable_percentage + p_usable_count
) +
  plot_annotation(
    title = "Motif/CpG usability QC across ChIP-seq experiments"
  )

# Preview in RStudio
print(p_qc_histograms_2x2)

# Save
ggsave(
  filename = file.path(qc_output_dir, "QC_histograms_motif_CpG_usability_2x2_threshold_highlighted.png"),
  plot = p_qc_histograms_2x2,
  width = 16,
  height = 11,
  dpi = 180,
  limitsize = FALSE
)

################################ Usable percentage split by QC reason ############

final_summary2_qc_wide_plot <- final_summary2_qc_wide_plot %>%
  dplyr::mutate(
    QC_reason_hierarchy = dplyr::case_when(
      Usable_1_18_Count < cutoff_usable_count ~ qc_label_low_usable_count,
      No_motif_in_peak_Percentage > cutoff_no_motif ~ qc_label_high_no_motif,
      Usable_1_18_Percentage < cutoff_usable_percentage ~ qc_label_low_usable_percentage,
      No_CG_in_motif_Percentage > cutoff_no_cg ~ qc_label_high_no_cg,
      TRUE ~ qc_label_ok
    ),
    QC_reason_hierarchy = factor(
      QC_reason_hierarchy,
      levels = qc_reason_levels
    )
  )

# 3. Histogram: Usable percentage, split by hierarchical QC reason
p_usable_percentage <- ggplot(
  final_summary2_qc_wide_plot,
  aes(
    x = Usable_1_18_Percentage,
    fill = Usable_percentage_QC_reason
  )
) +
  geom_histogram(
    bins = 40,
    color = "black",
    position = "stack"
  ) +
  geom_vline(
    xintercept = cutoff_usable_percentage,
    color = "#00664F",
    linewidth = 0.6
  ) +
  scale_fill_manual(
    values = c(
      "Low usable count (<500)" = "#004C73",
      "High no motif in peak (>50%)" = "#8B1C62",
      "High no CpG in motif (>75%)" = "#9A6700",
      "OK" = "#80CDC1"
    ),
    name = "QC reason"
  ) +
  labs(
    title = "Usable motif-CpG sites",
    subtitle = paste0(
      "Bars split by hierarchical QC reason; vertical line = ",
      cutoff_usable_percentage,
      "%"
    ),
    x = "Usable motif-CpG percentage",
    y = "Number of experiments"
  ) +
  theme_minimal(base_size = 16)

################################ Shared hierarchical QC reason ###################

final_summary2_qc_wide_plot <- final_summary2_qc_wide_plot %>%
  dplyr::mutate(
    QC_reason_hierarchy = dplyr::case_when(
      Usable_1_18_Count < cutoff_usable_count ~ qc_label_low_usable_count,
      No_motif_in_peak_Percentage > cutoff_no_motif ~ qc_label_high_no_motif,
      Usable_1_18_Percentage < cutoff_usable_percentage ~ qc_label_low_usable_percentage,
      No_CG_in_motif_Percentage > cutoff_no_cg ~ qc_label_high_no_cg,
      TRUE ~ qc_label_ok
    ),
    QC_reason_hierarchy = factor(
      QC_reason_hierarchy,
      levels = qc_reason_levels
    )
  )

################################ QC histograms with shared QC coloring ###########

# 1. Histogram: No motif in peak percentage
p_no_motif <- ggplot(
  final_summary2_qc_wide_plot,
  aes(
    x = No_motif_in_peak_Percentage,
    fill = QC_reason_hierarchy
  )
) +
  geom_histogram(
    bins = 40,
    color = "black",
    position = "stack"
  ) +
  geom_vline(
    xintercept = cutoff_no_motif,
    color = "#8B1C62",
    linewidth = 0.7
  ) +
  scale_fill_manual(
    values = qc_reason_colors,
    name = "QC reason"
  ) +
  labs(
    title = "No motif in peak",
    subtitle = paste0("Threshold: > ", cutoff_no_motif, "%"),
    x = "Percentage",
    y = "Number of experiments"
  ) +
  theme_minimal(base_size = 16)

# 2. Histogram: No CpG in motif percentage
p_no_cg <- ggplot(
  final_summary2_qc_wide_plot,
  aes(
    x = No_CG_in_motif_Percentage,
    fill = QC_reason_hierarchy
  )
) +
  geom_histogram(
    bins = 40,
    color = "black",
    position = "stack"
  ) +
  geom_vline(
    xintercept = cutoff_no_cg,
    color = "#9A6700",
    linewidth = 0.7
  ) +
  scale_fill_manual(
    values = qc_reason_colors,
    name = "QC reason"
  ) +
  labs(
    title = "No CpG in motif",
    subtitle = paste0("Warning threshold: > ", cutoff_no_cg, "%"),
    x = "Percentage",
    y = "Number of experiments"
  ) +
  theme_minimal(base_size = 16)

# 3. Histogram: Usable percentage
p_usable_percentage <- ggplot(
  final_summary2_qc_wide_plot,
  aes(
    x = Usable_1_18_Percentage,
    fill = QC_reason_hierarchy
  )
) +
  geom_histogram(
    bins = 40,
    color = "black",
    position = "stack"
  ) +
  geom_vline(
    xintercept = cutoff_usable_percentage,
    color = "#00664F",
    linewidth = 0.7
  ) +
  scale_fill_manual(
    values = qc_reason_colors,
    name = "QC reason"
  ) +
  labs(
    title = "Usable motif-CpG sites",
    subtitle = paste0("Threshold: < ", cutoff_usable_percentage, "%"),
    x = "Percentage",
    y = "Number of experiments"
  ) +
  theme_minimal(base_size = 16)

# 4. Barplot: Usable count, capped at >10000
p_usable_count <- ggplot(
  final_summary2_qc_wide_plot,
  aes(
    x = Usable_1_18_Count_bin,
    fill = QC_reason_hierarchy
  )
) +
  geom_bar(
    color = "black",
    position = "stack"
  ) +
  geom_vline(
    xintercept = 1.5,
    color = "#004C73",
    linewidth = 0.7
  ) +
  scale_fill_manual(
    values = qc_reason_colors,
    name = "QC reason"
  ) +
  labs(
    title = "Number of usable motif-CpG sites",
    subtitle = paste0(
      "Threshold: < ", cutoff_usable_count,
      "; experiments with >10,000 usable sites grouped into one bin"
    ),
    x = "Usable_1_18 count bin",
    y = "Number of experiments"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_qc_histograms_2x2 <- (
  p_no_motif + p_no_cg
) / (
  p_usable_percentage + p_usable_count
) +
  plot_annotation(
    title = "Motif/CpG usability QC across ChIP-seq experiments"
  )

print(p_qc_histograms_2x2)

ggsave(
  filename = file.path(qc_output_dir, "QC_histograms_motif_CpG_usability_2x2_shared_QC_coloring.png"),
  plot = p_qc_histograms_2x2,
  width = 16,
  height = 11,
  dpi = 180,
  limitsize = FALSE
)

################################ Top problematic candidates by QC criterion ######
################################ Top problematic candidates by QC criterion ######

################################ Check how many problematic candidates ###########

final_summary2_qc_wide_problematic <- final_summary2_qc_wide %>%
  dplyr::mutate(
    flag_low_usable_count = Usable_1_18_Count < cutoff_usable_count,
    flag_high_no_motif = No_motif_in_peak_Percentage > cutoff_no_motif,
    flag_high_no_cg = No_CG_in_motif_Percentage > cutoff_no_cg,
    flag_low_usable_percentage = Usable_1_18_Percentage < cutoff_usable_percentage
  )

qc_flag_counts <- final_summary2_qc_wide_problematic %>%
  dplyr::summarise(
    n_total_experiments = dplyr::n(),
    n_low_usable_count = sum(flag_low_usable_count, na.rm = TRUE),
    n_high_no_motif = sum(flag_high_no_motif, na.rm = TRUE),
    n_high_no_cg = sum(flag_high_no_cg, na.rm = TRUE),
    n_low_usable_percentage = sum(flag_low_usable_percentage, na.rm = TRUE)
  )

print(qc_flag_counts)
#View(qc_flag_counts)

################################ Top problematic candidates ######################

top_n_problematic <- 25

top_low_usable_count <- final_summary2_qc_wide_problematic %>%
  dplyr::arrange(Usable_1_18_Count) %>%
  dplyr::select(
    Protein,
    ExperimentID,
    Biosample,
    Usable_1_18_Count,
    Usable_1_18_Percentage,
    No_motif_in_peak_Percentage,
    No_CG_in_motif_Percentage,
    No_State_assignment_for_CG_Percentage
  ) %>%
  dplyr::slice_head(n = top_n_problematic)

top_high_no_motif <- final_summary2_qc_wide_problematic %>%
  dplyr::arrange(dplyr::desc(No_motif_in_peak_Percentage)) %>%
  dplyr::select(
    Protein,
    ExperimentID,
    Biosample,
    No_motif_in_peak_Percentage,
    Usable_1_18_Count,
    Usable_1_18_Percentage,
    No_CG_in_motif_Percentage,
    No_State_assignment_for_CG_Percentage
  ) %>%
  dplyr::slice_head(n = top_n_problematic)

top_high_no_cg <- final_summary2_qc_wide_problematic %>%
  dplyr::arrange(dplyr::desc(No_CG_in_motif_Percentage)) %>%
  dplyr::select(
    Protein,
    ExperimentID,
    Biosample,
    No_CG_in_motif_Percentage,
    Usable_1_18_Count,
    Usable_1_18_Percentage,
    No_motif_in_peak_Percentage,
    No_State_assignment_for_CG_Percentage
  ) %>%
  dplyr::slice_head(n = top_n_problematic)

top_low_usable_percentage <- final_summary2_qc_wide_problematic %>%
  dplyr::arrange(Usable_1_18_Percentage) %>%
  dplyr::select(
    Protein,
    ExperimentID,
    Biosample,
    Usable_1_18_Percentage,
    Usable_1_18_Count,
    No_motif_in_peak_Percentage,
    No_CG_in_motif_Percentage,
    No_State_assignment_for_CG_Percentage
  ) %>%
  dplyr::slice_head(n = top_n_problematic)

#View(top_low_usable_count)
#View(top_high_no_motif)
#View(top_high_no_cg)
#View(top_low_usable_percentage)


################################ Jaccard plot colored by motif/CpG QC ############

################################ Add pair-level QC reason to Jaccard table #######

all_jaccard_pairs_with_qc <- all_jaccard_pairs_with_qc %>%
  dplyr::mutate(
    pair_low_usable_count = dplyr::coalesce(
      Usable_1_18_Count_1 < cutoff_usable_count |
        Usable_1_18_Count_2 < cutoff_usable_count,
      FALSE
    ),
    
    pair_high_no_motif = dplyr::coalesce(
      No_motif_in_peak_Percentage_1 > cutoff_no_motif |
        No_motif_in_peak_Percentage_2 > cutoff_no_motif,
      FALSE
    ),
    
    pair_low_usable_percentage = dplyr::coalesce(
      Usable_1_18_Percentage_1 < cutoff_usable_percentage |
        Usable_1_18_Percentage_2 < cutoff_usable_percentage,
      FALSE
    ),
    
    pair_high_no_cg = dplyr::coalesce(
      No_CG_in_motif_Percentage_1 > cutoff_no_cg |
        No_CG_in_motif_Percentage_2 > cutoff_no_cg,
      FALSE
    ),
    
    pair_QC_reason = dplyr::case_when(
      pair_low_usable_count ~ qc_label_low_usable_count,
      pair_high_no_motif ~ qc_label_high_no_motif,
      pair_low_usable_percentage ~ qc_label_low_usable_percentage,
      pair_high_no_cg ~ qc_label_high_no_cg,
      TRUE ~ qc_label_ok
    ),
    
    pair_QC_reason = factor(
      pair_QC_reason,
      levels = qc_reason_levels
    )
  )

p_jaccard_biosample_qc <- ggplot(
  all_jaccard_pairs_with_qc,
  aes(
    x = reorder(biosample_pair, jaccard_index, FUN = median),
    y = jaccard_index
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    fill = "grey85",
    color = "black"
  ) +
  geom_jitter(
    aes(color = pair_QC_reason),
    width = 0.18,
    height = 0,
    alpha = 0.75,
    size = 2.3
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 4,
    fill = "white",
    color = "black"
  ) +
  scale_color_manual(
    values = qc_reason_colors,
    name = "Pair QC reason"
  ) +
  coord_flip() +
  labs(
    title = "Jaccard value distribution by cell-line pair",
    subtitle = "Each point is one pairwise experiment comparison; diamond = mean",
    x = "Cell-line pair",
    y = "Jaccard index"
  ) +
  theme_minimal(base_size = 18)+
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.7)

print(p_jaccard_biosample_qc)

ggsave(
  filename = file.path(
    qc_output_dir,
    "jaccard_distribution_by_cell_line_pair_colored_by_pair_QC_reason.png"
  ),
  plot = p_jaccard_biosample_qc,
  width = 14,
  height = 10,
  dpi = 180,
  limitsize = FALSE
)