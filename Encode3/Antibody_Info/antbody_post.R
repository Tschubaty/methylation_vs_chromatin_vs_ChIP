#################################################################
##  ENCODE Metadata Aggregator
##
##  Input:  final_summary.rds from analysis_Chromatin_State/
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
analysis_dir <- file.path(base_dir, "analysis_Chromatin_State")
antibody_dir <- this.dir
output_dir_plots <- file.path(antibody_dir, "plots")

dir.create(output_dir_plots, recursive = TRUE, showWarnings = FALSE)


# Load final chromatin summary
final_summary_file <- file.path(analysis_dir, "final_summary2.rds")

if (!file.exists(final_summary_file)) {
  stop("final_summary.rds not found at: ", final_summary_file)
}


final_summary <- readRDS(final_summary_file)

cat("Loaded final_summary from:\n", final_summary_file, "\n")
cat("Rows:", nrow(final_summary), "\n")
cat("Columns:", paste(colnames(final_summary), collapse = ", "), "\n")

# final_summary <- final_summary %>%
#   tidyr::separate(
#     protein_hit,
#     into = c("Protein", "File_Accession", "Motif"),
#     sep = "_",
#     remove = FALSE,
#     extra = "merge"
#   )

# load antibody mapping
antibody_df <- read_csv("ENCODE_antibody_complete_mapping.csv")

antibody_min <- antibody_df %>%
  group_by(File_Accession) %>%
  summarise(
    Experiment_Accession = paste(unique(Experiment_Accession), collapse = "; "),
    Antibody_Accession = paste(unique(Antibody_Accession), collapse = "; "),
    .groups = "drop"
  )

final_summary_merged <- final_summary %>%
  left_join(
    antibody_min,
    by = c("ExperimentID" = "File_Accession")
  )

group.colors <- c(
  HepG2    = "#F8766D",
  # A warm reddish color
  K562     = "#00BFC4",
  # A cool cyan color
  GM12878  = "#A3A500",
  # A yellow-green color, good contrast with the others
  A549     = "#E76BF3"   # A vibrant purple color
)

# colors exactly as you want
special_colors <- c(
  "No motif in peak" = "#C2185B",
  "No CG in motif" = "#E91E63",
  "No State assignment for CG" = "#F48FB1",
  "Usable_1_18" = "#1E90FF"
)

# Define chromatin state names and colors
Chromatin_State_colors <- c(
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

Chromatin_State_labels <- c(
  "1_TssA", "2_TssFlnk", "3_TssFlnkU", "4_TssFlnkD",
  "5_Tx", "6_TxWk", "7_EnhG1", "8_EnhG2",
  "9_EnhA1", "10_EnhA2", "11_EnhWk", "12_ZNF_Rpts",
  "13_Het", "14_TssBiv", "15_EnhBiv", "16_ReprPC",
  "17_ReprPCWk", "18_Quies", 
  "No motif in peak", "No CG in motif", "No State assignment for CG","Usable_1_18"
)


names(Chromatin_State_colors) <- c(
  as.character(1:18),
  "No motif in peak",
  "No CG in motif",
  "No State assignment for CG",
  "Usable_1_18"
)

#################################
# create subfolder once before loop
protein_plot_dir <- file.path(this.dir, "protein_barplots")
dir.create(protein_plot_dir, recursive = TRUE, showWarnings = FALSE)

proteins <- unique(final_summary_merged$Protein)

#protein_to_plot <- "ARNT"
if(FALSE){# plot toggel
for (protein_to_plot in proteins) {
  
  plot_order <- c("Chromatin states 1-18", "Special categories")
  
  protein_bar_data <- final_summary_merged %>%
    mutate(
      Chromatin_State = as.character(Chromatin_State),
      Group = paste(
        paste0("Biosample: ", Biosample),
        paste0("File ID: ", Experiment_Accession),
        paste0("Antibody: ", ifelse(is.na(Antibody_Accession), "N/A", Antibody_Accession)),
        sep = "\n"
      ),
      PlotCategory = case_when(
        Chromatin_State %in% as.character(1:18) ~ "Chromatin states 1-18",
        Chromatin_State %in% c(
          "No motif in peak",
          "No CG in motif",
          "No State assignment for CG"
        ) ~ "Special categories",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      Protein == protein_to_plot,
      !is.na(PlotCategory)
    ) %>%
    mutate(
      PlotCategory = factor(
        PlotCategory,
        levels = c("Chromatin states 1-18", "Special categories")
      ),
      Chromatin_State = factor(
        Chromatin_State,
        levels = c(
          as.character(1:18),
          "No motif in peak",
          "No CG in motif",
          "No State assignment for CG"
        )
      ),
      XGroup = paste(PlotCategory, Group, sep = "___")
    )
  
  if (nrow(protein_bar_data) == 0) next
  
  x_levels <- protein_bar_data %>%
    distinct(PlotCategory, Group, XGroup) %>%
    arrange(PlotCategory, Group) %>%
    pull(XGroup)
  
  protein_bar_data <- protein_bar_data %>%
    mutate(XGroup = factor(XGroup, levels = x_levels))
  
  x_labels_df <- protein_bar_data %>%
    distinct(XGroup, PlotCategory, Group) %>%
    arrange(PlotCategory, Group) %>%
    mutate(Label = paste0(as.character(PlotCategory), "\n", Group))
  
  x_labels <- x_labels_df$Label
  names(x_labels) <- x_labels_df$XGroup
  
  p <- ggplot(
    protein_bar_data,
    aes(x = XGroup, y = Count, fill = Chromatin_State)
  ) +
    geom_col(color = "black", width = 0.8) +
    scale_fill_manual(
      values = Chromatin_State_colors[c(
        as.character(1:18),
        "No motif in peak",
        "No CG in motif",
        "No State assignment for CG"
      )],
      drop = FALSE
    ) +
    scale_x_discrete(labels = x_labels) +
    labs(
      title = paste("Absolute counts for protein:", protein_to_plot),
      x = NULL,
      y = "Count"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "right"
    )
  
  print(p)
  
  ggsave(
    filename = file.path(
      protein_plot_dir,
      paste0("barplot_", protein_to_plot, "_absolute_counts.png")
    ),
    plot = p,
    width = 20,
    height = 10,
    limitsize = FALSE
  )
  
  cat("Saved plot for protein:", protein_to_plot, "\n")
}
} 

############################# sum tables
# 
# target_antibody_biosample <- ab_exp %>%
#   group_by(Target, Antibody_Accession) %>%
#   summarise(
#     n_biosamples = n_distinct(Biosample),
#     biosamples = paste(sort(unique(Biosample)), collapse = "; "),
#     n_experiments = n_distinct(Experiment_Accession),
#     .groups = "drop"
#   ) %>%
#   arrange(desc(n_biosamples), desc(n_experiments))
# 
# print(target_antibody_biosample, n = Inf)
# 
# target_antibody_biosample %>%
#   summarise(
#     total_target_antibody_pairs = n(),
#     pairs_with_2plus_celllines = sum(n_biosamples >= 2),
#     pairs_with_3plus_celllines = sum(n_biosamples >= 3),
#     pairs_with_4plus_celllines = sum(n_biosamples >= 4)
#   )


# PCA 
library(dplyr)
library(tidyr)

pca_input_all <- final_summary_merged %>%
  mutate(Chromatin_State = as.character(Chromatin_State)) %>%
  group_by(ExperimentID, Protein, Biosample, Antibody_Accession) %>%
  mutate(total_peaks = sum(Count, na.rm = TRUE)) %>%
  ungroup() %>%
  select(
    ExperimentID, Protein, Biosample, Antibody_Accession,
    Chromatin_State, Percentage, total_peaks
  ) %>%
  distinct() %>%
  pivot_wider(
    names_from = Chromatin_State,
    values_from = Percentage,
    values_fill = 0
  ) %>%
  mutate(
    log_total_peaks = log10(total_peaks + 1)
  )

library(ggplot2)

mat <- pca_input_all %>%
  select(-ExperimentID, -Protein, -Biosample, -Antibody_Accession)

pca <- prcomp(mat, scale. = TRUE)

pca_df <- as.data.frame(pca$x) %>%
  bind_cols(
    pca_input_all %>%
      select(ExperimentID, Protein, Biosample, Antibody_Accession)
  )

ggplot(pca_df, aes(PC1, PC2, color = Biosample)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "Global PCA of all ChIP-seq experiments",
    subtitle = "Color = cell line",
    x = "PC1",
    y = "PC2"
  )


# PCA with plot of specific protein lables 

library(ggrepel)

protein_to_label <- "ZBTB33"   # change this

p <- ggplot(pca_df, aes(PC1, PC2, color = Biosample)) +
  geom_point(size = 2, alpha = 0.45) +
  geom_point(
    data = pca_df %>% filter(Protein == protein_to_label),
    size = 3.5,
    alpha = 1
  ) +
  geom_text_repel(
    data = pca_df %>% filter(Protein == protein_to_label),
    aes(label = Antibody_Accession),
    size = 3,
    color = "black",
    max.overlaps = Inf
  ) +
  scale_color_manual(values = group.colors) +
  theme_minimal() +
  labs(
    title = paste("Global PCA highlighted target:", protein_to_label),
    subtitle = "Highlighted points labeled by antibody accession",
    x = "PC1",
    y = "PC2",
    color = "Biosample"
  )

print(p)