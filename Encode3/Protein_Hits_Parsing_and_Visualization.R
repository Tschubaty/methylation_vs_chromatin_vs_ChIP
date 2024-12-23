#################################################################
##  Description: Analyze CpGs covered by proteins across biosamples.
##  Summarizes counts by collapsing different experiments/binding motifs of the same protein.
##
##  Input: Methylation and protein occupancy data (HG38)
##         combined_methylation_data_compact_withpeaks.rds
##
##  Output: Single dodging bar plot (SVG) of CpG counts by summarized proteins across all biosamples.
##          Directory: Encode3/Protein_Hits_Parsing_and_Visualization/plots/
##
##  Version: 22.12.2024
##  Author: Daniel Batyrev (777634015)
#################################################################

# Clear R working environment
rm(list = ls())
cluster <- FALSE
if (cluster) {
  this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/Encode3/"
  picuture_file_extension <- "pdf"
} else {
  this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  picuture_file_extension <- "png"
}

setwd(this.dir)

#################################### Libs ########################################
library(foreach)
library(doParallel)
library(readr)
library(dplyr)
library(rlang)
library(ggplot2)
library(tidyr)
library(stringr)

##################################### INPUT ########################################
# Define input and output directories
output_analysis_dir <- file.path(this.dir, "Protein_Hits_Parsing_and_Visualization")
output_dir_plots <- file.path(output_analysis_dir, "plots")

if (!dir.exists(output_analysis_dir)) dir.create(output_analysis_dir, recursive = TRUE)
if (!dir.exists(output_dir_plots)) dir.create(output_dir_plots, recursive = TRUE)

# Load the data
combined_methylation_data_compact <- readRDS(file.path(this.dir, "analysis_chromatin_state", "combined_methylation_data_compact_withpeaks.rds"))

#################################### CONSTANTS #####################################
start_script <- Sys.time()

# Define biosamples
biosamples <- c("A549", "GM12878", "HepG2", "K562")

#################################### PROCESS DATA ##################################
# Initialize an empty list to store data for each biosample
protein_counts_list <- list()

# Iterate through each biosample to parse protein hits and count occurrences
for (biosample in biosamples) {
  start_biosample <- Sys.time()
  
  protein_hits_col <- paste0("protein_hits_", biosample)
  
  # Process data for the current biosample
  biosample_data <- combined_methylation_data_compact %>%
    filter(!is.na(!!sym(protein_hits_col))) %>%
    mutate(Proteins = strsplit(!!sym(protein_hits_col), ";")) %>%
    unnest(Proteins) %>%
    mutate(Summarized_Proteins = sapply(Proteins, function(x) strsplit(x, "_")[[1]][1])) %>%  # Extract protein name
    group_by(Summarized_Proteins) %>%
    summarize(CpG_Count = n(), .groups = "drop") %>%
    mutate(Biosample = biosample)
  
  # Add to the list
  protein_counts_list[[biosample]] <- biosample_data
  
  end_biosample <- Sys.time()
  print(paste("Processed biosample:", biosample, "in", round(difftime(end_biosample, start_biosample, units = "secs"), 2), "seconds"))
}

# Combine all biosamples into a single data frame
protein_counts <- bind_rows(protein_counts_list)

#################################### PLOT DATA #####################################
# Plot the dodging bars for CpG counts by summarized proteins
cpg_plot <- ggplot(protein_counts, aes(x = reorder(Summarized_Proteins, -CpG_Count), y = CpG_Count, fill = Biosample)) +
  geom_bar(stat = "identity", position = "dodge") +  # Dodging bars for each biosample
  coord_flip() +
  scale_fill_manual(values = c("A549" = "#E76BF3", "GM12878" = "#A3A500", "HepG2" = "#F8766D", "K562" = "#00BFC4")) +
  labs(
    title = "Number of CpGs Covered by Each Summarized Protein Across Biosamples",
    x = "Proteins (Summarized)",
    y = "Number of CpGs",
    fill = "Biosample"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),
    axis.text.y = element_text(size = 7)
  )

# Save the plot
ggsave(
  filename = file.path(output_dir_plots, "CpG_coverage_by_summarized_protein_dodged_bars.svg"),
  plot = cpg_plot,
  width = 10,
  height = 80,
  limitsize = FALSE
)

#################################### END ###########################################
end_script <- Sys.time()
print(paste("Script completed in", round(difftime(end_script, start_script, units = "mins"), 2), "minutes"))

#################################### PIE PLOTS FOR CHROMATIN STATES #####################################
#################################### PIE PLOTS FOR CHROMATIN STATES #####################################

# Number of top proteins to highlight
top_n <- 8  # Adjust this variable to control how many top proteins are highlighted

# Initialize an empty list to store pie chart data for each chromatin state and biosample
chromatin_pie_list <- list()

# Define chromatin states
chromatin_states <- 1:18

# Iterate through each chromatin state and biosample
for (chromatin_state in chromatin_states) {
  for (biosample in biosamples) {
    protein_hits_col <- paste0("protein_hits_", biosample)
    chromatin_col <- paste0("Chromatin_State_", biosample)
    
    # Filter data for the current chromatin state and biosample
    state_data <- combined_methylation_data_compact %>%
      filter(!is.na(!!sym(protein_hits_col))) %>%
      filter(!!sym(chromatin_col) == chromatin_state) %>%
      mutate(Proteins = strsplit(!!sym(protein_hits_col), ";")) %>%
      unnest(Proteins) %>%
      mutate(Summarized_Proteins = sapply(Proteins, function(x) strsplit(x, "_")[[1]][1])) %>%
      group_by(Summarized_Proteins) %>%
      summarize(CpG_Count = n(), .groups = "drop") %>%
      arrange(desc(CpG_Count)) %>%
      mutate(
        Color_Group = ifelse(row_number() <= top_n, Summarized_Proteins, "Other")
      )  # Assign a "Color_Group" for top_n proteins and others
    
    # Store the pie chart data with additional metadata
    chromatin_pie_list[[paste0("State_", chromatin_state, "_", biosample)]] <- state_data %>%
      mutate(Chromatin_State = chromatin_state, Biosample = biosample)
  }
}

# Combine all pie chart data
chromatin_pie_data <- bind_rows(chromatin_pie_list)

# Define a custom color palette for proteins
unique_proteins <- unique(chromatin_pie_data$Color_Group)
top_colors <- scales::hue_pal()(top_n)  # Generate distinct colors for top_n
other_color <- "gray"
protein_colors <- setNames(c(top_colors, rep(other_color, length(unique_proteins) - top_n)), unique_proteins)

# Create pie charts as a grid
pie_plot <- ggplot(chromatin_pie_data, aes(x = "", y = CpG_Count, fill = Color_Group)) +
  geom_bar(width = 1, stat = "identity", position = "stack") +
  coord_polar(theta = "y") +
  facet_grid(rows = vars(Chromatin_State), cols = vars(Biosample), scales = "free") +
  scale_fill_manual(values = protein_colors) +
  labs(
    title = paste("Distribution of Protein Counts Across Chromatin States and Biosamples (Top", top_n, "Proteins)"),
    fill = "Proteins"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

# Save the pie chart grid plot
ggsave(
  filename = file.path(output_dir_plots, paste0("Chromatin_State_Protein_Pie_Plots_Top", top_n, "_Gray_Others.svg")),
  plot = pie_plot,
  width = 20,
  height = 30,
  limitsize = FALSE
)
