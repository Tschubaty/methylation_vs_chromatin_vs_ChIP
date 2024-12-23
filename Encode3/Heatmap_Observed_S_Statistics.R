#################################################################
##  Script: Heatmap_Observed_S_Statistics
##  Description: Generate a heatmap to visualize different p-values 
##               (`p_value_S_bigger`, `p_value_S_smaller`, 
##                `p_value_S_bigger_same_state`, `p_value_S_smaller_same_state`) 
##               as tiles along the x-axis. Rows are distinct by 
##               `experiment_ids` or `motif`.
##
##  Input: RDS file containing `summary_test` data frame.
##         File path: 
##         C:/Users/Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/simulation/summary_test.rds
##
##  Output: Heatmap saved in a directory named after the script.
##
##  Version: 22.12.2024
##  Author: Daniel Batyrev (777634015)
#################################################################

# Clear R working environment
rm(list = ls())
cluster <- FALSE

# Set script name and working directory
script_name <- "Heatmap_Observed_S_Statistics"
this.dir <- if (cluster) "/ems/elsc-labs/meshorer-e/daniel.batyrev/Encode3/" else dirname(rstudioapi::getSourceEditorContext()$path)

setwd(this.dir)

#################################### LIBRARIES ########################################
library(ggplot2)
library(tidyr)
library(dplyr)

##################################### INPUT ########################################
# Define input and output paths
input_file <- "C:/Users/Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/simulation/summary_test.rds"
output_dir <- file.path(this.dir, script_name)
output_plots_dir <- file.path(output_dir, "plots")

# Create directories if they don't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(output_plots_dir)) dir.create(output_plots_dir, recursive = TRUE)

# Load the data
summary_test <- readRDS(input_file)

#################################### PROCESS DATA #####################################
# Prepare the data for the heatmap
heatmap_data <- summary_test %>%
  select(protein, biosample1, biosample2, experiment_ids, motif, 
         p_value_S_bigger, p_value_S_smaller, 
         p_value_S_bigger_same_state, p_value_S_smaller_same_state) %>%
  # Pivot data longer to create tiles for different p-value types
  pivot_longer(
    cols = starts_with("p_value_"),
    names_to = "P_Value_Type",
    values_to = "P_Value"
  ) %>%
  # Combine columns for y-axis labels
  mutate(
    Row_Label = paste(protein, biosample1, biosample2, experiment_ids, motif, sep = " | ")
  )

#################################### VISUALIZE DATA ###################################
# Create the heatmap
heatmap_plot <- ggplot(heatmap_data, aes(x = P_Value_Type, y = Row_Label, fill = P_Value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",   # Low p-values (significant)
    mid = "white", 
    high = "red",   # High p-values (non-significant)
    midpoint = 0.5,
    name = "p-value"
  ) +
  labs(
    title = "Heatmap of Observed S Statistics Across Tests",
    x = "P-Value Type",
    y = "Protein | Biosample1 | Biosample2 | Experiment IDs | Motif"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

#################################### SAVE PLOT ########################################
# Save the heatmap as an SVG
output_file <- file.path(output_plots_dir, paste0(script_name, "_heatmap.svg"))
ggsave(
  filename = output_file,
  plot = heatmap_plot,
  width = 14,
  height = 10 + nrow(summary_test) * 0.1,
  limitsize = FALSE # Adjust height based on the number of rows
)

print(paste("Heatmap saved to:", output_file))

#################################### END ##############################################
print("Script execution completed successfully.")
