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
#input_file <- "C:/Users/Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/simulation/summary_test.rds"
input_file <-  "D:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/permutation_test/stratified_test_100k.rds"
output_dir <- file.path(this.dir, script_name)
output_plots_dir <- file.path(output_dir, "plots")

# Create directories if they don't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(output_plots_dir)) dir.create(output_plots_dir, recursive = TRUE)

# Load the data
summary_test <- readRDS(input_file)

#stratified_test_100k.rds

#################################### PROCESS DATA #####################################
# # Prepare the data for the heatmap
# heatmap_data <- summary_test %>%
#   select(protein, biosample1, biosample2, experiment_ids, motif, 
#          p_value_S_bigger, p_value_S_smaller, 
#          p_value_S_bigger_same_state, p_value_S_smaller_same_state) %>%
#   # Pivot data longer to create tiles for different p-value types
#   pivot_longer(
#     cols = starts_with("p_value_"),
#     names_to = "P_Value_Type",
#     values_to = "P_Value"
#   ) %>%
#   # Combine columns for y-axis labels
#   mutate(
#     Row_Label = paste(protein, biosample1, biosample2, experiment_ids, motif, sep = " | ")
#   )

# Prepare the data for the heatmap
heatmap_data <- summary_test %>%
  select(protein, biosample1, biosample2, motif, 
         p_value_S_bigger, p_value_S_smaller) %>%
  # Pivot data longer to create tiles for different p-value types
  pivot_longer(
    cols = starts_with("p_value_"),
    names_to = "P_Value_Type",
    values_to = "P_Value"
  ) %>%
  # Combine columns for y-axis labels
  mutate(
    Row_Label = paste(protein, biosample1, biosample2, motif, sep = " | ")
  )

library(dplyr)

# Combine p_value_S_bigger and p_value_S_smaller into two-sided p-values and remember the direction
heatmap_data <- heatmap_data %>%
  group_by(protein, biosample1, biosample2, motif, Row_Label) %>%
  summarise(
    P_Value_Two_Sided = min(1, 2 * min(P_Value[P_Value_Type == "p_value_S_bigger"], 
                                       P_Value[P_Value_Type == "p_value_S_smaller"])),
    Direction = ifelse(
      P_Value[P_Value_Type == "p_value_S_bigger"] < P_Value[P_Value_Type == "p_value_S_smaller"],
      "S_bigger", "S_smaller"
    )
  ) %>%
  ungroup()

# Apply FDR correction to the two-sided p-values
heatmap_data <- heatmap_data %>%
  mutate(FDR_Adjusted_P_Value = p.adjust(P_Value_Two_Sided, method = "BH"))

# Display the results
print(heatmap_data)


#################################### VISUALIZE DATA ###################################
# # Create the heatmap
# heatmap_plot <- ggplot(heatmap_data, aes(x = P_Value_Type, y = Row_Label, fill = P_Value)) +
#   geom_tile(color = "white") +
#   scale_fill_gradient2(
#     low = "blue",   # Low p-values (significant)
#     mid = "white", 
#     high = "red",   # High p-values (non-significant)
#     midpoint = 0.5,
#     name = "p-value"
#   ) +
#   labs(
#     title = "Heatmap of Observed S Statistics Across Tests",
#     x = "P-Value Type",
#     y = "Protein | Biosample1 | Biosample2 | Experiment IDs | Motif"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     axis.text.y = element_text(size = 7),
#     legend.title = element_text(size = 10),
#     legend.text = element_text(size = 8)
#   )
#################################### TRANSFORM AND VISUALIZE 1D HEATMAP ###################################
# Add transformed values based on BH-corrected p-values
heatmap_data <- heatmap_data %>%
  mutate(
    Transformed_P_Value = ifelse(Direction == "S_bigger", 1 - FDR_Adjusted_P_Value, -1 + FDR_Adjusted_P_Value)
  )

# Create the 1D heatmap
heatmap_1d_plot <- ggplot(heatmap_data, aes(x = 1, y = Row_Label, fill = Transformed_P_Value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",   # Significant negative deviations
    mid = "white",  # Non-significant values close to zero
    high = "red",   # Significant positive deviations
    midpoint = 0,   # Zero as the midpoint
    name = "Transformed BH p-value"
  ) +
  labs(
    title = "1D Heatmap of Transformed BH-Corrected S Statistics",
    x = NULL,
    y = "Protein | Biosample1 | Biosample2 | Experiment IDs | Motif"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove X-axis text
    axis.ticks.x = element_blank(), # Remove X-axis ticks
    axis.text.y = element_text(size = 7),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

# Display the heatmap
print(heatmap_1d_plot)


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
