#################################################################
##  Motif Analysis classification
##  
##  input: fimo_single_experiments_on_known_motif_summary_table.txt
##  
##  output: Scatter plot of total ChIP Peaks vs motif matches     
##  v03 - 04.09.2024
##  Author: Daniel Batyrev 777634015
#################################################################

# Clear R working environment 
 rm(list=ls())

# Set the working directory to the directory where the script is located
this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this.dir)

# Load necessary libraries
library(ggplot2)
# Load the necessary library for pivoting data
library(tidyr)
library(dplyr)


# Load the data
df <- read.delim(file = file.path(this.dir, "fimo_single_experiments_on_known_motif_summary_table.txt"))

# Drop rows with NA values
df <- na.omit(df)

# Create a new column for the x-axis labels
df$label <- paste(df$protein, df$biosample, df$experiment, df$motif, sep = " ")

# Sort the data first by protein, then by biosample
df <- df %>%
  arrange(protein, biosample)

# Create the output directory if it doesn't exist
output_dir <- file.path(this.dir, "fimo_single_experiments_on_known_motif_plots")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}


# Create a scatter plot of total ChIP Peaks vs motif matches with x = y line
ggplot(df, aes(x = X..of.total.ChIP.Peaks, y = X..of.motif.matches, color = biosample)) +
  geom_point() +  # Scatter plot points colored by protein
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +  # x = y line
  theme_minimal() +             # Minimal theme
  labs(title = "Total ChIP Peaks vs Motif Matches",
       x = "Total ChIP Peaks",
       y = "Motif Matches") #+   # Axis labels
  #theme(plot.title = element_text(hjust = 0.5),  # Center the title
  #      legend.position = "none")  # Remove the legend


# Create a scatter plot of total ChIP Peaks vs motif matches with x = y line
ggplot(df, aes(x = X..of.motif.matches, y = df$X..of.motif.CpG.matches, color = biosample)) +
  geom_point() +  # Scatter plot points colored by protein
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +  # x = y line
  theme_minimal() +             # Minimal theme
  labs(title = "Total ChIP Peaks vs Motif Matches",
       x = "Total ChIP Peaks",
       y = "Motif Matches") #+   # Axis labels
#theme(plot.title = element_text(hjust = 0.5),  # Center the title
#      legend.position = "none")  # Remove the legend


# debug
# protein_name <- proteins[1]
# data <- df


#Function to generate plots per protein
generate_protein_plot <- function(protein_name, data) {
  # Filter data for the specific protein
  protein_data <- data %>% filter(protein == protein_name)

  # Create the plot
  p <- ggplot(protein_data, aes(x = reorder(label, interaction(biosample, experiment)), y = X..of.motif.matches)) +
    geom_bar(aes(fill = biosample), stat = "identity", alpha = 0.3) +  # Lighter fill for total motif matches
    geom_bar(aes(y = X..of.motif.CpG.matches, fill = biosample), stat = "identity", alpha = 1) +  # Darker fill for CpG matches
    theme_minimal() +  # Minimal theme
    labs(title = paste("Motif Matches for", protein_name),
         x = "Biosample, Experiment, and Motif",
         y = "Number of Matches") +  # Axis labels
    theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for readability
          plot.title = element_text(hjust = 0.5)) +  # Center the title
    scale_fill_manual(values = alpha(rainbow(length(unique(protein_data$biosample))), 0.8)) +  # Different colors for biosample
    coord_flip()  # Flip the coordinates to make it easier to read

  # Save the plot to a PNG file
  ggsave(filename = paste0(protein_name, "_motif_matches.png"), plot = p, path = output_dir, width = 10, height = 6)
}

# Loop through each protein and generate a plot
proteins <- unique(df$protein)
for (protein in proteins) {
  print(protein)
  generate_protein_plot(protein, df)
}



# # Function to generate bar plots per protein comparing total ChIP Peaks and motif matches
# generate_chip_motif_plot <- function(protein_name, data) {
#   # Filter data for the specific protein
#   protein_data <- data %>% filter(protein == protein_name)
# 
#   # Reshape the data for side-by-side bars (ChIP Peaks and Motif Matches)
#   protein_data_long <- protein_data %>%
#     tidyr::pivot_longer(cols = c(X..of.total.ChIP.Peaks, X..of.motif.matches),
#                  names_to = "Metric",
#                  values_to = "Value")
# 
#   # Create the plot
#   p <- ggplot(protein_data_long, aes(x = reorder(label, interaction(biosample, experiment)), y = Value, fill = biosample)) +
#     geom_bar(aes(alpha = Metric), position = "dodge", stat = "identity") +  # Dodge to make bars side by side
#     scale_alpha_manual(values = c(0.3, 1),
#                        labels = c("Total ChIP Peaks", "Motif Matches"),  # Add labels for light and dark bars
#                        name = "Metric") +  # Set legend title
#     theme_minimal() +  # Minimal theme
#     labs(title = paste("Total ChIP Peaks and Motif Matches for", protein_name),
#          x = "Biosample, Experiment, and Motif",
#          y = "Number of Peaks / Matches") +  # Axis labels
#     theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for readability
#           plot.title = element_text(hjust = 0.5)) +  # Center the title
#     coord_flip()  # Flip the coordinates to make it easier to read
# 
#   # Save the plot to a PNG file
#   ggsave(filename = paste0(protein_name, "_total_vs_motif_matches.png"), plot = p, path = output_dir, width = 10, height = 6)
# }
# 
# # Loop through each protein and generate a plot for total ChIP Peaks vs motif matches
# for (protein in proteins) {
#   print(paste("Generating total ChIP vs motif matches plot for:", protein))
#   generate_chip_motif_plot(protein, df)
# }



# Function to generate bar plots per protein comparing total ChIP Peaks, motif matches, and CpG motif matches
generate_chip_motif_plot <- function(protein_name, data) {
  # Filter data for the specific protein
  protein_data <- data %>% filter(protein == protein_name)
  
  # Reshape the data for side-by-side bars (ChIP Peaks and Motif Matches)
  protein_data_long <- protein_data %>%
    tidyr::pivot_longer(cols = c(X..of.total.ChIP.Peaks, X..of.motif.matches),
                        names_to = "Metric",
                        values_to = "Value")
  
  # Create the overlay data for CpG Matches
  # Set 0 for ChIP Peaks, and set the actual CpG matches for Motif Matches
  protein_data_long$CpG_Overlay <- ifelse(protein_data_long$Metric == "X..of.motif.matches", 
                                          protein_data$X..of.motif.CpG.matches, 0)
  
  # # Create the plot
  # p <- ggplot(protein_data_long, aes(x = reorder(label, interaction(biosample, experiment)),  fill = biosample)) +
  #   geom_bar(aes(alpha = Metric,y = CpG_Overlay), position = "dodge", stat = "identity") +  # Dodge to make bars side by side
  #   geom_bar(aes(alpha = Metric,y = Value), position = "dodge", stat = "identity") + 
  #   #geom_bar(aes(y = Value,aes(alpha = Metric)), position = "dodge", stat = "identity",fill = "black") +  # Dodge to make bars side by side
  #   #Overlay CpG matches on top of both bars, but only make the motif bar visible
  #   # geom_bar(aes(y = CpG_Overlay), position = position_dodge(width = 0.9), stat = "identity",
  #   #          fill = "black", alpha = 0.6) +  # CpG overlay with transparency
  #   scale_alpha_manual(values = c(0.3, 1),
  #                      labels = c("Total ChIP Peaks", "Motif Matches"),  # Add labels for light and dark bars
  #                      name = "Metric") +  # Set legend title
  #   theme_minimal() +  # Minimal theme
  #   labs(title = paste("Total ChIP Peaks and Motif Matches for", protein_name),
  #        x = "Biosample, Experiment, and Motif",
  #        y = "Number of Peaks / Matches") +  # Axis labels
  #   theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for readability
  #         plot.title = element_text(hjust = 0.5)) +  # Center the title
  #   coord_flip()  # Flip the coordinates to make it easier to read
  
  # Create the plot
  p <- ggplot(protein_data_long, aes(x = reorder(label, interaction(biosample, experiment)), fill = biosample)) +
    
    # Total motif matches and ChIP peaks
    geom_bar(aes(alpha = Metric, y = Value), position = "dodge", stat = "identity") +
    
    # CpG overlay bar
    geom_bar(aes(alpha = Metric, y = CpG_Overlay), position = "dodge", stat = "identity") +  # Dodge to make bars side by side

    # Custom alpha scale for shading between Total ChIP Peaks and Motif Matches
    scale_alpha_manual(values = c(0.3, 1),
                       labels = c( "# Motif Matches","# ChIP Peaks"),  # Add labels for light and dark bars
                       name = "Metric") +  # Set legend title
    
    # Add a manual legend entry for CpG matches (overlay inside motif matches)
    guides(fill = guide_legend(title = "Sample"),
           alpha = guide_legend(override.aes = list(fill = "gray"))) +  # Custom legend for CpG matches
    
    theme_minimal() +  # Minimal theme
    
    # Custom labels for the plot
    labs(title = paste("Total ChIP Peaks and Motif Matches for", protein_name),
         subtitle = "Motif Matches contain two shades: CpG matches (dark) and non-CpG matches (light)",  # Add a descriptive subtitle
         x = "Biosample, Experiment, and Motif",
         y = "Number of Peaks / Matches") +  # Axis labels
    
    # Rotate x-axis labels for readability
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5)) +  # Center the title and subtitle
    
    # Flip the coordinates for better readability
    coord_flip()
  
  # Save the plot to a PNG file
  ggsave(filename = paste0(protein_name, "_total_vs_motif_matches_with_CpG.png"), plot = p, path = output_dir, width = 10, height = 6)
}

# Loop through each protein and generate a plot for total ChIP Peaks vs motif matches and CpG matches
for (protein in proteins) {
  print(paste("Generating total ChIP vs motif matches plot for:", protein))
  generate_chip_motif_plot(protein, df)
}

