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
# rm(list=ls())

# Set the working directory to the directory where the script is located
this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this.dir)

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

# Load necessary libraries
library(ggplot2)
# Load the necessary library for pivoting data
library(tidyr)

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



# Function to generate bar plots per protein comparing total ChIP Peaks and motif matches
generate_chip_motif_plot <- function(protein_name, data) {
  # Filter data for the specific protein
  protein_data <- data %>% filter(protein == protein_name)

  # Reshape the data for side-by-side bars (ChIP Peaks and Motif Matches)
  protein_data_long <- protein_data %>%
    pivot_longer(cols = c(X..of.total.ChIP.Peaks, X..of.motif.matches),
                 names_to = "Metric",
                 values_to = "Value")

  # Create the plot
  p <- ggplot(protein_data_long, aes(x = reorder(label, interaction(biosample, experiment)), y = Value, fill = biosample)) +
    geom_bar(aes(alpha = Metric), position = "dodge", stat = "identity") +  # Dodge to make bars side by side
    scale_alpha_manual(values = c(0.3, 1),
                       labels = c("Total ChIP Peaks", "Motif Matches"),  # Add labels for light and dark bars
                       name = "Metric") +  # Set legend title
    theme_minimal() +  # Minimal theme
    labs(title = paste("Total ChIP Peaks and Motif Matches for", protein_name),
         x = "Biosample, Experiment, and Motif",
         y = "Number of Peaks / Matches") +  # Axis labels
    theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for readability
          plot.title = element_text(hjust = 0.5)) +  # Center the title
    coord_flip()  # Flip the coordinates to make it easier to read

  # Save the plot to a PNG file
  ggsave(filename = paste0(protein_name, "_total_vs_motif_matches.png"), plot = p, path = output_dir, width = 10, height = 6)
}

# Loop through each protein and generate a plot for total ChIP Peaks vs motif matches
for (protein in proteins) {
  print(paste("Generating total ChIP vs motif matches plot for:", protein))
  generate_chip_motif_plot(protein, df)
}


# Function to generate bar plots per protein comparing total ChIP Peaks, motif matches, and CpG motif matches
generate_chip_motif_plot <- function(protein_name, data) {
  # Filter data for the specific protein
  protein_data <- data %>% filter(protein == protein_name)
  
  # Reshape the data for side-by-side bars (ChIP Peaks and Motif Matches)
  protein_data_long <- protein_data %>%
    pivot_longer(cols = c(X..of.total.ChIP.Peaks, X..of.motif.matches), 
                 names_to = "Metric", 
                 values_to = "Value")
  
  # Create stacked bar data for motif matches and CpG matches
  protein_data_long_stacked <- protein_data %>%
    mutate(Non_CpG_Matches = X..of.motif.matches - X..of.motif.CpG.matches) %>%  # Calculate non-CpG matches
    pivot_longer(cols = c(X..of.motif.CpG.matches, Non_CpG_Matches), 
                 names_to = "CpG_Status", 
                 values_to = "CpG_Matches_Value")
  
  # Create the plot
  p <- ggplot(protein_data_long, aes(x = reorder(label, interaction(biosample, experiment)), y = Value, fill = biosample)) +
    geom_bar(aes(alpha = Metric), position = "dodge", stat = "identity") +  # Side-by-side bars for total ChIP Peaks and total Motif Matches
    geom_bar(data = protein_data_long_stacked, 
             aes(x = reorder(label, interaction(biosample, experiment)), y = CpG_Matches_Value, fill = biosample, alpha = CpG_Status), 
             stat = "identity", position = "stack") +  # Stacked bar for CpG and Non-CpG motif matches
    scale_alpha_manual(values = c(0.2, 0,0,0),  # Lighter for non-CpG Matches, darker for CpG matches
                       breaks = c("Non_CpG_Matches", "X..of.motif.CpG.matches"),  # Match CpG and non-CpG matches
                       labels = c("Non-CpG Matches", "CpG Matches"),  # Label the stacked sections
                       name = "Match Type") +  # Set legend title
    theme_minimal() +  # Minimal theme
    labs(title = paste("Total ChIP Peaks and Motif Matches for", protein_name), 
         x = "Biosample, Experiment, and Motif", 
         y = "Number of Peaks / Matches") +  # Axis labels
    theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels for readability
          plot.title = element_text(hjust = 0.5)) +  # Center the title
    coord_flip()  # Flip the coordinates to make it easier to read
  
  # Save the plot to a PNG file
  ggsave(filename = paste0(protein_name, "_total_vs_motif_matches_with_CpG.png"), plot = p, path = output_dir, width = 10, height = 6)
}

# Loop through each protein and generate a plot for total ChIP Peaks vs motif matches and CpG matches
for (protein in proteins) {
  print(paste("Generating total ChIP vs motif matches plot for:", protein))
  generate_chip_motif_plot(protein, df)
}
