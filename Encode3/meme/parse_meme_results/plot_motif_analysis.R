# Documentation -----------------------------------------------------------
##
##
##  input: 
##
##
##  output: 
##
##
##  v_01 21.07.2024
##  Author: Daniel Batyrev (HUJI 777634015)
##
# Set up Work Environment --------------------------------------------------
# Clear R working environment
rm(list = ls())
cluster <- FALSE
if (cluster) {
  this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/HumanEvo/HumanEvo/"
  picuture_file_extension <- "pdf"
} else{
  this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  picuture_file_extension <- "png"
}
setwd(this.dir)
detachAllPackages <- function() {
  basic.packages <-
    c(
      "package:stats",
      "package:graphics",
      "package:grDevices",
      "package:utils",
      "package:datasets",
      "package:methods",
      "package:base"
    )
  package.list <-
    search()[ifelse(unlist(gregexpr("package:", search())) == 1, TRUE, FALSE)]
  package.list <- setdiff(package.list, basic.packages)
  if (length(package.list) > 0)
    for (package in package.list)
      detach(package, character.only = TRUE)
}
detachAllPackages()
#################################################################################

# Install necessary packages if not already installed
# install.packages(c("dplyr", "ggplot2", "reshape2", "dendextend", "ggdendro"))

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(reshape2)
# library(dendextend)
# library(ggdendro)

# Set the input folder and file name
input_folder <- "C:\\Users\\Daniel Batyrev\\Documents\\GitHub\\methylation_vs_chromatin_vs_ChIP\\Encode3\\meme\\meme-chip_results"
file_name <- "all_motifs_groups_all_with_CG.csv"
file_path <- file.path(input_folder, file_name)

# Load the data
data <- read.csv(file_path)

# Plot and save the histogram of Distribution
p <- ggplot(data = data, mapping = aes(x = Distribution)) +
  geom_histogram(binwidth = 0.001) +
  xlim(0, 0.1) +
  labs(title = "Histogram of Distribution",
       x = "Distribution of concentration values",
       y = "Frequency") +
  theme_minimal()
ggsave(filename = file.path(this.dir, paste0("Histogram_of_Distribution.", picuture_file_extension)), plot = p)

# Filter the data
data_concentrated <- data[data$Distribution > 0.05,]

# Replace NA values with zero and select the highest concentration value for each protein
highest_concentration_data <- data %>%
  mutate(Distribution = ifelse(is.na(Distribution), 0, Distribution)) %>%
  group_by(Protein) %>%
  slice_max(order_by = Distribution, n = 1, with_ties = FALSE)

# Plot and save the histogram of highest concentration distribution by protein
p2 <- ggplot(data = highest_concentration_data, mapping = aes(x = Distribution)) +
  geom_histogram(binwidth = 0.01) +
  labs(title = "Histogram of Highest Concentration Distribution by Protein",
       x = "Distribution",
       y = "Frequency") +
  theme_minimal()
ggsave(filename = file.path(this.dir, paste0("Histogram_of_Highest_Concentration_Distribution_by_Protein.", picuture_file_extension)), plot = p2)

# Categorize the data
highest_concentration_data <- highest_concentration_data %>%
  mutate(Category = ifelse(Distribution > 0.05, "Above 0.05", "0.05 or Below"))

# Create a summary for the pie chart
pie_data <- highest_concentration_data %>%
  group_by(Category) %>%
  summarise(Count = n())

# Plot and save the pie chart with absolute numbers
p3 <- ggplot(data = pie_data, mapping = aes(x = "", y = Count, fill = Category)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  labs(title = "Proteins with Highest Concentration Values",
       x = "",
       y = "") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5))
ggsave(filename = file.path(this.dir, paste0("Proteins_with_Highest_Concentration_Values.", picuture_file_extension)), plot = p3)


# Ensure meme.chip.E.value is present and log-transform it
highest_concentration_data <- highest_concentration_data %>%
  filter(!is.na(meme.chip.E.value)) %>%
  mutate(log_E_value = log10(meme.chip.E.value))

# Plot and save the histogram of log-transformed meme.chip.E.value with x-axis limits set from -50 to 0
p4 <- ggplot(data = highest_concentration_data, mapping = aes(x = log_E_value)) +
  geom_histogram(binwidth = 1) +
  xlim(-50, 0) +
  labs(title = "Histogram of Log-transformed E-values",
       x = "Log10(E-value)",
       y = "Frequency") +
  geom_vline(xintercept = log(0.05/333), color = "red")+
  theme_minimal()

ggsave(filename = file.path(this.dir, paste0("Histogram_of_Log_transformed_E_values.", picuture_file_extension)), plot = p4)


# Define the threshold
threshold <- 0.05 / 333

# Categorize the proteins based on the threshold
highest_concentration_data <- highest_concentration_data %>%
  mutate(Category = ifelse(meme.chip.E.value < threshold, "Below Threshold", "Above Threshold"))


# Create a summary for the pie chart
pie_data <- highest_concentration_data %>%
  group_by(Category) %>%
  summarise(Count = n())

# Plot and save the pie chart with absolute numbers
p5 <- ggplot(data = pie_data, mapping = aes(x = "", y = Count, fill = Category)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  labs(title = "Proteins with E-values",
       x = "",
       y = "") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5))
ggsave(filename = file.path(this.dir, paste0("Proteins_with_E_values_Below_and_Above_Threshold.", picuture_file_extension)), plot = p5)



# Replace NA values with zero and select the highest concentration value for each protein
highest_concentration_data <- data %>%
  mutate(Distribution = ifelse(is.na(Distribution), 0, Distribution)) %>%
  group_by(Protein) %>%
  slice_max(order_by = Distribution, n = 1, with_ties = FALSE) %>%
  filter(!is.na(meme.chip.E.value))

# Categorize the proteins based on the threshold and concentration value
highest_concentration_data <- highest_concentration_data %>%
  mutate(Category = case_when(
    meme.chip.E.value < threshold & Distribution > 0.05 ~ "E-value < Threshold & Concentration > 0.05",
    meme.chip.E.value < threshold & Distribution <= 0.05 ~ "E-value < Threshold & Concentration <= 0.05",
    meme.chip.E.value >= threshold & Distribution > 0.05 ~ "E-value >= Threshold & Concentration > 0.05",
    TRUE ~ "E-value >= Threshold & Concentration <= 0.05"
  ))

# Create a summary for the pie chart
pie_data <- highest_concentration_data %>%
  group_by(Category) %>%
  summarise(Count = n())

# Plot and save the pie chart with absolute numbers
p6 <- ggplot(data = pie_data, mapping = aes(x = "", y = Count, fill = Category)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  labs(title = "Proteins by E-value and Concentration",
       x = "",
       y = "") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank()) +
  geom_text(aes(label = Count), position = position_stack(vjust = 0.5))

ggsave(filename = file.path(this.dir, paste0("Proteins_by_E_value_and_Concentration.", picuture_file_extension)), plot = p6)

library(tidyr)

##############################################
# Define the threshold
threshold <- 0.05 / 333

# Define the list of target motifs
target_motifs <- c("CTCF_HUMAN.H11MO.0.A", "MYCN_HUMAN.H11MO.0.A", "FOSL1_HUMAN.H11MO.0.A","ELK4_HUMAN.H11MO.0.A","CREB1_HUMAN.H11MO.0.A")

# Filter significant data
significant_data <- data[data$Distribution > 0.05 & data$meme.chip.E.value < threshold,]

# Split the `Known.or.Similar.Motifs` into individual motifs and create a long format dataframe
motifs_long <- significant_data %>%
  select(Protein, Group.ID, Known.or.Similar.Motifs) %>%
  mutate(Known.or.Similar.Motifs = strsplit(as.character(Known.or.Similar.Motifs), ";")) %>%
  unnest(Known.or.Similar.Motifs)

# Initialize a color palette and a counter for assigning colors
color_palette <- rainbow(length(target_motifs))
motif_colors <- setNames(color_palette, target_motifs)

# Initialize an empty column for motif relationship
motifs_long$is_related <- "Not Related"

# Iterate through the target motifs to mark related groups and assign colors
for (motif in target_motifs) {
  related_groups <- motifs_long %>%
    filter(Known.or.Similar.Motifs == motif) %>%
    pull(Group.ID)
  motifs_long <- motifs_long %>%
    mutate(is_related = ifelse(Group.ID %in% related_groups, motif, is_related))
}

# Count unique occurrences of each motif in different proteins
motif_counts <- motifs_long %>%
  group_by(Known.or.Similar.Motifs, is_related) %>%
  summarise(Protein_Count = n_distinct(Protein)) %>%
  arrange(desc(Protein_Count))

# Filter the top 50 motifs
top_50_motifs <- motif_counts %>%
  top_n(50, Protein_Count) %>%
  arrange(desc(Protein_Count))

# Define a color mapping function
get_color <- function(is_related) {
  if (is_related %in% target_motifs) {
    return(motif_colors[is_related])
  } else {
    return("gray")
  }
}

top_50_motifs$color <- sapply(top_50_motifs$is_related, get_color)

# Plot the top 50 motifs as a bar plot with different colors for related motifs
p <- ggplot(top_50_motifs[1:50,], aes(x = reorder(Known.or.Similar.Motifs, -Protein_Count), y = Protein_Count, fill = is_related)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(motif_colors, "Not Related" = "gray")) +
  coord_flip() +
  labs(title = "Top 50 Motifs by Protein Count",
       x = "Motif",
       y = "Protein Count",
       fill = "Motif Relationship") +
  theme_minimal()

# Print the plot
print(p)

# Save the plot
ggsave(filename = file.path(this.dir, paste0("Top_50_Motifs_by_Protein_Count.", picuture_file_extension)),
       plot = p,width = 12,height = 8)