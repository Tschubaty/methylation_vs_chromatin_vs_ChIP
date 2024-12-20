#################################################################
##  beta
##
##  input: everything in HG38
##
##          Encode3\
##
##
##  output:   Encode3\
##
##  v_1 29.09.2024
##  Author: Daniel Batyrev 777634015
#################################################################
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
#################################### Libs ########################################
library(foreach)
library(doParallel)
library(readr)
library(dplyr)
library(rlang)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(tidyr)
library(readr)
##################################### INPUT ########################################
# Define the input directory for chromosome files
input_WGBS_dir <- file.path(this.dir,"WGBS/byChr/concat_methylation")


input_ChIP_dir <- file.path(
  this.dir,
  "meme",
  "fimo_single_experiments_on_known_motif_and_peaks_add_methylation_v3"
)

# Set up parallel processing
n_cores <- detectCores() - 1
registerDoParallel(cores = n_cores)


# Define the directory where you want to save the plots
output_dir_plots <- file.path(this.dir,"plots")

# Ensure the output directory exists
if (!dir.exists(output_dir_plots)) {
  dir.create(output_dir_plots, recursive = TRUE)
}

# Define the output directory name
output_analysis_dir <- file.path(this.dir, "analysis_chromatin_state")

# Create the directory if it doesn't already exist
if (!dir.exists(output_analysis_dir)) {
  dir.create(output_analysis_dir, recursive = TRUE)
}

# Save the directory location in a variable
out_dir <- output_analysis_dir

################################## constants #####################################
start_script <- Sys.time()
# Suggested colors for the biosamples
group.colors <- c(
  HepG2    = "#F8766D",
  # A warm reddish color
  K562     = "#00BFC4",
  # A cool cyan color
  GM12878  = "#A3A500",
  # A yellow-green color, good contrast with the others
  A549     = "#E76BF3"   # A vibrant purple color
)


# set up parallel processing
n_cores <- detectCores() - 1

CHR_NAMES <- paste0("chr", c(1:22))

# Define chromatin state names
chromatin_state_names <- c(
  "1_TssA",       # Active TSS
  "2_TssFlnk",    # Flanking Active TSS
  "3_TssFlnkU",   # Flanking TSS Upstream
  "4_TssFlnkD",   # Flanking TSS Downstream
  "5_Tx",         # Strong Transcription
  "6_TxWk",       # Weak Transcription
  "7_EnhG1",      # Genic Enhancer 1
  "8_EnhG2",      # Genic Enhancer 2
  "9_EnhA1",      # Active Enhancer 1
  "10_EnhA2",     # Active Enhancer 2
  "11_EnhWk",     # Weak Enhancer
  "12_ZNF_Rpts",  # ZNF Genes & Repeats
  "13_Het",       # Heterochromatin
  "14_TssBiv",    # Bivalent TSS
  "15_EnhBiv",    # Bivalent Enhancer
  "16_ReprPC",    # Repressed Polycomb
  "17_ReprPCWk",  # Weak Repressed Polycomb
  "18_Quies"      # Quiescent/Low Activity
)

# Define corresponding hex color codes for each chromatin state
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
  "#000000"   # 18_Quies     -> Black
)

# Define chromatin state colors
chromatin_state_colors_short <- c(
  "1" = "#FF0000",  "2" = "#FF4500",  "3" = "#FF9900",  "4" = "#FFCC00",
  "5" = "#00CC00",  "6" = "#006400",  "7" = "#FFD700",  "8" = "#FFD700",
  "9" = "#FFFF00",  "10" = "#FFDD00", "11" = "#FFEA73", "12" = "#9370DB",
  "13" = "#C0C0C0", "14" = "#FF4500", "15" = "#FFDD00", "16" = "#808080",
  "17" = "#A9A9A9", "18" = "#000000"
)


# Create a named vector to map chromatin states to their colors
chromatin_state_map <- setNames(chromatin_state_colors, chromatin_state_names)

# View the mapping (optional)
chromatin_state_map

# Function to extract and clean column names from the first row of a text file
extract_colnames <- function(file_path) {
  # Read the first row
  first_row <- readLines(file_path, n = 1)
  
  # Remove the leading "#" and replace double tabs with single tab
  cleaned_row <- gsub("#", "", first_row)          # Remove leading #
  
  # Split the cleaned string into individual column names
  colnames <- unlist(strsplit(cleaned_row, "\t"))
  
  return(colnames)
}

# Function to read .bed files with dynamically extracted column names and correct column types
read_bed_file <- function(file_path) {
  # Extract the correct column names from the first row
  colnames <- extract_colnames(file_path)
  
  # Read the BED file, skipping the first line (header row)
  df <- readr::read_delim(
    file_path,
    delim = "\t",
    skip = 1,
    col_names = colnames,
    show_col_types = FALSE
  )
  
  return(df)
}


#######################################################################################

# v01 unknown only one categroy 
# 
# summarize_chromatin_states <- function(file_path, biosample) {
#   df_ChIP <- read_bed_file(file_path)
#   
#   # Extract the column for the biosample's chromatin state
#   chromatin_col <- paste0("Chromatin_State_", biosample)
#   
#   # Replace missing or "NA" states with "Unknown"
#   df_ChIP[[chromatin_col]] <- ifelse(
#     is.na(df_ChIP[[chromatin_col]]) | df_ChIP[[chromatin_col]] == ".",
#     "Unknown",
#     df_ChIP[[chromatin_col]]
#   )
#   
#   # Summarize the distribution of chromatin states
#   summary <- df_ChIP %>%
#     mutate(!!sym(chromatin_col) := ifelse(
#       is.na(!!sym(chromatin_col)) | !!sym(chromatin_col) == ".",
#       "Unknown",  # Assign a label for missing or NA states
#       !!sym(chromatin_col)
#     )) %>%
#     group_by(!!sym(chromatin_col)) %>%
#     summarize(Count = n(), .groups = "drop") %>%
#     arrange(desc(Count))
#   
#   # Ensure all chromatin states are present, including "Unknown"
#   chromatin_states <- tibble(
#     Chromatin_State = c(names(chromatin_state_colors_short), "Unknown")
#   )
#   
#   summary <- chromatin_states %>%
#     left_join(summary %>% rename(Chromatin_State = !!sym(chromatin_col)), 
#               by = "Chromatin_State") %>%
#     mutate(Count = replace_na(Count, 0)) %>%
#     arrange(match(Chromatin_State, c(names(chromatin_state_colors_short), "Unknown")))
#   
#   # Add percentages column to the summary
#   summary <- summary %>%
#     mutate(
#       Percentage = (Count / sum(Count)) * 100
#     )
#   
#   
#   return(summary)
# }
# 
# 
# # Initialize a list to store summaries for all experiments
# all_experiment_summaries <- list()
# 
# # Iterate through all proteins and their respective experiments
# proteins <- list.files(path = input_ChIP_dir)
# for (protein in proteins) {
#   print(protein)
#   # skil log folder 
#   if (protein == "log") next
#   # Get all files for the protein
#   protein_dir <- file.path(input_ChIP_dir, protein)
#   experiment_files <- list.files(path = protein_dir, full.names = TRUE)
#   
#   for (file_path in experiment_files) {
#     # Extract experiment details from file name
#     file_name <- basename(file_path)
#     experiment_id <- strsplit(file_name, "_")[[1]][3]
#     biosample <- strsplit(file_name, "_")[[1]][1]
#     
#     # Summarize chromatin states for this experiment
#     summary <- summarize_chromatin_states(file_path, biosample)
#     
#     # Add experiment details
#     summary <- summary %>%
#       mutate(Protein = protein, ExperimentID = experiment_id, Biosample = biosample)
#     
#     # Append to the list
#     all_experiment_summaries[[paste0(protein, "_", experiment_id)]] <- summary
#   }
# }
# 
# # Combine all summaries into a single data frame
# final_summary <- bind_rows(all_experiment_summaries)
# 
# # Ensure Chromatin_State is ordered
# final_summary <- final_summary %>%
#   mutate(
#     Chromatin_State = factor(
#       Chromatin_State,
#       levels = c(names(chromatin_state_colors_short), "Unknown")
#     )
#   )
# 
# # Save the summary to a CSV file
# output_csv <- file.path(out_dir, "chromatin_state_summary.csv")
# write_csv(final_summary, output_csv)

################################################################################
# details unknown categtory 
summarize_chromatin_states <- function(file_path, biosample) {
  df_ChIP <- read_bed_file(file_path)
  
  # Extract the column for the biosample's chromatin state
  chromatin_col <- paste0("Chromatin_State_", biosample)
  
  # Categorize rows based on start_motif, start_cg, and chromatin state
  df_ChIP <- df_ChIP %>%
    mutate(
      !!sym(chromatin_col) := case_when(
        is.na(!!sym(chromatin_col)) | !!sym(chromatin_col) == "." & start_motif == -1 ~ "No motif in peak",
        is.na(!!sym(chromatin_col)) | !!sym(chromatin_col) == "." & start_motif != -1 & start_cg == 0 ~ "No CG in motif",
        is.na(!!sym(chromatin_col)) | !!sym(chromatin_col) == "." & start_motif != -1 & start_cg != 0 ~ "No State assignment for CG",
        TRUE ~ !!sym(chromatin_col)  # Keep original value otherwise
      )
    )
  
  # Summarize the distribution of chromatin states
  summary <- df_ChIP %>%
    group_by(!!sym(chromatin_col)) %>%
    summarize(Count = n(), .groups = "drop") %>%
    arrange(desc(Count))
  
  # Ensure all chromatin states are present
  chromatin_states <- tibble(
    Chromatin_State = c(names(chromatin_state_colors_short), 
                        "No motif in peak", "No CG in motif", "No State assignment for CG")
  )
  
  summary <- chromatin_states %>%
    left_join(summary %>% rename(Chromatin_State = !!sym(chromatin_col)), 
              by = "Chromatin_State") %>%
    mutate(Count = replace_na(Count, 0)) %>%
    arrange(match(Chromatin_State, c(names(chromatin_state_colors_short), 
                                     "No motif in peak", "No CG in motif", "No State assignment for CG")))
  
  # Add percentages column to the summary
  summary <- summary %>%
    mutate(
      Percentage = (Count / sum(Count)) * 100
    )
  
  return(summary)
}

# Initialize a list to store summaries for all experiments
all_experiment_summaries <- list()

# Iterate through all proteins and their respective experiments
proteins <- list.files(path = input_ChIP_dir)
for (protein in proteins) {
  print(protein)
  # Skip log folder
  if (protein == "log") next
  # Get all files for the protein
  protein_dir <- file.path(input_ChIP_dir, protein)
  experiment_files <- list.files(path = protein_dir, full.names = TRUE)
  
  for (file_path in experiment_files) {
    # Extract experiment details from file name
    file_name <- basename(file_path)
    experiment_id <- strsplit(file_name, "_")[[1]][3]
    biosample <- strsplit(file_name, "_")[[1]][1]
    
    # Summarize chromatin states for this experiment
    summary <- summarize_chromatin_states(file_path, biosample)
    
    # Add experiment details
    summary <- summary %>%
      mutate(Protein = protein, ExperimentID = experiment_id, Biosample = biosample)
    
    # Append to the list
    all_experiment_summaries[[paste0(protein, "_", experiment_id)]] <- summary
  }
}

# Combine all summaries into a single data frame
final_summary <- bind_rows(all_experiment_summaries)

# Ensure Chromatin_State is ordered
final_summary <- final_summary %>%
  mutate(
    Chromatin_State = factor(
      Chromatin_State,
      levels = c(names(chromatin_state_colors_short), 
                 "No motif in peak", "No CG in motif", "No State assignment for CG")
    )
  )

# Save the summary to a CSV file
output_csv <- file.path(out_dir, "chromatin_state_summary.csv")
write_csv(final_summary, output_csv)

###############################################################################
# detailed unknown cide 

library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)

# Path to the saved file
output_csv <- file.path(out_dir, "chromatin_state_summary.csv")

# Load the CSV file into a dataframe
final_summary <- read_csv(output_csv)

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
  "#FF69B4",  # No motif in peak -> Pink
  "#FF69B4",  # No CG in motif -> Dodger Blue
  "#FF69B4"   # No State assignment for CG -> Lime Green
)

chromatin_state_labels <- c(
  "1_TssA", "2_TssFlnk", "3_TssFlnkU", "4_TssFlnkD",
  "5_Tx", "6_TxWk", "7_EnhG1", "8_EnhG2",
  "9_EnhA1", "10_EnhA2", "11_EnhWk", "12_ZNF_Rpts",
  "13_Het", "14_TssBiv", "15_EnhBiv", "16_ReprPC",
  "17_ReprPCWk", "18_Quies", 
  "No motif in peak", "No CG in motif", "No State assignment for CG"
)

# Ensure Chromatin_State is ordered
final_summary <- final_summary %>%
  mutate(
    Chromatin_State = factor(
      Chromatin_State,
      levels = c(as.character(1:18), "No motif in peak", "No CG in motif", "No State assignment for CG")
    )
  )

# Calculate z-scores and identify outliers
outlier_summary <- final_summary %>%
  group_by(Chromatin_State) %>%
  mutate(
    Mean_Percentage = mean(Percentage),
    SD_Percentage = sd(Percentage),
    Z_Score = (Percentage - Mean_Percentage) / SD_Percentage,
    Is_Outlier = abs(Z_Score) > 3,  # Outlier if z-score > threshold
    Label = ifelse(Is_Outlier, paste(Protein, Biosample, sep = " "), NA)
  ) %>%
  ungroup()

# Create final plot with outlier labels using ggrepel
ggplot(outlier_summary, aes(x = Chromatin_State, y = Percentage, color = Chromatin_State)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +  # Jitter points
  geom_text_repel(
    aes(label = Label), 
    size = 3, color = "black", na.rm = TRUE, 
    box.padding = 0.3, point.padding = 0.2, max.overlaps = Inf
  ) +  # Repelling labels for outliers
  labs(
    title = "Chromatin State Percentage Distribution with Outliers",
    x = "Chromatin State",
    y = "Percentage"
  ) +
  scale_x_discrete(labels = chromatin_state_labels) +  # Set discrete x-axis labels with full names
  scale_color_manual(values = chromatin_state_colors) +  # Apply colors to Chromatin States
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "none"  # Drop legend
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1))  # Format y-axis as percentage

################################################################################
# library(ggplot2)
# library(ggrepel)
# library(dplyr)
# library(readr)
# 
# # Path to the saved file
# output_csv <- file.path(out_dir, "chromatin_state_summary.csv")
# 
# # Load the CSV file into a dataframe
# final_summary <- read_csv(output_csv)
# 
# # Define chromatin state names and colors
# chromatin_state_colors <- c(
#   "#FF0000",  # 1_TssA       -> Red
#   "#FF4500",  # 2_TssFlnk    -> Orange-Red
#   "#FF9900",  # 3_TssFlnkU   -> Orange
#   "#FFCC00",  # 4_TssFlnkD   -> Yellow-Orange
#   "#00CC00",  # 5_Tx         -> Green
#   "#006400",  # 6_TxWk       -> Dark Green
#   "#FFD700",  # 7_EnhG1      -> Gold
#   "#FFD700",  # 8_EnhG2      -> Gold (same as EnhG1)
#   "#FFFF00",  # 9_EnhA1      -> Yellow
#   "#FFDD00",  # 10_EnhA2     -> Yellow-Orange
#   "#FFEA73",  # 11_EnhWk     -> Light Yellow
#   "#9370DB",  # 12_ZNF_Rpts  -> Purple
#   "#C0C0C0",  # 13_Het       -> Light Gray
#   "#FF4500",  # 14_TssBiv    -> Orange-Red (similar to TssFlnk)
#   "#FFDD00",  # 15_EnhBiv    -> Yellow-Orange (similar to EnhA2)
#   "#808080",  # 16_ReprPC    -> Dark Gray
#   "#A9A9A9",  # 17_ReprPCWk  -> Light Gray
#   "#000000",  # 18_Quies     -> Black
#   "#808080"   # Unknown      -> Light Gray (default for unknown)
# )
# 
# chromatin_state_labels <- c(
#   "1_TssA", "2_TssFlnk", "3_TssFlnkU", "4_TssFlnkD",
#   "5_Tx", "6_TxWk", "7_EnhG1", "8_EnhG2",
#   "9_EnhA1", "10_EnhA2", "11_EnhWk", "12_ZNF_Rpts",
#   "13_Het", "14_TssBiv", "15_EnhBiv", "16_ReprPC",
#   "17_ReprPCWk", "18_Quies", "Unknown"
# )
# 
# # Ensure Chromatin_State is ordered
# final_summary <- final_summary %>%
#   mutate(
#     Chromatin_State = factor(
#       Chromatin_State,
#       levels = c(as.character(1:18), "Unknown")
#     )
#   )
# 
# # Calculate z-scores and identify outliers
# outlier_summary <- final_summary %>%
#   group_by(Chromatin_State) %>%
#   mutate(
#     Mean_Percentage = mean(Percentage),
#     SD_Percentage = sd(Percentage),
#     Z_Score = (Percentage - Mean_Percentage) / SD_Percentage,
#     Is_Outlier = abs(Z_Score) > 3,  # Outlier if z-score > threshold
#     Label = ifelse(Is_Outlier, paste(Protein, Biosample, sep = " "), NA)
#   ) %>%
#   ungroup()
# 
# # Create final plot with outlier labels using ggrepel
# ggplot(outlier_summary, aes(x = Chromatin_State, y = Percentage, color = Chromatin_State)) +
#   geom_jitter(width = 0.2, size = 2, alpha = 0.7) +  # Jitter points
#   geom_text_repel(
#     aes(label = Label), 
#     size = 3, color = "black", na.rm = TRUE, 
#     box.padding = 0.3, point.padding = 0.2, max.overlaps = Inf
#   ) +  # Repelling labels for outliers
#   labs(
#     title = "Chromatin State Percentage Distribution with Outliers",
#     x = "Chromatin State",
#     y = "Percentage"
#   ) +
#   scale_x_discrete(labels = chromatin_state_labels) +  # Set discrete x-axis labels with full names
#   scale_color_manual(values = chromatin_state_colors) +  # Apply colors to Chromatin States
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
#     legend.position = "none"  # Drop legend
#   ) +
#   scale_y_continuous(labels = scales::percent_format(scale = 1))  # Format y-axis as percentage
###################################### sibnglke protein ###########################
# Define output directory for plots
protein_plot_dir <- file.path(output_dir_plots, "protein_plots")

# Get unique proteins
proteins <- unique(final_summary$Protein)

# Loop through each protein and generate a plot
for (protein in proteins) {
  
  # Filter data for the current protein
  protein_data <- final_summary %>%
    mutate(Label = ifelse(Protein == protein, paste(Protein, Biosample, sep = " "), NA))  # Label points for this protein
  
  # Create a base jittered plot
  jitter_plot <- ggplot(protein_data, aes(x = Chromatin_State, y = Percentage, color = Chromatin_State)) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.7)
  
  # Extract jittered positions
  jitter_data <- ggplot_build(jitter_plot)$data[[1]]
  
  # Combine jittered positions with the original data
  protein_data <- protein_data %>%
    mutate(
      Jittered_X = jitter_data$x,
      Jittered_Y = jitter_data$y
    )
  
  # Create final plot with lines always connecting to the center of the point
  final_plot <- ggplot(protein_data, aes(x = Jittered_X, y = Jittered_Y, color = Chromatin_State)) +
    geom_point(size = 2, alpha = 0.7) +  # Use jittered positions
    geom_segment(
      aes(x = Jittered_X, xend = Jittered_X, y = Jittered_Y, yend = Jittered_Y),
      linetype = "solid", color = "black", alpha = 0.5
    ) +  # Lines connecting to the center of the points
    geom_text_repel(
      aes(label = Label),
      size = 3, color = "black", na.rm = TRUE,
      box.padding = 0.3, point.padding = 0.2, max.overlaps = Inf
    ) +  # Repelling labels for protein-specific points
    labs(
      title = paste("Chromatin State Distribution for Protein:", protein),
      x = "Chromatin State",
      y = "Percentage"
    ) +
    scale_x_continuous(breaks = 1:19, labels = chromatin_state_labels) +  # Set discrete x-axis labels with full names
    scale_color_manual(values = chromatin_state_colors) +  # Apply colors to Chromatin States
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
      legend.position = "none"  # Drop legend
    ) +
    scale_y_continuous(labels = scales::percent_format(scale = 1))  # Format y-axis as percentage
  
  # Save the plot
  plot_file <- file.path(protein_plot_dir, paste0(protein, "_chromatin_state_plot.png"))
  ggsave(plot_file, final_plot, width = 24, height = 12)
  cat("Plot saved for protein:", protein, "\n")
  
}
