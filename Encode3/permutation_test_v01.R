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
  "meme",
  "fimo_single_experiments_on_known_motif_and_peaks_add_methylation_v3"
)
output_folder <- file.path(this.dir, "permutation_test")

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

# Helper function to extract column names from BED files
extract_colnames <- function(file_path) {
  first_row <- readLines(file_path, n = 1)
  cleaned_row <- gsub("#", "", first_row)  # Remove leading #
  colnames <- unlist(strsplit(cleaned_row, "\t"))
  return(colnames)
}

# Function to read BED files with dynamically extracted column names
read_bed_file <- function(file_path) {
  colnames <- extract_colnames(file_path)
  df <- readr::read_delim(
    file_path,
    delim = "\t",
    skip = 1,
    col_names = colnames,
    show_col_types = FALSE
  )
  return(df)
}

# List all protein folders, excluding irrelevant directories
protein_folders <- list.dirs(path = input_folder,
                             recursive = FALSE,
                             full.names = FALSE)
protein_folders <- protein_folders[!grepl("sbatch_scripts|log", protein_folders)]

################################# Main Process ##################################
# Check if `summary_df.rds` exists; if not, generate the summary
summary_df_path <- file.path(output_folder, "summary_df.rds")
if (!file.exists(summary_df_path)) {
  # Initialize an empty list to store summary results
  summary_list <- list()
  
  # Loop over each protein folder
  for (protein in protein_folders) {
    print(paste("Processing protein:", protein))
    
    # Get all files for the current protein
    files <- list.files(path = file.path(input_folder, protein))
    
    df_protein <- data.frame()
    
    # Process each BED file
    for (bed_file in files) {
      print(paste("Processing file:", bed_file))
      
      # Extract metadata from file name
      parts <- strsplit(bed_file, "_")[[1]]
      biosample <- parts[1]  # Example: GM12878
      experiment_id <- parts[3]  # Example: ENCFF320KXO
      motif <- sub("\\.bed$", "", parts[4])  # Extract motif
      
      # Read the BED file
      df <- read_bed_file(file_path = file.path(input_folder, protein, bed_file))
      df_protein <- rbind(df_protein, df)
      
      # Create a data frame with only unique peaks
      df_unique_peaks <- df %>% distinct(chr, start_peak, end_peak, .keep_all = TRUE)
      
      # Calculate categories based on unique peaks
      n_peak_with_CpG <- sum(df_unique_peaks$start_motif != -1 &
                               df_unique_peaks$start_cg != 0)
      n_peak_no_CpG <- sum(df_unique_peaks$start_motif != -1 &
                             df_unique_peaks$start_cg == 0)
      n_peak_no_motif_no_CpG <- sum(df_unique_peaks$start_motif == -1 &
                                      df_unique_peaks$start_cg == 0)
      n_peaks <- nrow(df_unique_peaks)
      
      # Append results to summary list
      summary_list[[length(summary_list) + 1]] <- data.frame(
        protein = protein,
        biosample = biosample,
        experiment_id = experiment_id,
        motif = motif,
        n_peaks = n_peaks,
        n_peak_with_CpG = n_peak_with_CpG,
        n_peak_no_CpG = n_peak_no_CpG,
        n_peak_no_motif_no_CpG = n_peak_no_motif_no_CpG,
        n_CpG = sum(df$start_cg != -1)
      )
    }
    protein_output_dir <- file.path(output_folder, protein)
    if (!dir.exists(protein_output_dir)) {
      dir.create(protein_output_dir, recursive = TRUE)
    }
    
    # Save combined data for the protein
    saveRDS(object = df_protein,
            file = file.path(output_folder, protein, paste0(protein, ".rds")))
  }
  
  # Combine all summaries into a single data frame
  summary_df <- do.call(rbind, summary_list)
  summary_df$no_motif_no_CpG_ratio <- summary_df$n_peak_no_motif_no_CpG / summary_df$n_peaks
  summary_df$n_motif <- (summary_df$n_peak_no_CpG +  summary_df$n_peak_with_CpG) / summary_df$n_peaks
  
  # Save the summary data frame
  saveRDS(object = summary_df, file = summary_df_path)
  print(paste("Summary data saved at:", summary_df_path))
} else {
  # Load existing summary data
  summary_df <- readRDS(summary_df_path)
  print(paste("Loaded existing summary data from:", summary_df_path))
}

################################# End Script ####################################
end_script <- Sys.time()
print(paste("Script completed in", round(
  difftime(end_script, start_script, units = "mins"), 2
), "minutes"))




# # Plot the histogram
# histogram_plot <- ggplot(summary_df, aes(x = no_motif_no_CpG_ratio)) +
#   geom_histogram(
#     binwidth = 0.05,
#     fill = "#00BFC4",
#     color = "black",
#     alpha = 0.7
#   ) +
#   labs(title = "Distribution of No Motif No CpG Peaks Ratio", x = "Ratio of No Motif No CpG Peaks", y = "Count") +
#   theme_minimal()
#
# # Save the absolute plot
# ggsave(
#   filename = file.path(output_folder, "histogram_plot_fractions.png"),
#   plot = histogram_plot,
#   width = 10,
#   height = 6
# )
#
# # Loop over protein folders
# for (protein in protein_folders) {
#   # protein <- "CTCF"
#   print(protein)
#
#   # Create output directory for each protein if it doesn't exist
#   protein_output_dir <- file.path(output_folder, protein)
#   if (!dir.exists(protein_output_dir)) {
#     dir.create(protein_output_dir, recursive = TRUE)
#   }
#
#   # Filter the data for the specific protein
#
#   protein_data <- summary_df[summary_df$protein == protein, ]
#
#   # # Create the bar plot
#   # bar_plot <- ggplot(protein_data, aes(x = interaction(biosample, experiment_id, motif),
#   #                                      y = no_motif_no_CpG_ratio,
#   #                                      fill = biosample)) +
#   #   geom_bar(stat = "identity") +
#   #   labs(
#   #     title = paste("Ratio of Peaks with No Motif and No CpG for", protein),
#   #     x = "Biosample x Experiment ID",
#   #     y = "Ratio of No Motif No CpG Peaks"
#   #   ) +
#   #   scale_fill_manual(values = group.colors) +  # Use the custom biosample colors
#   #   theme_minimal() +
#   #   theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
#
#   # Print the plot
#   print(bar_plot)
#
#   ########### continue
#
#   # Now run the code
#   # summary_df_long <- summary_df %>%
#   #   mutate(matched_peaks = n_peaks - n_peaks_unmatched) %>%  # Calculate matched peaks
#   #   select(protein, biosample, experiment_id, motif, n_peaks, n_peaks_unmatched, matched_peaks) %>%
#   #   pivot_longer(cols = c(matched_peaks, n_peaks_unmatched),
#   #                names_to = "peak_type",
#   #                values_to = "peak_count")  # Reshape data to long format for stacked bars
#   #
#   # # Create a new variable that combines experiment_id and motif
#   # summary_df_long <- summary_df_long %>%
#   #   mutate(experiment_motif = paste(biosample,experiment_id, motif, sep = "_"))
#   #
#   #
#   # # Generate plot for absolute number of peaks
#   # p_absolute <- ggplot(data = summary_df_long[summary_df_long$protein == protein,],
#   #                      aes(x = experiment_motif, y = peak_count, fill = biosample, alpha = peak_type)) +
#   #   geom_bar(stat = "identity", position = "stack") +  # Stack matched and unmatched peaks
#   #   scale_alpha_manual(values = c(matched_peaks = 1, n_peaks_unmatched = 0.5)) +  # Different alpha for unmatched peaks
#   #   scale_fill_manual(values = group.colors) +  # Use the custom color palette
#   #   labs(title = paste("Number of Peaks per Biosample for", protein),
#   #        x = "Experiment ID and Motif",
#   #        y = "Number of Peaks") +
#   #   theme_minimal() +
#   #   theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for readability
#   #
#   # # Save the absolute plot
#   # ggsave(filename = file.path(protein_output_dir, paste0(protein, "_absolute_peaks.png")),
#   #        plot = p_absolute, width = 10, height = 6)
#   #
#   # # Calculate percentage for each row in summary_df_long
#   # summary_df_long <- summary_df_long %>%
#   #   group_by(experiment_motif, biosample) %>%  # Group by experiment_motif and biosample
#   #   mutate(percentage = peak_count / sum(peak_count) * 100) %>%  # Calculate percentage
#   #   ungroup()
#   #
#   # # Generate plot for percentage of peaks
#   # p_percentage <- ggplot(data = summary_df_long[summary_df_long$protein == protein,],
#   #                        aes(x = experiment_motif, y = percentage, fill = biosample, alpha = peak_type)) +
#   #   geom_bar(stat = "identity", position = "stack") +  # Stack matched and unmatched peaks
#   #   scale_alpha_manual(values = c(matched_peaks = 1, n_peaks_unmatched = 0.5)) +  # Different alpha for unmatched peaks
#   #   scale_fill_manual(values = group.colors) +  # Use the custom color palette
#   #   labs(title = paste("Percentage of Peaks per Biosample for", protein),
#   #        x = "Experiment ID and Motif",
#   #        y = "Percentage of Peaks") +
#   #   theme_minimal() +
#   #   theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for readability
#   #
#   # # Save the percentage plot
#   # ggsave(filename = file.path(protein_output_dir, paste0(protein, "_percentage_peaks.png")),
#   #        plot = p_percentage, width = 10, height = 6)
# }

# # Add the three ratios as new columns in summary_df
# summary_df <- summary_df %>%
#   mutate(ratio_CpG_to_peaks = n_CpG / n_peaks,  # First ratio
#          ratio_CpG_to_matched_peaks = n_CpG / (n_peaks - n_peaks_unmatched),  # Second ratio
#          ratio_matched_peaks_to_total = (n_peaks - n_peaks_unmatched) / n_peaks)  # New ratio
#
# # Calculate the maximum value for each ratio column
# max_ratio1 <- max(summary_df$ratio_CpG_to_peaks, na.rm = TRUE)
# max_ratio2 <- max(summary_df$ratio_CpG_to_matched_peaks, na.rm = TRUE)
# max_ratio3 <- max(summary_df$ratio_matched_peaks_to_total, na.rm = TRUE)  # New max value
#
# # Plot histogram for the first ratio (n_CpG / n_peaks)
# p_histogram_ratio1 <- ggplot(summary_df, aes(x = ratio_CpG_to_peaks)) +
#   geom_histogram(binwidth = 0.05, alpha = 0.7, color = "black") +
#   labs(title = "Histogram of n_CpG / n_peaks",
#        x = "Ratio of CpG to Peaks",
#        y = "Frequency") +
#   scale_x_continuous(breaks = seq(0, ceiling(max_ratio1), by = 1),  # Major ticks at whole numbers
#                      minor_breaks = seq(0, ceiling(max_ratio1), by = 0.5)) +  # Minor ticks at 0.5 intervals
#   theme_minimal()
#
# # Plot histogram for the second ratio (n_CpG / (n_peaks - n_peaks_unmatched))
# p_histogram_ratio2 <- ggplot(summary_df, aes(x = ratio_CpG_to_matched_peaks)) +
#   geom_histogram(binwidth = 0.05, alpha = 0.7, color = "black") +
#   labs(title = "Histogram of n_CpG / (n_peaks - n_peaks_unmatched)",
#        x = "Ratio of CpG to Matched Peaks",
#        y = "Frequency") +
#   scale_x_continuous(breaks = seq(0, ceiling(max_ratio2), by = 1),  # Major ticks at whole numbers
#                      minor_breaks = seq(0, ceiling(max_ratio2), by = 0.5)) +  # Minor ticks at 0.5 intervals
#   theme_minimal()
#
# # Plot histogram for the new ratio ((n_peaks - n_peaks_unmatched) / n_peaks)
# p_histogram_ratio3 <- ggplot(summary_df, aes(x = ratio_matched_peaks_to_total)) +
#   geom_histogram(binwidth = 0.01, alpha = 0.7, color = "black") +
#   labs(title = "Histogram of Matched Peaks / Total Peaks",
#        x = "Ratio of Matched Peaks to Total Peaks",
#        y = "Frequency") +
#   scale_x_continuous(breaks = seq(0, 1, by = 0.1),  # Major ticks at whole numbers
#                      minor_breaks = seq(0, 1, by = 0.01)) +  # Minor ticks at 0.5 intervals
#   theme_minimal()
#
# # Display all three histograms
# print(p_histogram_ratio1)
# print(p_histogram_ratio2)
# print(p_histogram_ratio3)
#
#
# # Add the ratios and a new column that categorizes entries based on ratio_matched_peaks_to_total
# summary_df <- summary_df %>%
#   mutate(ratio_CpG_to_peaks = n_CpG / n_peaks,  # First ratio
#          ratio_CpG_to_matched_peaks = n_CpG / (n_peaks - n_peaks_unmatched),  # Second ratio
#          ratio_matched_peaks_to_total = (n_peaks - n_peaks_unmatched) / n_peaks,  # Third ratio
#          low_matched_ratio = ifelse(ratio_matched_peaks_to_total < 0.25, "< 0.25", ">= 0.25"))  # New category
#
# # Calculate the maximum value for each ratio column
# max_ratio1 <- max(summary_df$ratio_CpG_to_peaks, na.rm = TRUE)
# max_ratio2 <- max(summary_df$ratio_CpG_to_matched_peaks, na.rm = TRUE)
# max_ratio3 <- max(summary_df$ratio_matched_peaks_to_total, na.rm = TRUE)
#
# # Plot histogram for the first ratio (n_CpG / n_peaks) with stacked colors
# p_histogram_ratio1 <- ggplot(summary_df, aes(x = ratio_CpG_to_peaks, fill = low_matched_ratio)) +
#   geom_histogram(binwidth = 0.05, alpha = 0.7, color = "black", position = "stack") +
#   scale_fill_manual(values = c("< 0.25" = "red", ">= 0.25" = "blue")) +
#   labs(title = "Histogram of n_CpG / n_peaks",
#        x = "Ratio of CpG to Peaks",
#        y = "Frequency") +
#   scale_x_continuous(breaks = seq(0, ceiling(max_ratio1), by = 1),  # Major ticks at whole numbers
#                      minor_breaks = seq(0, ceiling(max_ratio1), by = 0.5)) +  # Minor ticks at 0.5 intervals
#   theme_minimal()
#
# # Plot histogram for the second ratio (n_CpG / (n_peaks - n_peaks_unmatched)) with stacked colors
# p_histogram_ratio2 <- ggplot(summary_df, aes(x = ratio_CpG_to_matched_peaks, fill = low_matched_ratio)) +
#   geom_histogram(binwidth = 0.05, alpha = 0.7, color = "black", position = "stack") +
#   scale_fill_manual(values = c("< 0.25" = "red", ">= 0.25" = "blue")) +
#   labs(title = "Histogram of n_CpG / (n_peaks - n_peaks_unmatched)",
#        x = "Ratio of CpG to Matched Peaks",
#        y = "Frequency") +
#   scale_x_continuous(breaks = seq(0, ceiling(max_ratio2), by = 1),  # Major ticks at whole numbers
#                      minor_breaks = seq(0, ceiling(max_ratio2), by = 0.5)) +  # Minor ticks at 0.5 intervals
#   theme_minimal()
#
# # Plot histogram for the new ratio ((n_peaks - n_peaks_unmatched) / n_peaks) with stacked colors
# p_histogram_ratio3 <- ggplot(summary_df, aes(x = ratio_matched_peaks_to_total, fill = low_matched_ratio)) +
#   geom_histogram(binwidth = 0.01, alpha = 0.7, color = "black", position = "stack") +
#   scale_fill_manual(values = c("< 0.25" = "red", ">= 0.25" = "blue")) +
#   labs(title = "Histogram of Matched Peaks / Total Peaks",
#        x = "Ratio of Matched Peaks to Total Peaks",
#        y = "Frequency") +
#   scale_x_continuous(breaks = seq(0, 1, by = 0.1),  # Major ticks at whole numbers
#                      minor_breaks = seq(0, 1, by = 0.01)) +  # Minor ticks at 0.5 intervals
#   theme_minimal()
#
# # Display all three histograms
# print(p_histogram_ratio1)
# print(p_histogram_ratio2)
# print(p_histogram_ratio3)

##########################################################
# library(ComplexUpset)
# library(ggplot2)
# library(reshape2)
#
# # Function to create intersection bar plot for a motif
# create_intersection_barplot <- function(cpg_list, experiment_biosample, motif, output_folder, protein) {
#
#   # Convert the CpG list into a binary matrix for ComplexUpset
#   binary_cpg_df <- fromList(cpg_list)  # List of CpG positions for each experiment
#
#   # Melt binary matrix to long format (for easier ggplot use)
#   binary_cpg_long <- melt(as.matrix(binary_cpg_df), varnames = c("CpG_Position", "Experiment"), value.name = "Presence")
#
#   # Filter only rows where the CpG is present
#   binary_cpg_long <- binary_cpg_long %>% filter(Presence == 1)
#
#   # Create a column representing the intersection combinations as a string (e.g., "Exp1_Exp2")
#   binary_cpg_long <- binary_cpg_long %>%
#     group_by(CpG_Position) %>%
#     summarize(Intersection = paste(Experiment, collapse = "_"))  # Collapses experiments involved in each CpG
#
#   # Count how many CpGs are in each intersection combination, split by experiment
#   intersection_counts <- binary_cpg_long %>%
#     separate(Intersection, into = paste0("Exp", 1:ncol(binary_cpg_df)), sep = "_", fill = "right", remove = FALSE) %>%
#     pivot_longer(cols = starts_with("Exp"), names_to = "experiment", values_to = "Experiment", values_drop_na = TRUE) %>%
#     group_by(Experiment, Intersection) %>%
#     summarise(CpG_Count = n())
#
#   # Add a biosample column based on experiment names
#   intersection_counts$biosample <- experiment_biosample$biosample[match(intersection_counts$Experiment, experiment_biosample$experiment_id)]
#
#   # Map colors for the intersections based on biosample and combinations
#   unique_intersections <- unique(intersection_counts$Intersection)
#   intersection_colors <- rainbow(length(unique_intersections))  # Customize this palette if needed
#
#   # Create the stacked bar plot for intersections
#   p <- ggplot(intersection_counts, aes(x = Experiment, y = CpG_Count, fill = Intersection)) +
#     geom_bar(stat = "identity", position = "stack") +
#     scale_fill_manual(values = intersection_colors) +  # Apply colors based on intersection
#     labs(title = paste("CpG Counts for Motif:", motif),
#          x = "Experiment",
#          y = "CpG Count") +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels
#
#   # Correct the file path for saving the plot
#   output_file <- file.path(output_folder, protein, paste0("IntersectionBarplot_", motif, ".png"))
#   if (!dir.exists(file.path(output_folder, protein))) {
#     dir.create(file.path(output_folder, protein), recursive = TRUE)  # Ensure the directory exists
#   }
#
#   # Save the plot to a file
#   ggsave(output_file, plot = p, width = 10, height = 6)
#
#   print(paste("Bar plot saved for motif:", motif, "in protein folder:", protein))
# }
#
# # Iterate over each motif and create a bar plot for it
# for (motif in unique_motifs) {
#   df_motif <- df_protein %>% filter(motif == !!motif)
#
#   # Create the CpG list (same as for Venn diagram)
#   cpg_list <- list()
#   for (exp_id in unique(df_motif$experiment_id)) {
#     cpg_positions <- df_motif %>%
#       filter(experiment_id == !!exp_id) %>%
#       select(chr, start, end) %>%
#       mutate(position = paste(chr, start, end, sep = ":")) %>%
#       pull(position)
#     cpg_list[[exp_id]] <- unique(cpg_positions)
#   }
#
#   # Get biosample information for experiments
#   experiment_biosample <- df_motif %>%
#     select(experiment_id, biosample) %>%
#     distinct()
#
#   # Call the function to create the bar plot
#   create_intersection_barplot(cpg_list, experiment_biosample, motif, output_folder, protein)
# }

# ############################################################################
# library(ComplexUpset)
# library(ggplot2)
# library(reshape2)
# library(colorspace)  # For lightening the colors
#
# # Function to create intersection bar plot for a motif with dynamic color scaling
# create_intersection_barplot <- function(cpg_list, experiment_biosample, motif, output_folder, protein, group_colors) {
#
#   # Convert the CpG list into a binary matrix for ComplexUpset
#   binary_cpg_df <- fromList(cpg_list)  # List of CpG positions for each experiment
#
#   # Melt binary matrix to long format (for easier ggplot use)
#   binary_cpg_long <- melt(as.matrix(binary_cpg_df), varnames = c("CpG_Position", "Experiment"), value.name = "Presence")
#
#   # Filter only rows where the CpG is present
#   binary_cpg_long <- binary_cpg_long %>% filter(Presence == 1)
#
#   # Create a column representing the intersection combinations as a string (e.g., "Exp1_Exp2")
#   binary_cpg_long <- binary_cpg_long %>%
#     group_by(CpG_Position) %>%
#     summarize(Intersection = paste(Experiment, collapse = "_"))  # Collapses experiments involved in each CpG
#
#   # Count how many CpGs are in each intersection combination, split by experiment
#   intersection_counts <- binary_cpg_long %>%
#     separate(Intersection, into = paste0("Exp", 1:ncol(binary_cpg_df)), sep = "_", fill = "right", remove = FALSE) %>%
#     pivot_longer(cols = starts_with("Exp"), names_to = "experiment", values_to = "Experiment", values_drop_na = TRUE) %>%
#     group_by(Experiment, Intersection) %>%
#     summarise(CpG_Count = n(), .groups = 'drop')
#
#   # Add a biosample column based on experiment names
#   intersection_counts$biosample <- experiment_biosample$biosample[match(intersection_counts$Experiment, experiment_biosample$experiment_id)]
#
#   # Determine whether each intersection involves multiple experiments or just one
#   intersection_counts <- intersection_counts %>%
#     group_by(Intersection) %>%
#     mutate(experiment_count = n_distinct(Experiment)) %>%  # Count how many experiments are involved in each intersection
#     ungroup()
#
#   # Assign colors based on biosample and experiment count (light for single, dark for multiple)
#   intersection_counts$color <- sapply(1:nrow(intersection_counts), function(i) {
#     base_color <- group_colors[intersection_counts$biosample[i]]
#     if (intersection_counts$experiment_count[i] > 1) {
#       base_color  # Darker version for multiple experiments
#     } else {
#       lighten(base_color, amount = 0.5)  # Lighter version for single experiments using colorspace::lighten()
#     }
#   })
#
#   # Create the stacked bar plot for intersections
#   p <- ggplot(intersection_counts, aes(x = Experiment, y = CpG_Count, fill = Intersection)) +
#     geom_bar(stat = "identity", position = "stack", aes(fill = color)) +
#     scale_fill_identity() +  # Use the assigned colors
#     labs(title = paste("CpG Counts for Motif:", motif),
#          x = "Experiment",
#          y = "CpG Count") +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels
#
#   # Correct the file path for saving the plot
#   output_file <- file.path(output_folder, protein, paste0("IntersectionBarplot_", motif, ".png"))
#   if (!dir.exists(file.path(output_folder, protein))) {
#     dir.create(file.path(output_folder, protein), recursive = TRUE)  # Ensure the directory exists
#   }
#
#   # Save the plot to a file
#   ggsave(output_file, plot = p, width = 10, height = 6)
#
#   print(paste("Bar plot saved for motif:", motif, "in protein folder:", protein))
# }
#
# # Iterate over each motif and create a bar plot for it
# for (motif in unique_motifs) {
#   df_motif <- df_protein %>% filter(motif == !!motif)
#
#   # Create the CpG list (same as for Venn diagram)
#   cpg_list <- list()
#   for (exp_id in unique(df_motif$experiment_id)) {
#     cpg_positions <- df_motif %>%
#       filter(experiment_id == !!exp_id) %>%
#       select(chr, start, end) %>%
#       mutate(position = paste(chr, start, end, sep = ":")) %>%
#       pull(position)
#     cpg_list[[exp_id]] <- unique(cpg_positions)
#   }
#
#   # Get biosample information for experiments
#   experiment_biosample <- df_motif %>%
#     select(experiment_id, biosample) %>%
#     distinct()
#
#   # Call the function to create the bar plot
#   create_intersection_barplot(cpg_list, experiment_biosample, motif, output_folder, protein, group.colors)
# }
#}

########################## copy data oif pairs of sampples #####################

if (!exists("summary_df")) {
  summary_df <- readRDS(file = file.path(output_folder, "summary_df.rds"))
  message("Loaded summary_df from file.")
} else {
  message("summary_df is already loaded in the environment.")
}



#Filter summary_df to keep only the rows with the lowest no_motif_no_CpG_ratio for each (protein, biosample, motif) combination
filtered_summary_df <- summary_df %>%
  group_by(protein, biosample, motif) %>%
  slice_min(no_motif_no_CpG_ratio, with_ties = FALSE) %>%
  ungroup()


# Loop over protein folders
for (protein in unique(summary_df$protein)) {
  # protein <- "CTCF"
  print(protein)
  
  # Create output directory for each protein if it doesn't exist
  protein_output_dir <- file.path(output_folder, protein)
  if (!dir.exists(protein_output_dir)) {
    dir.create(protein_output_dir, recursive = TRUE)
  }
  protein_data <- filtered_summary_df[filtered_summary_df$protein == protein, ]
  # Filter the data for the specific protein
  
  # Loop over unique motifs
  unique_motifs <- unique(protein_data$motif)
  
  for (motif in unique_motifs) {
    # Filter the data for the current motif
    # debug # motif <- unique_motifs[1]
    print(motif)
    motif_data <- protein_data[protein_data$motif == motif, ]
    
    # Generate all combinations of two biosamples
    biosample_combinations <- combn(unique(motif_data$biosample), 2, simplify = FALSE)
    
    # Loop through each combination of two biosamples
    for (biosample_comb in biosample_combinations) {
      print(biosample_comb)
      
      merge_df_file_name <- file.path(
        output_folder,
        protein,
        paste(
          unlist(biosample_comb)[1],
          unlist(biosample_comb)[2],
          protein,
          motif,
          "best_comparison_all_data_points.rds",
          sep = "_"
        )
      )
      if (!exists(merge_df_file_name)) {
        next
      }
      
      # Read the BED file
      df1 <- read_bed_file(file_path = file.path(
        input_folder,
        protein,
        paste(
          biosample_comb[1],
          protein,
          motif_data$experiment_id[motif_data$biosample == biosample_comb[1]],
          motif,
          "ChIP_with_methylation.bed",
          sep = "_"
        )
      ))
      
      # Read the BED file
      df2 <- read_bed_file(file_path = file.path(
        input_folder,
        protein,
        paste(
          biosample_comb[2],
          protein,
          motif_data$experiment_id[motif_data$biosample == biosample_comb[2]],
          motif,
          "ChIP_with_methylation.bed",
          sep = "_"
        )
      ))
      
      # Discard unmatched values
      df1 <- df1[df1$start_cg != 0, ]
      df2 <- df2[df2$start_cg != 0, ]
      
      # Convert all "." values to NA in df1 and df2
      df1 <- df1 %>%
        mutate_all(~ replace(., . == ".", NA))
      
      df2 <- df2 %>%
        mutate_all(~ replace(., . == ".", NA))
      
      # List of columns to be converted to integers where applicable for both df1 and df2
      integer_columns <- c(
        "start_cg",
        "end_cg",
        "start_motif",
        "end_motif",
        "start_peak",
        "end_peak",
        "mRead_A549",
        "nRead_A549",
        "mRead_GM12878",
        "nRead_GM12878",
        "mRead_HepG2",
        "nRead_HepG2",
        "mRead_K562",
        "nRead_K562"
      )
      
      # Convert relevant columns to integers for df1
      df1[integer_columns] <- lapply(df1[integer_columns], function(x)
        as.integer(as.character(x)))
      
      # Convert relevant columns to integers for df2
      df2[integer_columns] <- lapply(df2[integer_columns], function(x)
        as.integer(as.character(x)))
      
      # Remove duplicates from both df1 and df2 before the join based on chr, start_cg, end_cg, strand
      df1 <- df1 %>%
        distinct(chr, start_cg, end_cg, strand, .keep_all = TRUE)
      
      df2 <- df2 %>%
        distinct(chr, start_cg, end_cg, strand, .keep_all = TRUE)
      
      # Recalculate all fRead_* columns in df1
      df1 <- df1 %>%
        mutate(across(
          starts_with("fRead_"),
          .fns = function(x) {
            # Extract the base name from the current column (e.g., HepG2 from fRead_HepG2)
            base_col <- sub("fRead_", "", cur_column())
            mRead_col <- paste0("mRead_", base_col)
            nRead_col <- paste0("nRead_", base_col)
            
            # Perform the calculation with NA handling for division by zero
            ifelse(df1[[nRead_col]] == 0, NA, df1[[mRead_col]] / df1[[nRead_col]])
          },
          .names = "{col}"
        ))
      
      # Recalculate all fRead_* columns in df2
      df2 <- df2 %>%
        mutate(across(
          starts_with("fRead_"),
          .fns = function(x) {
            # Extract the base name from the current column (e.g., HepG2 from fRead_HepG2)
            base_col <- sub("fRead_", "", cur_column())
            mRead_col <- paste0("mRead_", base_col)
            nRead_col <- paste0("nRead_", base_col)
            
            # Perform the calculation with NA handling for division by zero
            ifelse(df2[[nRead_col]] == 0, NA, df2[[mRead_col]] / df2[[nRead_col]])
          },
          .names = "{col}"
        ))
      
      # First, drop all unneeded columns from df1 and df2, keeping only the relevant ones
      df1 <- df1 %>%
        select(
          chr,
          start_cg,
          end_cg,
          strand,
          sample,
          experiment,
          paste0("Chromatin_State_", biosample_comb[1]),
          paste0("fRead_", biosample_comb[1]),
          paste0("Chromatin_State_", biosample_comb[2]),
          paste0("fRead_", biosample_comb[2])
        )
      
      df2 <- df2 %>%
        select(
          chr,
          start_cg,
          end_cg,
          strand,
          sample,
          experiment,
          paste0("Chromatin_State_", biosample_comb[1]),
          paste0("fRead_", biosample_comb[1]),
          paste0("Chromatin_State_", biosample_comb[2]),
          paste0("fRead_", biosample_comb[2])
        )
      
      # Perform the full join on the relevant columns
      merged_df <- full_join(
        df1,
        df2,
        by = c(
          "chr",
          "start_cg",
          "end_cg",
          "strand",
          paste0("Chromatin_State_", biosample_comb[1]),
          paste0("fRead_", biosample_comb[1]),
          paste0("Chromatin_State_", biosample_comb[2]),
          paste0("fRead_", biosample_comb[2])
        ),
        #      suffix = c(paste0("_",biosample_comb[1]), paste0("_",biosample_comb[1]))
      )
      
      # Combine the relevant fields (sample, experiment, etc.)
      merged_df <- merged_df %>%
        mutate(
          # Combine samples: if both are present, concatenate them; if only one is present, keep the non-NA value
          sample = case_when(
            !is.na(sample.x) &
              !is.na(sample.y) ~ paste0(sample.x, "_", sample.y),
            # Both are present!is.na(sample.x) ~ sample.x,
            # Only sample.x is present!is.na(sample.y) ~ sample.y                                             # Only sample.y is present
          ),
          
          # Combine experiments: if both are present, concatenate them; if only one is present, keep the non-NA value
          experiment = case_when(
            !is.na(experiment.x) &
              !is.na(experiment.y) ~ paste0(experiment.x, "_", experiment.y),
            # Both are present!is.na(experiment.x) ~ experiment.x,
            # Only experiment.x is present!is.na(experiment.y) ~ experiment.y                                                    # Only experiment.y is present
          )
        ) %>%
        # Select the final set of relevant columns
        select(
          chr,
          start_cg,
          end_cg,
          strand,
          sample,
          experiment,
          paste0("Chromatin_State_", biosample_comb[1]),
          paste0("fRead_", biosample_comb[1]),
          paste0("Chromatin_State_", biosample_comb[2]),
          paste0("fRead_", biosample_comb[2])
        ) %>%
        distinct()  # Ensure no duplicates remain
      
      # Check the structure of the cleaned merged data
      #str(merged_df)
      
      saveRDS(object = merged_df, file =  merge_df_file_name)
      
    }
  }
  
}

############################## calculate statistric #############################

# Function to calculate S_statistic for a given data frame
calculate_S_statistic <- function(df, biosample_comb) {
  mean_diff_per_sample <- df %>%
    mutate(diff_fRead = !!sym(paste0("fRead_", biosample_comb[1])) -!!sym(paste0("fRead_", biosample_comb[2]))) %>%  # Calculate difference
    group_by(sample) %>%  # Group by sample
    summarize(mean_diff = mean(diff_fRead, na.rm = TRUE))  # Calculate mean difference per sample
  
  # Calculate the S statistic
  S_statistic <- mean_diff_per_sample$mean_diff[biosample_comb[1] == mean_diff_per_sample$sample] -
    mean_diff_per_sample$mean_diff[biosample_comb[2] == mean_diff_per_sample$sample]
  
  return(S_statistic)
}

library(doParallel)
library(foreach)

# Set up the parallel backend (before the outer loop)
num_cores <- detectCores() - 1  # Use one fewer core than available
cl <- makeCluster(num_cores)
registerDoParallel(cl)

summary_test <- list()
n_permutations <- 100000

start <- which("NR3C1" == protein_folders) - 1
stop <- which("NR3C1" == protein_folders) + 1
# [start:stop]
# Loop over protein folders
for (protein in protein_folders) {
  # [1:1]
  # protein <- "CTCF"
  print(protein)
  
  protein_plot_list <- list()
  
  protein_data <- filtered_summary_df[filtered_summary_df$protein == protein, ]
  # Filter the data for the specific protein
  
  # Loop over unique motifs
  unique_motifs <- unique(protein_data$motif)
  
  protein_output_dir <- file.path(output_folder, protein)
  
  for (motif in unique_motifs) {
    # Filter the data for the current motif
    # debug # motif <- unique_motifs[1]
    
    
    print(motif)
    motif_data <- protein_data[protein_data$motif == motif, ]
    
    # Generate all combinations of two biosamples
    biosample_combinations <- combn(unique(motif_data$biosample), 2, simplify = FALSE)
    
    motif_out_dir <- file.path(protein_output_dir, motif)
    
    # Loop through each combination of two biosamples
    for (biosample_comb in biosample_combinations) {
      print(biosample_comb)
      biosample_both <- paste(unlist(biosample_comb)[1], unlist(biosample_comb)[2], sep = "_")
      
      biosample_plot_list <- list()
      # Create output directory for each protein if it doesn't exist
      sim_output_dir <- file.path(motif_out_dir, biosample_both)
      if (!dir.exists(sim_output_dir)) {
        dir.create(sim_output_dir, recursive = TRUE)
      }
      
      # Dynamically generate the chromatin state column names for both biosamples
      chromatin_state_col1 <- paste0("Chromatin_State_", biosample_comb[1])
      chromatin_state_col2 <- paste0("Chromatin_State_", biosample_comb[2])
      
      
      merged_df <- readRDS(file = file.path(
        protein_output_dir,
        paste(
          unlist(biosample_comb)[1],
          unlist(biosample_comb)[2],
          protein,
          motif,
          "best_comparison_all_data_points.rds",
          sep = "_"
        )
      ))
      # Initialize a list to store the data frames
      df_list <- list()
      
      # Add the complete data frame
      df_list[["all_states"]] <- merged_df[complete.cases(merged_df), ]
      
      # Add the same-state data frame
      df_list[["same_states"]] <- df_list[["all_states"]] %>%
        filter(!!sym(chromatin_state_col1) == !!sym(chromatin_state_col2))
      
      
      # Add data frames for each unique chromatin state in chromatin_state_col1
      chromatin_states <- names(chromatin_state_colors_short)
      
      for (state in chromatin_states) {
        df_list[[paste0("state_", state)]] <- df_list[["all_states"]] %>%
          filter(!!sym(chromatin_state_col1) == state &
                   !!sym(chromatin_state_col2) == state)
      }
      
      saveRDS(object = df_list,
              file = file.path(sim_output_dir, "df_list.rds"))
      
      for (data_set in names(df_list)) {
        #debug data_set <- names(df_list)[1]
        current_df <- df_list[[data_set]]
        
        # Calculate sample counts for the current data frame
        sample_counts <- table(current_df$sample)
        
        # Skip processing if there are fewer than 3 unique samples
        if (length(sample_counts) < 3) {
          summary_test[[length(summary_test) + 1]] <- data.frame(
            protein = protein,
            biosample1 = unlist(biosample_comb)[1],
            biosample2 = unlist(biosample_comb)[2],
            experiment_ids = df_complete$experiment[which.max(nchar(df_complete$experiment))],
            motif = motif,
            data_set = data_set,
            n_1 = as.numeric(ifelse(
              !is.na(sample_counts[biosample_comb[1]]), sample_counts[biosample_comb[1]], 0
            )),
            n_2 = as.numeric(ifelse(
              !is.na(sample_counts[biosample_comb[2]]), sample_counts[biosample_comb[2]], 0
            )),
            n_both = as.numeric(ifelse(
              !is.na(sample_counts[biosample_both]), sample_counts[biosample_both], 0
            )),
            p_value_S_bigger = NA,
            p_value_S_smaller = NA,
            observed_S_statistic = NA
          )
          
          # Save the plots in protein_plot_list
          protein_plot_list[[protein]][[motif]][[biosample_both]][[data_set]] <- list(
            scatter_plot = ggplot() +
              theme_void() +  # Removes all plot elements
              labs(
                paste(
                  motif,
                  "Scatterplot: CpG Methylation in",
                  biosample_comb[1],
                  "vs",
                  biosample_comb[2]
                )
              ),
            # Optional: Add a title or label,
            histogram_plot = ggplot() +
              theme_void() +  # Removes all plot elements
              labs(title = paste(
                "Histogram of Permuted S Values for", data_set
              ))  # Optional: Add a title or label
          )
          
          print(paste("skippede", data_set, ":", "no 3 lable categories"))
          next
        }
        
        
        # Calculate the observed S statistic
        observed_S_statistic <- calculate_S_statistic(current_df, biosample_comb)
        print(paste("Observed S for", data_set, ":", observed_S_statistic))
        
        # Perform the permutation test
        
        permuted_S_values <- foreach(
          i = 1:n_permutations,
          .combine = 'c',
          .packages = c('dplyr', 'rlang')
        ) %dopar% {
          df_permuted <- current_df %>%
            mutate(sample = sample(sample))  # Randomly permute the sample column
          
          calculate_S_statistic(df_permuted, biosample_comb)
        }
        
        # Calculate p-values
        p_value_S_bigger <- mean(permuted_S_values <= observed_S_statistic)
        p_value_S_smaller <- mean(permuted_S_values >= observed_S_statistic)
        
        # # Save results to the results list
        # results_list[[data_set]] <- list(
        #   observed_S_statistic = observed_S_statistic,
        #   p_value_S_bigger = p_value_S_bigger,
        #   p_value_S_smaller = p_value_S_smaller
        # )
        
        # Convert permuted S values to a data frame for visualization
        permuted_S_df <- data.frame(S_statistic = permuted_S_values)
        
        # Save the permutation results as a CSV
        write.csv(
          permuted_S_df,
          file = file.path(
            sim_output_dir,
            paste0("permuted_S_values_", data_set, ".csv")
          ),
          row.names = FALSE
        )
        
        # Create a histogram of permuted S statistics
        histogram_plot <- ggplot(permuted_S_df, aes(x = S_statistic)) +
          geom_histogram(
            bins = 100,
            fill = "lightblue",
            color = "black"
          ) +
          geom_vline(
            aes(xintercept = observed_S_statistic),
            color = "red",
            linetype = "dashed",
            linewidth = 1
          ) +
          labs(
            title = paste("Histogram of Permuted S Values for", data_set),
            x = "S_statistic",
            y = "Frequency"
          ) +
          theme_minimal()
        
        # Save the histogram plot
        ggsave(
          file = file.path(
            sim_output_dir,
            paste0("S_statistic_histogram_", data_set, ".svg")
          ),
          plot = histogram_plot,
          device = "svg",
          limitsize = FALSE
        )  # Ensure larger plots can be saved
        
        
        # Define the two samples dynamically
        fRead_sample1 <- paste0("fRead_", biosample_comb[1])
        fRead_sample2 <- paste0("fRead_", biosample_comb[2])
        
        # Assign colors dynamically
        current_colors <- group.colors
        current_colors["A549_K562"] <- "gray"
        
        # Modify the `sample` column to include "present in both samples" where appropriate
        current_df <- current_df %>%
          mutate(
            sample_label = case_when(
              sample == biosample_comb[1] ~ biosample_comb[1],
              sample == biosample_comb[2] ~ biosample_comb[2],
              TRUE ~ "present in both samples"
            )
          )
        
        scatter_plot <- ggplot(current_df) +
          geom_point(aes(
            x = !!sym(paste0("fRead_", biosample_comb[1])),
            y = !!sym(paste0("fRead_", biosample_comb[2])),
            color = sample_label
          ), size = 1.5) +
          scale_color_manual(
            values = c(group.colors, "present in both samples" = "gray"),
            name = "Biosample"
          ) +
          labs(
            title = paste(
              motif,
              "Scatterplot: CpG Methylation in",
              biosample_comb[1],
              "vs",
              biosample_comb[2]
            ),
            x = paste("CpG Methylation in", biosample_comb[1]),
            y = paste("CpG Methylation in", biosample_comb[2])
          ) +
          theme_minimal()
        
        # Ensure the directory for saving the plot exists
        scatter_output_dir <- file.path(sim_output_dir)
        if (!dir.exists(scatter_output_dir)) {
          dir.create(scatter_output_dir, recursive = TRUE)
        }
        
        # Save the scatter plot as an SVG file
        scatter_plot_file <- file.path(
          scatter_output_dir,
          paste0(
            motif,
            "_Scatterplot_CpG_Methylation_",
            biosample_comb[1],
            "_vs_",
            biosample_comb[2],
            ".svg"
          )
        )
        
        ggsave(
          filename = scatter_plot_file,
          plot = scatter_plot,
          width = 8,
          # Width of the plot in inches
          height = 6,
          # Height of the plot in inches
          device = "svg"  # Specify the SVG format
        )
        
        
        # Save the plots in protein_plot_list
        protein_plot_list[[protein]][[motif]][[biosample_both]][[data_set]] <- list(scatter_plot = scatter_plot, histogram_plot = histogram_plot)
        
        summary_test[[length(summary_test) + 1]] <- data.frame(
          protein = protein,
          biosample1 = unlist(biosample_comb)[1],
          biosample2 = unlist(biosample_comb)[2],
          experiment_ids = df_complete$experiment[which.max(nchar(df_complete$experiment))],
          motif = motif,
          data_set = data_set,
          n_1 = sample_counts[biosample_comb[1]],
          n_2 = sample_counts[biosample_comb[2]],
          n_both = sample_counts[biosample_both],
          p_value_S_bigger = p_value_S_bigger,
          p_value_S_smaller = p_value_S_smaller,
          observed_S_statistic = observed_S_statistic
        )
        
      }
      
      
    } # finish biosample loop
    
    
  }# finish motif loop
  
}
# Stop the cluster after all loops are done
stopCluster(cl)
# Combine the list into a single data frame
summary_test <- do.call(rbind, summary_test)
# Remove row names explicitly
rownames(summary_test) <- NULL
saveRDS(object = summary_test,
        file = file.path(output_folder, "summary_test_100k.rds"))
# summary_test[ !is.na(summary_test$p_value_S_smaller_state1) & (1/10000 > summary_test$p_value_S_smaller_state1 | 1/10000 > summary_test$p_value_S_bigger_state1),]



################## new versiomstrtisfied stat  ###########################################

# Function to calculate S_statistic for a given data frame
calculate_S_statistic <- function(df, biosample_comb) {
  mean_diff_per_sample <- df %>%
    mutate(diff_fRead = !!sym(paste0("fRead_", biosample_comb[1])) -!!sym(paste0("fRead_", biosample_comb[2]))) %>%  # Calculate difference
    group_by(sample) %>%  # Group by sample
    summarize(mean_diff = mean(diff_fRead, na.rm = TRUE))  # Calculate mean difference per sample
  
  # Calculate the S statistic
  S_statistic <- mean_diff_per_sample$mean_diff[biosample_comb[1] == mean_diff_per_sample$sample] -
    mean_diff_per_sample$mean_diff[biosample_comb[2] == mean_diff_per_sample$sample]
  
  return(S_statistic)
}

library(doParallel)
library(foreach)

# Set up the parallel backend (before the outer loop)
num_cores <- detectCores() - 1  # Use one fewer core than available
cl <- makeCluster(num_cores)
registerDoParallel(cl)

stratified_test <- list()
n_permutations <- 100000

start <- 1 #which("NR3C1" == protein_folders) -1
stop <- #start +1
  stop <- length(protein_folders)

# Loop over protein folders
for (protein in protein_folders[start:stop]) {
  # [start:stop]
  # debug protein <- "CTCF"
  # debug protein <- protein_folders[1]
  print(protein)
  
  protein_plot_list <- list()
  protein_data <- filtered_summary_df[filtered_summary_df$protein == protein, ]
  
  # Loop over unique motifs
  unique_motifs <- unique(protein_data$motif)
  protein_output_dir <- file.path(output_folder, protein)
  
  for (motif in unique_motifs) {
    # debug motif<- unique_motifs[1]
    
    print(motif)
    motif_data <- protein_data[protein_data$motif == motif, ]
    
    # Generate all combinations of two biosamples
    biosample_combinations <- combn(unique(motif_data$biosample), 2, simplify = FALSE)
    motif_out_dir <- file.path(protein_output_dir, motif)
    
    # Loop through each combination of two biosamples
    for (biosample_comb in biosample_combinations) {
      # debug biosample_comb <- biosample_combinations[[1]]
      print(biosample_comb)
      biosample_both <- paste(unlist(biosample_comb), collapse = "_")
      sim_output_dir <- file.path(motif_out_dir, biosample_both)
      dir.create(sim_output_dir,
                 recursive = TRUE,
                 showWarnings = FALSE)
      
      chromatin_state_col1 <- paste0("Chromatin_State_", biosample_comb[1])
      chromatin_state_col2 <- paste0("Chromatin_State_", biosample_comb[2])
      
      # Load merged data
      rds_file <- file.path(
        protein_output_dir,
        paste(
          biosample_both,
          protein,
          motif,
          "best_comparison_all_data_points.rds",
          sep = "_"
        )
      )
      
      merged_df <- readRDS(rds_file)
      merged_df <- merged_df[complete.cases(merged_df), ]
      
      # Filter CpGs with same chromatin state
      current_df <- merged_df %>%
        filter(!!sym(chromatin_state_col1) == !!sym(chromatin_state_col2))
      
      sample_counts <- table(current_df$sample)
      
      if (length(sample_counts) < 3) {
        stratified_test[[length(stratified_test) + 1]] <- data.frame(
          protein = protein,
          biosample1 = biosample_comb[1],
          biosample2 = biosample_comb[2],
          motif = motif,
          data_set = "stratified_test",
          # Placeholder; update as needed
          n_1 = ifelse(
            biosample_comb[1] %in% names(sample_counts),
            as.numeric(sample_counts[biosample_comb[1]]),
            0
          ),
          n_2 = ifelse(
            biosample_comb[2] %in% names(sample_counts),
            as.numeric(sample_counts[biosample_comb[2]]),
            0
          ),
          n_both = ifelse(
            biosample_both %in% names(sample_counts),
            as.numeric(sample_counts[biosample_both]),
            0
          ),
          p_value_S_bigger = NA,
          p_value_S_smaller = NA,
          observed_S_statistic = NA
        )
        
        next
      }
      
      
      # Calculate observed S statistic
      observed_S_statistic <- calculate_S_statistic(current_df, biosample_comb)
      #print(paste("Observed S for", protein, motif, ":", observed_S_statistic))
      
      # Define the file path for the saved permutation results
      permuted_file <- file.path(sim_output_dir,
                                 paste0("permuted_stratified_S_df", ".csv"))
      
      # Check if the file already exists
      if (file.exists(permuted_file)) {
        # Load the existing permutation results
        permuted_stratified_S_df <- read.csv(permuted_file)
        permuted_S_values <- permuted_stratified_S_df$S_statistic
      } else {
        # Perform permutation test
        permuted_S_values <- foreach(
          i = 1:n_permutations,
          .combine = 'c',
          .packages = c('dplyr', 'rlang')
        ) %dopar% {
          current_df %>%
            group_by(!!sym(chromatin_state_col1),!!sym(chromatin_state_col2)) %>%
            mutate(sample = sample(sample)) %>%
            ungroup() %>%
            calculate_S_statistic(biosample_comb)
        }
      }
      p_value_S_bigger <- mean(permuted_S_values >= observed_S_statistic)
      p_value_S_smaller <- mean(permuted_S_values <= observed_S_statistic)
      
      cat(
        sprintf(
          "Observed S statistic: %.4f | Out of %d permutations, %d >= observed (p=%s), %d <= observed (p=%s), two-sided p=%s\n",
          observed_S_statistic,
          length(permuted_S_values),
          sum(permuted_S_values >= observed_S_statistic),
          ifelse(
            p_value_S_bigger == 0,
            "< 1/length(permuted_S_values)",
            sprintf("%.4f", p_value_S_bigger)
          ),
          sum(permuted_S_values <= observed_S_statistic),
          ifelse(
            p_value_S_smaller == 0,
            "< 1/length(permuted_S_values)",
            sprintf("%.4f", p_value_S_smaller)
          ),
          ifelse(
            2 * min(p_value_S_bigger, p_value_S_smaller) == 0,
            "< 1/length(permuted_S_values)",
            sprintf("%.4f", 2 * min(
              p_value_S_bigger, p_value_S_smaller
            ))
          )
        )
      )
      
      plot_flag <- FALSE
      if (plot_flag) {
        # Convert permuted S values to a data frame for visualization
        permuted_stratified_S_df <- data.frame(S_statistic = permuted_S_values)
        
        # Save the permutation results as a CSV
        write.csv(
          permuted_stratified_S_df,
          file = file.path(
            sim_output_dir,
            paste0("permuted_stratified_S_df", ".csv")
          ),
          row.names = FALSE
        )
        
        # Create a histogram of permuted S statistics
        histogram_plot <- ggplot(permuted_stratified_S_df, aes(x = S_statistic)) +
          geom_histogram(
            bins = 100,
            fill = "lightblue",
            color = "black"
          ) +
          geom_vline(
            aes(xintercept = observed_S_statistic),
            color = "red",
            linetype = "dashed",
            linewidth = 1
          ) +
          labs(
            title = paste("Histogram of Permuted S Values for", data_set),
            x = "S_statistic",
            y = "Frequency"
          ) +
          theme_minimal()
        
        # Save the histogram plot
        ggsave(
          file = file.path(
            sim_output_dir,
            paste0("S_stratified_statistic_histogram_", ".svg")
          ),
          plot = histogram_plot,
          device = "svg",
          limitsize = FALSE
        )  # Ensure larger plots can be saved
        
        
        # Define the two samples dynamically
        fRead_sample1 <- paste0("fRead_", biosample_comb[1])
        fRead_sample2 <- paste0("fRead_", biosample_comb[2])
        
        # Assign colors dynamically
        current_colors <- group.colors
        current_colors[biosample_both] <- "gray"
        
        # Modify the `sample` column to include "present in both samples" where appropriate
        current_df <- current_df %>%
          mutate(
            sample_label = case_when(
              sample == biosample_comb[1] ~ biosample_comb[1],
              sample == biosample_comb[2] ~ biosample_comb[2],
              TRUE ~ "present in both samples"
            )
          )
        
        scatter_plot <- ggplot(current_df) +
          geom_point(aes(
            x = !!sym(paste0("fRead_", biosample_comb[1])),
            y = !!sym(paste0("fRead_", biosample_comb[2])),
            color = sample_label
          ), size = 1.5) +
          scale_color_manual(
            values = c(group.colors, "present in both samples" = "gray"),
            name = "Biosample"
          ) +
          labs(
            title = paste(
              motif,
              "Scatterplot: CpG Methylation in",
              biosample_comb[1],
              "vs",
              biosample_comb[2]
            ),
            x = paste("CpG Methylation in", biosample_comb[1]),
            y = paste("CpG Methylation in", biosample_comb[2])
          ) +
          theme_minimal()
        
        # Ensure the directory for saving the plot exists
        scatter_output_dir <- file.path(sim_output_dir)
        if (!dir.exists(scatter_output_dir)) {
          dir.create(scatter_output_dir, recursive = TRUE)
        }
        
        # Save the scatter plot as an SVG file
        scatter_plot_file <- file.path(
          scatter_output_dir,
          paste0(
            motif,
            "_Scatterplot_CpG_Methylation_",
            biosample_comb[1],
            "_vs_",
            biosample_comb[2],
            "_stratified.svg"
          )
        )
        
        ggsave(
          filename = scatter_plot_file,
          plot = scatter_plot,
          width = 8,
          # Width of the plot in inches
          height = 6,
          # Height of the plot in inches
          device = "svg"  # Specify the SVG format
        )
      }
      
      stratified_test[[length(stratified_test) + 1]] <- data.frame(
        protein = protein,
        biosample1 = biosample_comb[1],
        biosample2 = biosample_comb[2],
        motif = motif,
        data_set = "stratified_test",
        # Placeholder; update as needed
        n_1 = sample_counts[biosample_comb[1]],
        n_2 = sample_counts[biosample_comb[2]],
        n_both = sample_counts[biosample_both],
        p_value_S_bigger = p_value_S_bigger,
        p_value_S_smaller = p_value_S_smaller,
        observed_S_statistic = observed_S_statistic
      )
      
    }
  }
}  
  # Stop the cluster after all loops are done
  stopCluster(cl)
  
  # Combine results into a single data frame and save
  stratified_test <- do.call(rbind, stratified_test)
  rownames(stratified_test) <- NULL
  saveRDS(stratified_test,
          file = file.path(output_folder, "stratified_test_100k.rds"))
  
  
######################## END STRAT TEST ###########################################
  
  
  
  
  
  
  
  # Save the dataframe as an Excel sheet using openxlsx functions
  openxlsx::write.xlsx(
    summary_test,
    file = file.path(output_folder, "summary_test.xlsx"),
    overwrite = TRUE
  )
  
  
  summary_test <- readRDS(file = file.path(output_folder, "summary_test.rds"))
  
  # Loop over protein folders
  for (protein in protein_folders) {
    # [1:1]
    # protein <- "CTCF"
    print(protein)
    
    protein_data <- filtered_summary_df[filtered_summary_df$protein == protein, ]
    # Filter the data for the specific protein
    
    # Loop over unique motifs
    unique_motifs <- unique(protein_data$motif)
    
    for (motif in unique_motifs) {
      # Filter the data for the current motif
      # debug # motif <- unique_motifs[1]
      print(motif)
      motif_data <- protein_data[protein_data$motif == motif, ]
      
      # Generate all combinations of two biosamples
      biosample_combinations <- combn(unique(motif_data$biosample), 2, simplify = FALSE)
      
      # Loop through each combination of two biosamples
      for (biosample_comb in biosample_combinations) {
        print(biosample_comb)
        
        # Dynamically generate the chromatin state column names for both biosamples
        chromatin_state_col1 <- paste0("Chromatin_State_", biosample_comb[1])
        chromatin_state_col2 <- paste0("Chromatin_State_", biosample_comb[2])
        
        
        merged_df <- readRDS(file = file.path(
          output_folder,
          protein,
          paste(
            unlist(biosample_comb)[1],
            unlist(biosample_comb)[2],
            protein,
            motif,
            "best_comparison_all_data_points.rds",
            sep = "_"
          )
        ))
        
        biosample_both <- paste(unlist(biosample_comb)[1],
                                unlist(biosample_comb)[2],
                                sep = "_")
        
        df_complete <- merged_df[complete.cases(merged_df), ]
        
        # Filter rows where chromatin states are the same for both samples
        df_same_state <- df_complete %>%
          filter(!!sym(chromatin_state_col1) == !!sym(chromatin_state_col2))
        
        # Filter rows where chromatin states are "1" in both samples
        df_state1 <- df_complete %>%
          filter(!!sym(chromatin_state_col1) == "1" &
                   !!sym(chromatin_state_col2) == "1")
        
        
        sample_counts <- table(df_complete$sample)
        sample_counts_same_state <- table(df_same_state$sample)
        sample_counts_state1 <- table(df_state1$sample)
        
        # Dynamically generate the column names for fRead values
        fRead_col1 <- paste0("fRead_", biosample_comb[1])
        fRead_col2 <- paste0("fRead_", biosample_comb[2])
        
        # Create the scatter plot
        ggplot(df_state1,
               aes_string(
                 x = fRead_col1,
                 y = fRead_col2,
                 color = "sample"
               )) +
          geom_point() +
          labs(
            title = paste(
              protein,
              motif,
              "Scatter Plot of fRead_",
              biosample_comb[1],
              " vs fRead_",
              biosample_comb[2]
            ),
            x = paste("fRead_", biosample_comb[1]),
            y = paste("fRead_", biosample_comb[2])
          ) +
          theme_minimal()
        
      }
    }
  }