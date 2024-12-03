#################################################################
##  Chip-seq vs methylation vs Chromatin
##
##  input: everything in HG38
##
##          Encode3\meme\fimo_single_experiments_on_known_motif_CpG_and_methylation_and_peaks
##
##
##  output:   Encode3\simulation
##
##  v_1 24.09.2024
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
# Load necessary libraries
library(stringr)
library(readr)
library(foreach)
library(doParallel)
library(dplyr)
library(rlang)
library(ggplot2)
library(tidyr)
library(VennDiagram)
##################################### INPUT ########################################



################################## constants #####################################
start_script <- Sys.time()
# Set input and output directories
input_folder <- file.path(
  this.dir,
  "meme",
  "fimo_single_experiments_on_known_motif_and_peaks_add_methylation_v3"
)
output_folder <- file.path(this.dir, "simulation")
RNA_seq_folder <- "RNA-seq"

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

# Define chromatin state colors
chromatin_state_colors <- c(
  "1" = "#FF0000",  "2" = "#FF4500",  "3" = "#FF9900",  "4" = "#FFCC00",
  "5" = "#00CC00",  "6" = "#006400",  "7" = "#FFD700",  "8" = "#FFD700",
  "9" = "#FFFF00",  "10" = "#FFDD00", "11" = "#FFEA73", "12" = "#9370DB",
  "13" = "#C0C0C0", "14" = "#FF4500", "15" = "#FFDD00", "16" = "#808080",
  "17" = "#A9A9A9", "18" = "#000000"
)


# set up parallel processing
n_cores <- detectCores() - 1

CHR_NAMES <- paste0("chr", c(1:22))

chromatin_state_names <-
  c(
    "1_TssA",
    "2_TssFlnk",
    "3_TssFlnkU",
    "4_TssFlnkD",
    "5_Tx",
    "6_TxWk",
    "7_EnhG1",
    "8_EnhG2",
    "9_EnhA1",
    "10_EnhA2",
    "11_EnhWk",
    "12_ZNF_Rpts",
    "13_Het",
    "14_TssBiv",
    "15_EnhBiv",
    "16_ReprPC",
    "17_ReprPCWk",
    "18_Quies"
  )

# readLines(bed_file, n = 1)
#file_path <- bed_file
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

# Get the list of all files matching the pattern *_methylation.chromatinstate.fimo.bed
protein_folders <- list.dirs(path = input_folder,
                             recursive = FALSE ,
                             full.names = FALSE)
#List all protein folders, excluding "sbatch_scripts" and "log"
protein_folders <- protein_folders[!grepl("sbatch_scripts|log", protein_folders)]

if (FALSE) {
  # debug
  #protein_folders <- protein_folders[1:4]
  
  # Initialize an empty list to store results for each file
  summary_list <- list()
  
  ##### Big data creation loop
  
  for (protein in protein_folders) {
    # [1:1]
    print(protein)  # protein <- protein_folders[1]
    
    # Get all files for the current protein
    files <- list.files(path = file.path(input_folder, protein))
    
    df_protein <- data.frame()
    
    for (bed_file in files) {
      print(bed_file)  # bed_file <- files[1]
      
      # Split the file name by "_"
      parts <- strsplit(bed_file, "_")[[1]]
      
      # Extract the components based on position
      biosample <- parts[1]            # Example: GM12878
      experiment_id <- parts[3]        # Example: ENCFF320KXO
      motif <- sub("\\.bed$", "", parts[4])  # Extract motif
      
      
      # Read the BED file
      df <- read_bed_file(file_path = file.path(input_folder, protein, bed_file))
      
      # Create a data frame with only unique peaks based on (chr_peak, start_peak, end_peak)
      df_unique_peaks <- df %>%
        distinct(chr, start_peak, end_peak, .keep_all = TRUE)
      
      # Calculate the different categories based on the unique peaks
      n_peak_with_CpG <- sum(df_unique_peaks$start_motif != -1 &
                               df_unique_peaks$start_cg != 0)
      n_peak_no_CpG <- sum(df_unique_peaks$start_motif != -1 &
                             df_unique_peaks$start_cg == 0)
      n_peak_no_motif_no_CpG <- sum(df_unique_peaks$start_motif == -1 &
                                      df_unique_peaks$start_cg == 0)
      
      # Count total unique peaks
      n_peaks <- nrow(df_unique_peaks)
      
      
      df_unique_CpG <- df_unique_peaks[df_unique_peaks$start_motif != -1 & df_unique_peaks$start_cg != 0,]
      
      state_column <- paste0("Chromatin_State_", biosample)
      
      # Define the states you want to count
      states <- c(".", as.character(1:18))
      
      # Count occurrences of each state dynamically
      state_counts <- setNames(rep(0, length(states)), states) # Initialize with 0
      state_counts_dynamic <- table(df_unique_CpG[[state_column]])
      state_counts[names(state_counts_dynamic)] <- as.numeric(state_counts_dynamic)
      
      # Convert to a named vector with all states accounted for
      state_counts <- state_counts
      
      
      
      # Collect the results in a list
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
      
      
      df_protein <- rbind(df_protein,df)
      
    }
    
    
    
    # Loop over unique motifs
    unique_motifs <- unique(df_protein$motif)
    
    for (motif in unique_motifs) {
      # Filter the data for the current motif
      # debug motif <- unique_motifs[1]
      print(motif)
      df_motif <- df_protein %>% filter(motif == !!motif)
      
      # Filter data for the specific sample (GM12878)
      data_filtered <- df_motif %>%
        filter(sample == "GM12878")
      
      # Calculate proportions of ChromHMM states for each protein, experiment, and motif
      state_proportions <- data_filtered %>%
        group_by(protein, experiment, motif, Chromatin_State_GM12878) %>%
        summarise(Count = n(), .groups = "drop") %>%
        group_by(protein, experiment, motif) %>%
        mutate(
          Total = sum(Count),
          Proportion = Count / Total
        )
      
      state_proportions <- data_filtered %>%
        group_by(protein, experiment, motif, Chromatin_State_GM12878) %>%
        summarise(Count = n(), .groups = "drop") %>%
        group_by(protein, experiment, motif) %>%
        mutate(
          Total = sum(Count),
          Proportion = Count / Total
        ) %>%
        mutate(
          Chromatin_State_GM12878 = factor(Chromatin_State_GM12878, levels = names(chromatin_state_colors))
        )
      
      # Plot the proportions with custom colors
      ggplot(state_proportions, aes(x = Chromatin_State_GM12878, y = Proportion, fill = Chromatin_State_GM12878)) +
        geom_bar(stat = "identity") +
        facet_wrap(~ protein + experiment + motif, scales = "free_y") +
        scale_fill_manual(values = chromatin_state_colors) +
        scale_y_continuous(labels = scales::percent_format()) +
        labs(
          title = "Proportions of Chromatin States (GM12878)",
          x = "Chromatin State",
          y = "Proportion of Entries",
          fill = "Chromatin State"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(size = 10)
        )
      
    }
    
    
  }
  
  # Combine the list into a single data frame
  summary_df <- do.call(rbind, summary_list)
  
  # Calculate the ratio of n_peak_no_motif_no_CpG to n_peaks
  summary_df$no_motif_no_CpG_ratio <- summary_df$n_peak_no_motif_no_CpG / summary_df$n_peaks
  
  # Display the summary DataFrame
  print(dim(summary_df))
  
  saveRDS(object = summary_df,
          file = file.path(output_folder, "summary_df.rds"))
  rm(summary_list)
  gc()
  
  
  
  # Plot the histogram
  histogram_plot <- ggplot(summary_df, aes(x = no_motif_no_CpG_ratio)) +
    geom_histogram(
      binwidth = 0.05,
      fill = "#00BFC4",
      color = "black",
      alpha = 0.7
    ) +
    labs(title = "Distribution of No Motif No CpG Peaks Ratio", x = "Ratio of No Motif No CpG Peaks", y = "Count") +
    theme_minimal()
  
  # Save the absolute plot
  ggsave(
    filename = file.path(output_folder, "histogram_plot_fractions.png"),
    plot = histogram_plot,
    width = 10,
    height = 6
  )
  
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
}
summary_df <- readRDS(file = file.path(output_folder, "summary_df.rds"))

# Filter summary_df to keep only the rows with the lowest no_motif_no_CpG_ratio for each (protein, biosample, motif) combination
filtered_summary_df <- summary_df %>%
  group_by(protein, biosample, motif) %>%
  slice_min(no_motif_no_CpG_ratio, with_ties = FALSE) %>%
  ungroup()


# Loop over protein folders
for (protein in protein_folders) {
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
        "fRead_A549",
        "mRead_GM12878",
        "nRead_GM12878",
        "fRead_GM12878",
        "mRead_HepG2",
        "nRead_HepG2",
        "fRead_HepG2",
        "mRead_K562",
        "nRead_K562",
        "fRead_K562"
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
        )
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
      
      saveRDS(object = merged_df,
              file =  file.path(
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
      
    }
  }
  
}


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

start <- which("NR3C1" == protein_folders)-1
stop <- which("NR3C1" == protein_folders)+1
# [start:stop]
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
      
      biosample_both <- paste(unlist(biosample_comb)[1], unlist(biosample_comb)[2], sep = "_")
      
      # Create output directory for each protein if it doesn't exist
      sim_output_dir <- file.path(output_folder, protein, motif, biosample_both)
      if (!dir.exists(sim_output_dir)) {
        dir.create(sim_output_dir, recursive = TRUE)
      }
      
      df_complete <- merged_df[complete.cases(merged_df), ]
      
      # Filter rows where chromatin states are the same for both samples
      df_same_state <- df_complete %>%
        filter(!!sym(chromatin_state_col1) == !!sym(chromatin_state_col2))
      
      # Filter rows where chromatin states are "1" in both samples
      df_state1 <- df_complete %>%
        filter(!!sym(chromatin_state_col1) == "1" & !!sym(chromatin_state_col2) == "1")
      
      
      sample_counts <- table(df_complete$sample)
      sample_counts_same_state <- table(df_same_state$sample)
      sample_counts_state1 <- table(df_state1$sample)
      set.seed(123)  # For reproducibility
      n_permutations <- 10000
      
      
      if (dim(sample_counts) < 3) {
        observed_S_statistic <- NA
        p_value_S_bigger <- NA
        p_value_S_smaller <- NA
        
      } else{
        # Calculate the observed S_statistic
        observed_S_statistic <- calculate_S_statistic(df_complete, biosample_comb)
        print(observed_S_statistic)
        
        # Perform the permutation test
        
        permuted_S_values <- numeric(n_permutations)
        # Parallelized permutation test (using the cluster initialized once)
        permuted_S_values <- foreach(
          i = 1:n_permutations,
          .combine = 'c',
          .packages = c('dplyr', 'rlang')
        ) %dopar% {
          df_permuted <- df_complete %>%
            mutate(sample = sample(sample))  # Randomly permute the sample column
          
          calculate_S_statistic(df_permuted, biosample_comb)
        }
        # Calculate the p-value (proportion of permuted S values greater or equal to observed S)
        p_value_S_bigger <- mean(permuted_S_values <= observed_S_statistic)
        p_value_S_smaller <- mean(permuted_S_values >= observed_S_statistic)
        # Convert permuted_S_values into a data frame for ggplot
        permuted_S_df <- data.frame(S_statistic = permuted_S_values)
        write.csv(
          permuted_S_df,
          file = file.path(
            sim_output_dir,
            "permuted_S_values_all_states_no_filter.csv"
          ),
          row.names = FALSE
        )
        # Create the histogram with a red vertical line for observed_S_statistic
        histogram_plot <- ggplot(permuted_S_df, aes(x = S_statistic)) +
          geom_histogram(
            binwidth = 0.01,
            fill = "lightblue",
            color = "black"
          ) +  # Histogram of permuted values
          geom_vline(
            aes(xintercept = observed_S_statistic),
            color = "red",
            linetype = "dashed",
            linewidth = 1
          ) +  # Red vertical line
          labs(
            title = paste(
              "Histogram of Permuted S_statistic Values all states",
              motif
            ),
            x = "S_statistic",
            y = "Frequency"
          ) +
          theme_minimal()  # Optional: Apply a clean theme
        
        # Save the plot as a .png file
        ggsave(file.path(
          sim_output_dir,
          "S_statistic_histogram_no_filter.png"
        ),
        plot = histogram_plot)
        
      }
      
      
      
      if (dim(sample_counts_same_state) < 3) {
        observed_S_statistic_same_state <- NA
        p_value_S_bigger_same_state <- NA
        p_value_S_smaller_same_state <- NA
      } else{
        observed_S_statistic_same_state <- calculate_S_statistic(df_same_state, biosample_comb)
        
        permuted_S_values_same_state <- numeric(n_permutations)
        permuted_S_values_same_state <- foreach(
          i = 1:n_permutations,
          .combine = 'c',
          .packages = c('dplyr', 'rlang')
        ) %dopar% {
          df_permuted <- df_same_state %>%
            mutate(sample = sample(sample))  # Randomly permute the sample column
          
          calculate_S_statistic(df_permuted, biosample_comb)
        }
        # Calculate the p-value (proportion of permuted S values greater or equal to observed S)
        p_value_S_bigger_same_state <- mean(permuted_S_values_same_state <= observed_S_statistic_same_state)
        p_value_S_smaller_same_state <- mean(permuted_S_values_same_state >= observed_S_statistic_same_state)
        
        # Convert permuted_S_values into a data frame for ggplot
        permuted_S_df_same_state <- data.frame(S_statistic = permuted_S_values_same_state)
        
        write.csv(
          permuted_S_df,
          file = file.path(
            sim_output_dir,
            "permuted_S_values_all_states_same_state_filter.csv"
          ),
          row.names = FALSE
        )
        
        # Create the histogram with a red vertical line for observed_S_statistic
        histogram_plot_same_state <- ggplot(permuted_S_df_same_state, aes(x = S_statistic)) +
          geom_histogram(
            binwidth = 0.01,
            fill = "lightblue",
            color = "black"
          ) +  # Histogram of permuted values
          geom_vline(
            aes(xintercept = observed_S_statistic_same_state),
            color = "red",
            linetype = "dashed",
            linewidth = 1
          ) +  # Red vertical line
          labs(
            title = paste(
              "Histogram of Permuted S_statistic Values same states",
              motif
            ),
            x = "S_statistic",
            y = "Frequency"
          ) +
          theme_minimal()  # Optional: Apply a clean theme
        
        # Save the plot as a .png file
        ggsave(
          file.path(
            sim_output_dir,
            "S_statistic_histogram_same_state_filter.png"
          ),
          plot = histogram_plot_same_state
        )
      }
      
      # state 1 only 
      #       sample_counts_state1 <- table(df_state1$sample)
      if (dim(sample_counts_state1) < 3) {
        observed_S_statistic_state1 <- NA
        p_value_S_bigger_state1 <- NA
        p_value_S_smaller_state1 <- NA
      } else{
        observed_S_statistic_state1 <- calculate_S_statistic(df_state1, biosample_comb)
        
        permuted_S_values_state1 <- numeric(n_permutations)
        permuted_S_values_state1 <- foreach(
          i = 1:n_permutations,
          .combine = 'c',
          .packages = c('dplyr', 'rlang')
        ) %dopar% {
          df_permuted <- df_state1 %>%
            mutate(sample = sample(sample))  # Randomly permute the sample column
          
          calculate_S_statistic(df_permuted, biosample_comb)
        }
        # Calculate the p-value (proportion of permuted S values greater or equal to observed S)
        p_value_S_bigger_state1 <- mean(permuted_S_values_state1 <= observed_S_statistic_state1)
        p_value_S_smaller_state1 <- mean(permuted_S_values_state1 >= observed_S_statistic_state1)
        
        # Convert permuted_S_values into a data frame for ggplot
        permuted_S_df_state1 <- data.frame(S_statistic = permuted_S_values_state1)
        
        write.csv(
          permuted_S_df,
          file = file.path(
            sim_output_dir,
            "permuted_S_values_all_states_state1_filter.csv"
          ),
          row.names = FALSE
        )
        
        # Create the histogram with a red vertical line for observed_S_statistic
        histogram_plot_state1 <- ggplot(permuted_S_df_state1, aes(x = S_statistic)) +
          geom_histogram(
            binwidth = 0.01,
            fill = "lightblue",
            color = "black"
          ) +  # Histogram of permuted values
          geom_vline(
            aes(xintercept = observed_S_statistic_state1),
            color = "red",
            linetype = "dashed",
            linewidth = 1
          ) +  # Red vertical line
          labs(
            title = paste(
              "Histogram of Permuted S_statistic Values same states",
              motif
            ),
            x = "S_statistic",
            y = "Frequency"
          ) +
          theme_minimal()  # Optional: Apply a clean theme
        
        # Save the plot as a .png file
        ggsave(
          file.path(
            sim_output_dir,
            "S_statistic_histogram_state1_filter.png"
          ),
          plot = histogram_plot_state1
        )
      }
      
      # Collect the results in a list
      summary_test[[length(summary_test) + 1]] <- data.frame(
        protein = protein,
        biosample1 = unlist(biosample_comb)[1],
        biosample2 = unlist(biosample_comb)[2],
        experiment_ids = df_complete$experiment[which.max(nchar(df_complete$experiment))],
        motif = motif,
        n_1 = sample_counts[biosample_comb[1]],
        n_2 = sample_counts[biosample_comb[2]],
        n_both = sample_counts[biosample_both],
        n_1_same_state = sample_counts_same_state [biosample_comb[1]],
        n_2_same_state  = sample_counts_same_state [biosample_comb[2]],
        n_both_same_state  = sample_counts_same_state [biosample_both],
        n_1_state1 = sample_counts_state1 [biosample_comb[1]],
        n_2_state1  = sample_counts_state1 [biosample_comb[2]],
        n_both_state1  = sample_counts_state1 [biosample_both],
        p_value_S_bigger = p_value_S_bigger,
        p_value_S_smaller = p_value_S_smaller,
        p_value_S_bigger_same_state = p_value_S_bigger_same_state,
        p_value_S_smaller_same_state = p_value_S_smaller_same_state,
        p_value_S_bigger_state1 = p_value_S_bigger_state1,
        p_value_S_smaller_state1 = p_value_S_smaller_state1,
        observed_S_statistic = observed_S_statistic,
        observed_S_statistic_same_state = observed_S_statistic_same_state,
        observed_S_statistic_state1 = observed_S_statistic_state1
      )
      
    }
  }
  
}
# Stop the cluster after all loops are done
stopCluster(cl)
# Combine the list into a single data frame
summary_test <- do.call(rbind, summary_test)
# Remove row names explicitly
rownames(summary_test) <- NULL
saveRDS(object = summary_test,
        file = file.path(output_folder, "summary_test.rds"))
# summary_test[ !is.na(summary_test$p_value_S_smaller_state1) & (1/10000 > summary_test$p_value_S_smaller_state1 | 1/10000 > summary_test$p_value_S_bigger_state1),]

# Save the dataframe as an Excel sheet using openxlsx functions
openxlsx::write.xlsx(summary_test, file = file.path(output_folder, "summary_test.xlsx"), overwrite = TRUE)


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
      
      biosample_both <- paste(unlist(biosample_comb)[1], unlist(biosample_comb)[2], sep = "_")
      
      df_complete <- merged_df[complete.cases(merged_df), ]
      
      # Filter rows where chromatin states are the same for both samples
      df_same_state <- df_complete %>%
        filter(!!sym(chromatin_state_col1) == !!sym(chromatin_state_col2))
      
      # Filter rows where chromatin states are "1" in both samples
      df_state1 <- df_complete %>%
        filter(!!sym(chromatin_state_col1) == "1" & !!sym(chromatin_state_col2) == "1")
      
      
      sample_counts <- table(df_complete$sample)
      sample_counts_same_state <- table(df_same_state$sample)
      sample_counts_state1 <- table(df_state1$sample)
     
      # Dynamically generate the column names for fRead values
      fRead_col1 <- paste0("fRead_", biosample_comb[1])
      fRead_col2 <- paste0("fRead_", biosample_comb[2])
      
      # Create the scatter plot
      ggplot(df_state1, aes_string(x = fRead_col1, y = fRead_col2, color = "sample")) +
        geom_point() +
        labs(
          title = paste(protein,motif,"Scatter Plot of fRead_", biosample_comb[1], " vs fRead_", biosample_comb[2]),
          x = paste("fRead_", biosample_comb[1]),
          y = paste("fRead_", biosample_comb[2])
        ) +
        theme_minimal()
      
    }
  }
}