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
  "fimo_single_experiments_on_known_motif_and_peaks_add_methylation_v2"
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
    
    # Start by subsetting df into df_compact and adding biosample, protein, experiment_id, motif
    df_compact <- df[, c("chr",
                         "start_cg",
                         "end_cg",
                         "Chromatin_State",
                         "fRead",
                         "strand")]
    df_compact$fRead <- as.numeric(df$mRead) / as.numeric(df$nRead)
    df_compact$biosample <- biosample
    df_compact$protein <- protein
    df_compact$experiment_id <- experiment_id
    df_compact$motif <- motif
    df_compact <- df_compact[df_compact$start_cg != 0, ]
    df_protein <- rbind(df_protein, df_compact)
    
    
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
  }
  
  #  Replace "." with NA in all columns (character and numeric-like)
  df_protein <- df_protein %>%
    mutate(across(everything(), ~ ifelse(. == ".", NA, .)))  # Replace "." with NA for all columns
  
  #  Convert numeric-like columns to integer or numeric
  df_protein <- df_protein %>%
    mutate(across(c(start_cg, end_cg, fRead), as.numeric), # Convert numeric-like columns to numeric
           across(
             c(
               chr,
               strand,
               Chromatin_State,
               biosample,
               protein,
               experiment_id,
               motif
             ),
             as.factor
           ))  # Convert character columns to factors)
  
  # Convert numeric columns to integer where applicable, and keep others as numeric
  df_protein <- df_protein %>%
    mutate(across(c(start_cg, end_cg), ~ ifelse(is.na(.), NA_integer_, as.integer(.))))  # Convert to integer if applicable)
  
  
  # Reorder Chromatin_State by numeric value
  df_protein$Chromatin_State <- factor(df_protein$Chromatin_State, levels = as.character(sort(as.numeric(
    levels(df_protein$Chromatin_State)
  ))))
  
  # Reorder levels for chr and chr_peak
  df_protein$chr <- factor(df_protein$chr, levels = CHR_NAMES)
  
  
  # Check the new levels
  #levels(df_protein$Chromatin_State)
  #levels(df_protein$chr)
  
  
  
  #  Print the structure to verify
  #str(df_protein)
  
  saveRDS(object = df_protein,
          file = file.path(output_folder, protein, paste0(protein, ".rds")))
  
  
  # Loop over unique motifs
  unique_motifs <- unique(df_protein$motif)
  
  for (motif in unique_motifs) {
    # Filter the data for the current motif
    # debug motif <- unique_motifs[1]
    print(motif)
    df_motif <- df_protein %>% filter(motif == !!motif)
    
    # # Get the unique experiment_ids and biosamples for the current motif
    # experiment_biosample <- df_motif %>%
    #   select(experiment_id, biosample) %>%
    #   distinct()
    #
    # # Create a list of CpG positions for each experiment_id
    # cpg_list <- list()
    # for (exp_id in unique(experiment_biosample$experiment_id)) {
    #   cpg_positions <- df_motif %>%
    #     filter(experiment_id == !!exp_id) %>%
    #     select(chr, start, end) %>%
    #     mutate(position = paste(chr, start, end, sep = ":")) %>%
    #     pull(position)  # Get the positions as a vector
    #   cpg_list[[exp_id]] <- unique(cpg_positions)  # Add to list
    # }
    #
    # # Set the alpha based on the number of experiments per biosample
    # biosample_alpha <- experiment_biosample %>%
    #   group_by(biosample) %>%
    #   mutate(alpha = ifelse(n() > 1, 0.5, 1))  # Multiple experiments have transparency
    #
    # # Convert biosample_alpha$alpha to a numeric vector (not list)
    # experiment_alpha <- biosample_alpha$alpha
    #
    # # Correctly map colors for each experiment based on the biosample
    # experiment_colors <- sapply(experiment_biosample$biosample, function(b) {
    #   group.colors[match(b, names(group.colors))]  # Ensure matching biosample and color order
    # })
    #
    # # Generate the Venn diagram for the current motif
    # venn_plot <- venn.diagram(
    #   x = cpg_list,
    #   category.names = names(cpg_list),
    #   filename = NULL,  # Don't save to a file yet
    #   output = TRUE,
    #   fill = experiment_colors,  # Use the biosample colors
    #   alpha = as.numeric(experiment_alpha),  # Ensure alpha is numeric
    #   cat.col = experiment_colors,  # Color the labels by biosample
    #   cat.cex = 1.5,
    #   cex = 1.5,
    #   main = paste("CpG Positions for Motif:", motif)
    # )
    #
    # # Save the Venn diagram
    # venn_filename <- file.path(output_folder,protein, paste0("VennDiagram_", motif, ".png"))
    # png(venn_filename, width = 800, height = 800)
    # grid.draw(venn_plot)
    # dev.off()
    #
    # # Optional: Print a message to indicate that the plot is saved
    # print(paste("Venn diagram saved for motif:", motif))
    
  }
  
  
}

# Combine the list into a single data frame
summary_df <- do.call(rbind, summary_list)

# Display the summary DataFrame
print(dim(summary_df))

saveRDS(object = summary_df,file = file.path(output_folder,"summary_df.rds"))
rm(summary_list)
gc()

# Calculate the ratio of n_peak_no_motif_no_CpG to n_peaks
summary_df$no_motif_no_CpG_ratio <- summary_df$n_peak_no_motif_no_CpG / summary_df$n_peaks

# Plot the histogram
histogram_plot <- ggplot(summary_df, aes(x = no_motif_no_CpG_ratio)) +
  geom_histogram(binwidth = 0.05, fill = "#00BFC4", color = "black", alpha = 0.7) +
  labs(
    title = "Distribution of No Motif No CpG Peaks Ratio",
    x = "Ratio of No Motif No CpG Peaks",
    y = "Count"
  ) +
  theme_minimal()

# Save the absolute plot
ggsave(filename = file.path(output_folder, "histogram_plot_fractions.png"),
       plot = histogram_plot, width = 10, height = 6)

# Loop over protein folders
for (protein in protein_folders) {
  # protein <- "CTCF"
  print(protein)
  
  # Create output directory for each protein if it doesn't exist
  protein_output_dir <- file.path(output_folder, protein)
  if (!dir.exists(protein_output_dir)) {
    dir.create(protein_output_dir, recursive = TRUE)
  }
  
  # Filter the data for the specific protein
  
  protein_data <- summary_df[summary_df$protein == protein, ]
  
  # Create a new column for the ratio of n_peak_no_motif_no_CpG to n_peaks
  protein_data$no_motif_no_CpG_ratio <- protein_data$n_peak_no_motif_no_CpG / protein_data$n_peaks
  
  # Create the bar plot
  bar_plot <- ggplot(protein_data, aes(x = interaction(biosample, experiment_id, motif), 
                                       y = no_motif_no_CpG_ratio, 
                                       fill = biosample)) +
    geom_bar(stat = "identity") +
    labs(
      title = paste("Ratio of Peaks with No Motif and No CpG for", protein),
      x = "Biosample x Experiment ID",
      y = "Ratio of No Motif No CpG Peaks"
    ) +
    scale_fill_manual(values = group.colors) +  # Use the custom biosample colors
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
  
  # Print the plot
  print(bar_plot)
  
  ########### continue 
  
  # Now run the code
  summary_df_long <- summary_df %>%
    mutate(matched_peaks = n_peaks - n_peaks_unmatched) %>%  # Calculate matched peaks
    select(protein, biosample, experiment_id, motif, n_peaks, n_peaks_unmatched, matched_peaks) %>%
    pivot_longer(cols = c(matched_peaks, n_peaks_unmatched), 
                 names_to = "peak_type", 
                 values_to = "peak_count")  # Reshape data to long format for stacked bars
  
  # Create a new variable that combines experiment_id and motif
  summary_df_long <- summary_df_long %>%
    mutate(experiment_motif = paste(biosample,experiment_id, motif, sep = "_"))
  
  
  # Generate plot for absolute number of peaks
  p_absolute <- ggplot(data = summary_df_long[summary_df_long$protein == protein,], 
                       aes(x = experiment_motif, y = peak_count, fill = biosample, alpha = peak_type)) +
    geom_bar(stat = "identity", position = "stack") +  # Stack matched and unmatched peaks
    scale_alpha_manual(values = c(matched_peaks = 1, n_peaks_unmatched = 0.5)) +  # Different alpha for unmatched peaks
    scale_fill_manual(values = group.colors) +  # Use the custom color palette
    labs(title = paste("Number of Peaks per Biosample for", protein),
         x = "Experiment ID and Motif", 
         y = "Number of Peaks") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for readability
  
  # Save the absolute plot
  ggsave(filename = file.path(protein_output_dir, paste0(protein, "_absolute_peaks.png")),
         plot = p_absolute, width = 10, height = 6)
  
  # Calculate percentage for each row in summary_df_long
  summary_df_long <- summary_df_long %>%
    group_by(experiment_motif, biosample) %>%  # Group by experiment_motif and biosample
    mutate(percentage = peak_count / sum(peak_count) * 100) %>%  # Calculate percentage
    ungroup()
  
  # Generate plot for percentage of peaks
  p_percentage <- ggplot(data = summary_df_long[summary_df_long$protein == protein,], 
                         aes(x = experiment_motif, y = percentage, fill = biosample, alpha = peak_type)) +
    geom_bar(stat = "identity", position = "stack") +  # Stack matched and unmatched peaks
    scale_alpha_manual(values = c(matched_peaks = 1, n_peaks_unmatched = 0.5)) +  # Different alpha for unmatched peaks
    scale_fill_manual(values = group.colors) +  # Use the custom color palette
    labs(title = paste("Percentage of Peaks per Biosample for", protein),
         x = "Experiment ID and Motif", 
         y = "Percentage of Peaks") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))  # Rotate x-axis labels for readability
  
  # Save the percentage plot
  ggsave(filename = file.path(protein_output_dir, paste0(protein, "_percentage_peaks.png")),
         plot = p_percentage, width = 10, height = 6)
}

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
library(ComplexUpset)
library(ggplot2)
library(reshape2)

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
  
  # Create a new column for the ratio of n_peak_no_motif_no_CpG to n_peaks
  protein_data$no_motif_no_CpG_ratio <- protein_data$n_peak_no_motif_no_CpG / protein_data$n_peaks
  
  df_protein <- readRDS(file = file.path(output_folder, protein, paste0(protein, ".rds")))
  
  # Loop over unique motifs
  unique_motifs <- unique(df_protein$motif)
  
  for (motif in unique_motifs) {
    # Filter the data for the current motif
    # debug motif <- unique_motifs[1]
    print(motif)
    df_motif <- df_protein %>% filter(motif == !!motif)
    motif_data <- protein_data[protein_data$motif == motif,]
    
    # Generate all combinations of two biosamples
    biosample_combinations <- combn(unique(motif_data$biosample), 2, simplify = FALSE)
    
    # Loop through each combination of two biosamples
    for (biosample_comb in biosample_combinations) {
      
      # Subset df_protein to contain only the rows for the two selected biosamples
      subset_df_protein <- df_motif %>%
        filter(biosample %in% biosample_comb)
      
      # Add a column to check if the combination of (chr, start_cg, end_cg) is in one, the other, or both biosamples
      subset_df_protein <- subset_df_protein %>%
        group_by(chr, start_cg, end_cg) %>%
        mutate(
          presence_in_biosamples = case_when(
            n_distinct(biosample) == 1 ~ biosample[1],  # If only one biosample, use its name
            n_distinct(biosample) == 2 ~ "both"         # If present in both, use "both"
          )
        ) %>%
        ungroup()  # Ensure to ungroup after mutation
      
      # Ensure there are at least two biosamples in the subset before processing
      if (nrow(subset_df_protein) >= 2) {
        # Print or process the subset dataframe for the current combination
        cat("Processing combination of biosamples:", paste(biosample_comb, collapse = ", "), "\n")
        print(subset_df_protein)
        
        # Add any further analysis or plotting here, for example:
        # You can summarize, create visualizations, or export the subset
      }
    }
    
    
  }
}
