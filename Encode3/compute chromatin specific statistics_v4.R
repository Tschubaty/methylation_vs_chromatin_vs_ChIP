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
##  v_x 16.12.2024
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
# reads WGBS data
combined_methylation_data_compact <- readRDS(file.path(input_WGBS_dir,"combined_methylation_data_compact.rds"))

####################### make WGBS plots and save in list ###########################
if(FALSE){

# Initialize an empty list to store plots
plot_list <- list()

for (state in names(chromatin_state_colors_short)) {
  print(state)
  
  # Initialize a list to store plots for this state
  plot_list[[state]] <- list()
  
  # Loop through each biosample and generate plots
  for (biosample in names(group.colors)) {
    print(biosample)
    
    # Access the correct columns for Chromatin State and Methylation
    chromatin_col <- paste0("Chromatin_State_", biosample)
    methylation_col <- paste0("fRead_", biosample)
    
    # Filter data for the current chromatin state
    state_data <- combined_methylation_data_compact %>%
      filter(!!sym(chromatin_col) == state) %>%
      filter(!is.na(!!sym(methylation_col))) %>%
      select(!!sym(methylation_col)) %>%
      rename(Methylation = !!sym(methylation_col))
    
    # Histogram for Methylation Values
    histogram_plot <- ggplot(state_data, aes(x = Methylation)) +
      geom_histogram(
        alpha = 0.7,
        position = "identity",
        breaks = seq(0, 1, length.out = 101),
        fill = chromatin_state_colors_short[state]
      ) +
      labs(
        title = paste(
          "Histogram of Methylation of ",
          biosample,
          " for Chromatin State:",
          state,
          "in Biosample:",
          biosample
        ),
        x = "Methylation Level",
        y = "Frequency"
      ) +
      theme_minimal()
    
    # ECDF for Methylation Values
    ecdf_plot <- ggplot(state_data,
                        aes(x = Methylation, color = chromatin_state_colors_short[state])) +
      stat_ecdf(geom = "step", size = 1.2) +
      labs(
        title = paste(
          "ECDF of Methylation of " ,
          biosample,
          " in Chromatin State:",
          state
        ),
        x = "Methylation Level",
        y = "ECDF",
        color = state
      ) +
      theme_minimal()
    
    # # Save the plots
    #ggsave(file.path(output_dir_plots, paste0("histogram_",biosample,"_", state, ".png")), histogram_plot, width = 10, height = 6)
    #ggsave(file.path(output_dir_plots, paste0("ecdf_", biosample,"_",state, ".png")), ecdf_plot, width = 10, height = 6)
    #
    # Optional: Display plots in RStudio viewer
    # print(histogram_plot)
    # print(ecdf_plot)
    
    # Save plots into the list with a clear naming structure
    plot_list[[state]][[biosample]] <- list(Histogram = histogram_plot, ECDF = ecdf_plot)
  }
}

saveRDS(object = plot_list,file = file.path(out_dir,"plot_list,rds"))
rm(plot_list)
rm(ecdf_plot)
rm(final_grid_plot)
rm(histogram_plot)
rm(combined_plot)
rm(state_data)
gc()
}
####################### make grid plots and save plots  ###########################
if(FALSE){
library(ggplot2)
library(dplyr)
library(gridExtra) # For combining plots

# Loop through each chromatin state
for (state in names(chromatin_state_colors_short)) {
  print(state)
  
  # List to store plots for the current chromatin state
  plot_list <- list()
  
  # Loop through each biosample and generate plots
  for (biosample in names(group.colors)) {
    print(biosample)
    
    # Access the correct columns for Chromatin State and Methylation
    chromatin_col <- paste0("Chromatin_State_", biosample)
    methylation_col <- paste0("fRead_", biosample)
    
    # Filter data for the current chromatin state
    state_data <- combined_methylation_data_compact %>%
      filter(!!sym(chromatin_col) == state) %>%
      filter(!is.na(!!sym(methylation_col))) %>%
      select(!!sym(methylation_col)) %>%
      rename(Methylation = !!sym(methylation_col))
    
    # Histogram for Methylation Values
    histogram_plot <- ggplot(state_data, aes(x = Methylation)) +
      geom_histogram(
        alpha = 0.7,
        position = "identity",
        breaks = seq(0, 1, length.out = 101),
        fill = chromatin_state_colors_short[state]
      ) +
      labs(
        title = paste("Histogram:", biosample),
        x = "Methylation Level",
        y = "Frequency"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(size = 10))
    
    # ECDF for Methylation Values
    ecdf_plot <- ggplot(state_data, aes(x = Methylation)) +
      stat_ecdf(geom = "step",
                size = 1.2,
                color = chromatin_state_colors_short[state]) +
      labs(title = paste("ECDF:", biosample),
           x = "Methylation Level",
           y = "ECDF") +
      theme_minimal() +
      theme(plot.title = element_text(size = 10))
    
    # Combine histogram and ECDF into a vertical grid (2 rows)
    combined_plot <- gridExtra::grid.arrange(histogram_plot, ecdf_plot, nrow = 2)
    
    # Append to plot list
    plot_list[[biosample]] <- combined_plot
  }
  
  # Combine all biosample plots into a 2x4 grid
  final_grid_plot <- gridExtra::grid.arrange(
    grobs = plot_list,
    ncol = 4,
    # 4 biosamples side-by-side
    top = paste("Chromatin State:", state)
  )
  
  # Save the combined grid plot
  ggsave(
    file.path(
      output_dir_plots,
      paste0("combined_hist_ecdf_", state, ".png")
    ),
    final_grid_plot,
    width = 16,
    height = 10
  )
  
  print(paste("Saved grid plot for Chromatin State:", state))
}
}
###################### mark occupancy ###############


library(dplyr)

# Initialize empty columns for protein hits per biosample
combined_methylation_data_compact <- combined_methylation_data_compact %>%
  mutate(
    protein_hits_A549 = NA_character_,
    protein_hits_GM12878 = NA_character_,
    protein_hits_HepG2 = NA_character_,
    protein_hits_K562 = NA_character_
  )

# List of biosamples and their respective ChIP directories
biosamples <- c("A549", "GM12878", "HepG2", "K562")

# Iterate through all proteins and their respective experiments
proteins <- list.files(path = input_ChIP_dir)

for (protein in proteins) {
  print(paste("Processing protein:", protein))
  
  # Skip log folder
  if (protein == "log") next
  
  # Get all files for the protein
  protein_dir <- file.path(input_ChIP_dir, protein)
  experiment_files <- list.files(path = protein_dir, full.names = TRUE)
  
  for (file_path in experiment_files) {
    cat("\nStarting experiment file:", file_path, "\n")
    start_total <- Sys.time()
    
    # Extract experiment details
    file_name <- basename(file_path)
    experiment_id <- strsplit(file_name, "_")[[1]][3]
    biosample <- strsplit(file_name, "_")[[1]][1]
    motif <- strsplit(file_name, "_")[[1]][4]
    
    # Skip if biosample is not in the list
    if (!(biosample %in% biosamples)) next
    
    # Read ChIP data
    start_read <- Sys.time()
    df_ChIP <- read_bed_file(file_path)
    end_read <- Sys.time()
    print(paste("Data loaded in:", round(difftime(end_read, start_read, units = "secs"), 2), "seconds"))
    
    # Identify matching CpGs in combined_methylation_data_compact
    start_match <- Sys.time()
    match_idx <- which(
      combined_methylation_data_compact$chr %in% df_ChIP$chr &
        combined_methylation_data_compact$Start %in% df_ChIP$start_cg &
        combined_methylation_data_compact$End %in% df_ChIP$end_cg
    )
    end_match <- Sys.time()
    print(paste("Index matching completed in:", round(difftime(end_match, start_match, units = "secs"), 2), "seconds"))
    
    # Add protein information to the respective biosample column
    start_update <- Sys.time()
    if (length(match_idx) > 0) {
      column_name <- paste0("protein_hits_", biosample)

      # # Use mapply for row-wise updates
      # combined_methylation_data_compact[match_idx, column_name] <- mapply(
      #   function(current_value, new_value) {
      #     if (is.na(current_value)) {
      #       new_value
      #     } else {
      #       paste(current_value, new_value, sep = ";")
      #     }
      #   },
      #   current_value = combined_methylation_data_compact[match_idx, column_name],
      #   new_value = paste0(protein, "_", experiment_id)
      # )
      
      # # Perform row-wise updates using Map (element-wise)
      # combined_methylation_data_compact[match_idx, column_name] <- Map(
      #   function(current_value, new_value) {
      #     if (is.na(current_value)) {
      #       new_value  # Insert new value if NA
      #     } else {
      #       paste(current_value, new_value, sep = ";")  # Append new value with semicolon
      #     }
      #   },
      #   current_value = combined_methylation_data_compact[match_idx, column_name],
      #   new_value = rep(paste0(protein, "_", experiment_id), length(match_idx))
      # )
      
      #Update the column using vapply
      combined_methylation_data_compact[match_idx, column_name] <- vapply(
        match_idx,
        function(i) {
          current_value <- combined_methylation_data_compact[i, column_name]
          new_value <- paste0(protein, "_", experiment_id,"_",motif)
          
          if (is.na(current_value)) {
            new_value  # Insert new value if NA
          } else {
            paste(current_value, new_value, sep = ";")  # Append new value
          }
        },
        FUN.VALUE = character(1)  # Ensure a character vector is returned
      )

    }
    end_update <- Sys.time()
    print(paste("Protein information added in:", round(difftime(end_update, start_update, units = "secs"), 2), "seconds"))
    
    # Total time for this iteration
    end_total <- Sys.time()
    print(paste("Total time for this file:", round(difftime(end_total, start_total, units = "secs"), 2), "seconds"))
  }
}


# Final output
print(head(combined_methylation_data_compact))
saveRDS(object = combined_methylation_data_compact,file = file.path(out_dir,"combined_methylation_data_compact_withpeaks.rds"))
test <- head(combined_methylation_data_compact[!is.na(combined_methylation_data_compact$protein_hits_GM12878),])
#################################################################################

# List of biosamples
biosamples <- c("A549", "GM12878", "HepG2", "K562")

# Iterate through each biosample
for (biosample in biosamples) {
  print(paste("Processing biosample:", biosample))
  
  # Identify relevant columns
  protein_hits_col <- paste0("protein_hits_", biosample)
  methylation_col <- paste0("fRead_", biosample)
  chromatin_state_col <- paste0("Chromatin_State_", biosample)
  
  # Subset only relevant columns and create 'is_occupied' column
  biosample_methylation_data <- combined_methylation_data_compact %>%
    select(!!sym(methylation_col), !!sym(protein_hits_col), !!sym(chromatin_state_col)) %>%
    rename(
      Methylation = !!sym(methylation_col),
      ProteinHits = !!sym(protein_hits_col),
      ChromatinState = !!sym(chromatin_state_col)
    ) %>%
    mutate(is_occupied = !is.na(ProteinHits))  # Logical for protein occupancy
  
  # Print summary of the processed data
  print(paste("Total CpGs in", biosample, ":", nrow(biosample_methylation_data)))
  print(paste("Occupied CpGs:", sum(biosample_methylation_data$is_occupied)))
  print(paste("Unoccupied CpGs:", sum(!biosample_methylation_data$is_occupied)))
  
  # Histogram plot
  histogram_plot <- ggplot(biosample_methylation_data, aes(x = Methylation, fill = is_occupied)) +
    geom_histogram(bins = 100, alpha = 0.7, position = "identity") +
    scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "grey")) +
    labs(
      title = paste("Methylation Histogram -", biosample),
      x = "Methylation Level",
      y = "Frequency",
      fill = "Occupancy"
    ) +
    theme_minimal()
  
  # Save histogram plot
  histogram_file <- paste0("occupancy_histogram_methylation_", biosample, ".png")
  ggsave(filename = file.path(output_dir_plots,histogram_file), plot = histogram_plot, width = 10, height = 6)
  print(paste("Histogram plot saved:", histogram_file))
  
  # Summarize data for the bar plot
  plot_data <- biosample_methylation_data %>%
    group_by(ChromatinState, is_occupied) %>%
    summarize(Count = n(), .groups = "drop") %>%
    mutate(is_occupied = ifelse(is_occupied, "Occupied", "Unoccupied"))
  
  # Bar plot
  bar_plot <- ggplot(plot_data, aes(x = ChromatinState, y = Count, fill = is_occupied)) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.8) +
    scale_fill_manual(values = c("Occupied" = "blue", "Unoccupied" = "grey")) +
    labs(
      title = paste("CpG Occupancy in Chromatin States -", biosample),
      x = "Chromatin State",
      y = "Number of CpGs",
      fill = "Occupancy"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
      panel.grid.major.x = element_blank()  # Remove vertical grid lines
    )
  
  # Save bar plot
  bar_plot_file <- paste0("occupancy_barplot_chromatin_state_", biosample, ".png")
  ggsave(filename = file.path(output_dir_plots,bar_plot_file), plot = bar_plot, width = 10, height = 6)
  print(paste("Bar plot saved:", bar_plot_file))
}


library(dplyr)
library(stringr)

for (protein in proteins) {
  print(paste("Processing protein:", protein))
  
  # Skip log folder
  if (protein == "log") next
  
  # Get all files for the protein
  protein_dir <- file.path(input_ChIP_dir, protein)
  experiment_files <- list.files(path = protein_dir, full.names = TRUE)
  
  for (file_path in experiment_files) {
    cat("\nStarting experiment file:", file_path, "\n")
    start_total <- Sys.time()
    
    # Extract experiment details
    file_name <- basename(file_path)
    experiment_id <- strsplit(file_name, "_")[[1]][3]
    biosample <- strsplit(file_name, "_")[[1]][1]
    motif <- strsplit(file_name, "_")[[1]][4]
    
    # Identify relevant columns
    protein_hits_col <- paste0("protein_hits_", biosample)
    methylation_col <- paste0("fRead_", biosample)
    chromatin_state_col <- paste0("Chromatin_State_", biosample)
    
    # Name for the current peak
    peak_name <- paste0(protein, "_", experiment_id, "_", motif)
    
    # Process biosample-specific data
    biosample_methylation_data <- combined_methylation_data_compact %>%
      select(!!sym(methylation_col), !!sym(protein_hits_col), !!sym(chromatin_state_col)) %>%
      rename(
        Methylation = !!sym(methylation_col),
        ProteinHits = !!sym(protein_hits_col),
        ChromatinState = !!sym(chromatin_state_col)
      ) %>%
      mutate(
        ProteinHits = case_when(
          is.na(ProteinHits) ~ "not occupied",  # If NA, not occupied
          str_detect(ProteinHits, fixed(peak_name)) ~ "peak_name",  # If peak_name matches
          TRUE ~ "other protein"  # Otherwise, other proteins
        ),
        is_occupied = ProteinHits != "not occupied"  # Logical for occupancy
      )
    
    # Define the output subfolder for the current protein
    protein_subfolder <- file.path(out_dir, "biosample_methylation_data", protein)
    if (!dir.exists(protein_subfolder)) {
      dir.create(protein_subfolder, recursive = TRUE)
    }
    
    # Create the filename with experiment_id, biosample, and motif
    output_file <- file.path(protein_subfolder, paste0("biosample_methylation_data_", biosample, "_", experiment_id, "_", motif, ".rds"))
    
    # Save the data as an RDS file
    #saveRDS(biosample_methylation_data, output_file)
    
    # Log message
    print(paste("Data saved for biosample:", biosample, "experiment:", experiment_id, "motif:", motif, "at", output_file))
    
    
    # Total time for this iteration
    end_total <- Sys.time()
    print(paste("Total time for this file:", round(difftime(end_total, start_total, units = "secs"), 2), "seconds"))
  }
}



# # details unknown categtory 
# permutation_test <- function(file_path, biosample) {
#   df_ChIP <- read_bed_file(file_path)
#   
#   # Extract the column for the biosample's chromatin state
#   chromatin_col <- paste0("Chromatin_State_", biosample)
#   
#   # Categorize rows based on start_motif, start_cg, and chromatin state
#   df_ChIP <- df_ChIP %>%
#     mutate(
#       !!sym(chromatin_col) := case_when(
#         is.na(!!sym(chromatin_col)) | !!sym(chromatin_col) == "." & start_motif == -1 ~ "No motif in peak",
#         is.na(!!sym(chromatin_col)) | !!sym(chromatin_col) == "." & start_motif != -1 & start_cg == 0 ~ "No CG in motif",
#         is.na(!!sym(chromatin_col)) | !!sym(chromatin_col) == "." & start_motif != -1 & start_cg != 0 ~ "No State assignment for CG",
#         TRUE ~ !!sym(chromatin_col)  # Keep original value otherwise
#       )
#     )
#   
# 
#   # return( )
# }

if(FALSE){
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
    print(file_path)
    # Extract experiment details from file name
    file_name <- basename(file_path)
    experiment_id <- strsplit(file_name, "_")[[1]][3]
    biosample <- strsplit(file_name, "_")[[1]][1]
    
    df_ChIP <- read_bed_file(file_path)
    
    # Extract the column for the biosample's chromatin state
    chromatin_col <- paste0("Chromatin_State_", biosample)
    

    # complete plot here make histogrammns and ECDF of df_ChIP & combined_methylation_data_compact for each chromatrin state 
    
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
output_csv <- file.path(out_dir, "chromatin_state_tests.csv")
write_csv(final_summary, output_csv)
}
###############################################################################
if(FALSE){

# Loop through each protein folder
for (protein_file in protein_files[1:10]) {
  # protein_file in protein_files[1:10]
  protein_name <- gsub(pattern = ".fimo_methylation_peaks.bed",replacement = "",x = basename(protein_file))  # Extract protein name
  df<- read.delim(header = FALSE,file = file.path(input_folder,protein_file),col.names = colnames_df)
  head(df)
  
  ChIP_files <- list.files(path = file.path(ChIP.seq_folder,protein_name),pattern = "*.bed")
  # for(ChIP_file in ChIP_files)
  # df_chip <- read.delim(file = 
  
  }
  
  
  bed_files <- list.files(protein_folder, pattern = "*_methylation.chromatinstate.fimo.bed", full.names = TRUE)
  print(protein_name)
  # Extract sample types for each file
  sample_types <- sapply(bed_files, extract_sample_type)
  
  # Generate all pairwise combinations of the sample types
  sample_pairs <- combn(sample_types, 2, simplify = FALSE)
  
  # # Process each .bed file
  # for (bed_file in bed_files) {
  #   print(basename(bed_file))
  #   #df <- read_bed_file(bed_file)
  #   
  #   # Perform further analysis on the dataframe `df`
  #   #print(head(df))
  # }
  
  for (pair in sample_pairs) {
    sample_1 <- pair[[1]]
    sample_2 <- pair[[2]]
    
    # Find the corresponding files for the two samples
    file_1 <- bed_files[grep(sample_1, bed_files)]
    file_2 <- bed_files[grep(sample_2, bed_files)]
    
    print(sprintf("Comparing samples: %s vs %s", sample_1, sample_2))
    
    # Read the data for both samples
    df_1 <- read_bed_file(file_1)
    df_2 <- read_bed_file(file_2)
    
    df_merged <- full_outer_join_samples(df_1, df_2, sample_1, sample_2)
    
    # Print the first few rows of the merged dataframe to verify
    print(head(df_merged))
    
  }
}



# df_rna <-
#   read.csv(file = file.path(RNA_seq_folder, "rna_seq_of_chip_sep_tartgets_summery.csv"))
# files <-
#   list.files(path = input_folder, pattern = "*.ChiP_vs_WGBS.*")
# 
# files <- file.path(input_folder , files)


# #debug files <- list.files(path = input_folder, pattern = "*random*")
# files[length(files) + 1] <-
#   "monte_carlo_results/data/script4.random.peaks_200bp_apart.RDS"
# 
# #debug
# files <- files[c(
#   which(grepl(x = files, pattern = "AGO1"))[1],
#   which(grepl(x = files, pattern = "CEBPB"))[1],
#   which(grepl(x = files, pattern = "IKZF1"))[1],
#   which(grepl(x = files, pattern = "MAFF"))[1],
#   which(grepl(x = files, pattern = "MAFK"))[1],
#   which(grepl(x = files, pattern = "NBN"))[1],
#   which(grepl(x = files, pattern = "PCBP2"))[1],
#   which(grepl(x = files, pattern = "PRPF4"))[1]
# )]

# files <- files[c(
#   which(grepl(x = files, pattern = "ZBTB33"))[1],
#   which(grepl(x = files, pattern = "AGO1"))[1],
#   which(grepl(x = files, pattern = "MAFF"))[1],
#   which(grepl(x = files, pattern = "PCBP2"))[1],
#   which(grepl(x = files, pattern = "CEBPB"))[1],
#   which(grepl(x = files, pattern = "MAFK"))[1],
#   which(grepl(x = files, pattern = "PRPF4"))[1],
#   which(grepl(x = files, pattern = "JUN"))[1],
#   which(grepl(x = files, pattern = "ATF1"))[1],
#   which(grepl(x = files, pattern = "YY1"))[1],
#   which(grepl(x = files, pattern = "CTCF"))[1]
# )]

if(FALSE){

processed_results_summery <- data.frame()


#n_cores <- 2
# path_file <- file.path(output_folder,sprintf("output_log start %s.txt",start_script))
# cl <- makeCluster(n_cores,outfile=path_file)
registerDoParallel(n_cores)

processed_results_summery <-
  foreach (f = 1:length(files), .combine = rbind) %dopar% {
    #for (f in 1:nrow(files)) {
    # f <- which(grepl(x = files,pattern = "ZBTB33"))[1]
    #f <- which(grepl(x = files,pattern = "GPBP1L1"))
    # f <- 1
    # f <- which(grepl(x = files,pattern = "ZZZ3"))[1]
    
    ## START execution
    start_time <- Sys.time()
    file_name <- files[f]
    protein_name <-
      gsub(pattern = "4. chromatin_annotated_v03/|.ChiP_vs_WGBS_vs_chromatin.RDS|monte_carlo_results/data/script4.|.peaks_200bp_apart.RDS|",
           replacement = "",
           x = file_name)
    
    # create folder
    target_folder <- file.path(output_folder, "targets")
    dir.create(file.path(output_folder, "targets"), showWarnings = FALSE)
    # create protein folder
    output_folder_protein <-
      file.path(target_folder, protein_name)
    dir.create(output_folder_protein, showWarnings = FALSE)
    # create plot folder
    output_folder_protein_plots <-
      file.path(output_folder_protein, "plots")
    dir.create(output_folder_protein_plots, showWarnings = FALSE)
    # create data folder
    output_folder_protein_data <-
      file.path(output_folder_protein, "data")
    dir.create(output_folder_protein_data, showWarnings = FALSE)
    
    summery_file_name <- file.path(output_folder_protein_data,
                                   paste("summery",
                                         protein_name,
                                         "RDS",
                                         sep = "."))
    
    if (file.exists(summery_file_name) & cluster) {
      df_temp <- readRDS(summery_file_name)
    } else{
      # read input
      df <- readRDS(file = file_name)
      # allocate output
      df_temp <- data.frame()
      
      ############# general plots #################
      
      
      # look only at data with same chromatin state
      df <-
        df[df$chromatin_state_HepG2 == df$chromatin_state_K562,]
      df$chromatin_state <-
        df$chromatin_state_HepG2
      df <-
        df[, !(names(df) %in% c("chromatin_state_HepG2" , "chromatin_state_K562"))]
      
      # look only at data with methylation values
      df <-
        df[!is.na(df$fRead_K562) & !is.na(df$fRead_HepG2) , ]
      #  restrict to unique peaks
      merge_colnames <-
        c(
          "chr",
          "merge_peak",
          "ChiP_state",
          "mRead_K562",
          "mRead_HepG2",
          "fRead_K562",
          "fRead_HepG2",
          "mean_n_K562",
          "mean_n_HepG2" ,
          "n_CPG" ,
          "chromatin_state",
          "hits",
          "name",
          "nRead_K562",
          "nRead_HepG2"
        )
      
      df <- df[!duplicated(df[, merge_colnames]), merge_colnames]
      
      
      # start plotting
      
      # ven diagram
      
      
      if (!cluster) {
        # if no data or random data
        if (protein_name == "random") {
          plot_ven <- ggplot2::ggplot()
        } else{
          ven_values <-  c(
            "HepG2" = sum(df$ChiP_state == "HepG2"),
            "K562" = sum(df$ChiP_state == "K562"),
            "K562&HepG2" = sum(df$ChiP_state == "both")
          )
          plot_ven <- plot(
            eulerr::euler(ven_values),
            fills = list(fill = group.colors),
            legend = list(side = "right"),
            quantities = list(cex = 2),
            labels = c()
          )
          
          ggplot2::ggsave(
            filename = paste(
              protein_name,
              "Venn_Diagramm_peaks",
              picuture_file_extension,
              sep = "."
            ),
            plot = plot_ven,
            device = picuture_file_extension,
            path =  output_folder_protein_plots,
            width = 1920,
            height = 1080,
            units = "px"
          )
        }
        
        # plot bio replicates
        plot_bio_rep <-
          ggplot2::ggplot(data = df,
                          mapping = ggplot2::aes(x = hits, fill = ChiP_state)) +
          ggplot2::geom_bar() +
          ggplot2::ggtitle(paste("peak classification according to sample: ", protein_name)) +
          ggplot2::scale_fill_manual(name = "ChiP_state",
                                     values = group.colors) +
          ggplot2::xlab("# replicated") +
          ggplot2::theme(
            legend.key = ggplot2::element_rect(fill = NA),
            text = ggplot2::element_text(size = 15),
            axis.line = ggplot2::element_line(colour = "black"),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
          )
        
        ggplot2::ggsave(
          filename = paste(
            protein_name,
            "peak classification according to sample",
            picuture_file_extension,
            sep = "."
          ),
          plot = plot_bio_rep,
          device = picuture_file_extension,
          path =  output_folder_protein_plots,
          width = 1920,
          height = 1080,
          units = "px"
        )
        
        # extract all values
        rna_p <- df_rna[df_rna$gene_name == protein_name, ]
        #check if rna data present
        is_rna_data_present  <- nrow(rna_p) > 0
        if (is_rna_data_present) {
          ex <-
            unique(gsub(
              x =  grep(
                x = colnames(df_rna),
                pattern = "*._TPM",
                value = TRUE
              ),
              pattern = "_TPM.*",
              replacement = ""
            ))
          
          TPM <- as.vector(t(rna_p[paste(ex, "TPM", sep = "_")]))
          names(TPM) <- NULL
          TPM_ci_lower_bound <-
            as.vector(t(rna_p[paste(ex, "TPM_ci_lower_bound", sep = "_")]))
          names(TPM_ci_lower_bound) <- NULL
          TPM_ci_upper_bound <-
            as.vector(t(rna_p[paste(ex, "TPM_ci_upper_bound", sep = "_")]))
          names(TPM_ci_upper_bound) <- NULL
          FPKM = as.vector(t(rna_p[paste(ex, "FPKM", sep = "_")]))
          names(FPKM) <- NULL
          FPKM_ci_lower_bound <-
            as.vector(t(rna_p[paste(ex, "FPKM_ci_lower_bound", sep = "_")]))
          names(FPKM_ci_lower_bound) <- NULL
          FPKM_ci_upper_bound <-
            as.vector(t(rna_p[paste(ex, "FPKM_ci_upper_bound", sep = "_")]))
          names(FPKM_ci_upper_bound) <- NULL
          
          # build data_frame
          df_rna_p <- data.frame(
            experiment = ex,
            TPM = TPM,
            TPM_ci_lower_bound = TPM_ci_lower_bound,
            TPM_ci_upper_bound = TPM_ci_upper_bound,
            FPKM = FPKM,
            FPKM_ci_lower_bound = FPKM_ci_lower_bound,
            FPKM_ci_upper_bound = FPKM_ci_upper_bound
          )
          # plot_TPM <-
          #   ggplot2::ggplot(
          #     data = df_rna_p,
          #     mapping = ggplot2::aes(x = experiment, y = TPM, fill = experiment)
          #   ) +
          #   ggplot2::geom_bar(stat = "identity") +
          #   ggplot2::geom_errorbar(
          #     mapping = ggplot2::aes(ymin = TPM_ci_lower_bound, ymax = TPM_ci_upper_bound),
          #     width = 0.2,
          #     size = 1,
          #     color = "blue"
          #   ) +
          #   ggplot2::ggtitle(paste(protein_name, "expression in TPM"))
          
          
          plot_FPKM <-
            ggplot2::ggplot(
              data = df_rna_p,
              mapping = ggplot2::aes(
                x = experiment,
                y = FPKM,
                fill = experiment
              )
            ) +
            ggplot2::geom_bar(stat = "identity") +
            ggplot2::geom_errorbar(
              mapping = ggplot2::aes(ymin = FPKM_ci_lower_bound, ymax = FPKM_ci_upper_bound),
              width = 0.2,
              size = 1,
              color = "black"
            ) +
            ggplot2::ggtitle(paste(protein_name, "expression in FPKM")) +
            ggplot2::theme(
              legend.key = ggplot2::element_rect(fill = NA),
              text = ggplot2::element_text(size = 15),
              axis.line = ggplot2::element_line(colour = "black"),
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              panel.background = ggplot2::element_blank(),
            )
        } else{
          plot_TPM <- ggplot2::ggplot()
          plot_FPKM <- ggplot2::ggplot()
          
        }
        
        ggplot2::ggsave(
          filename = paste(
            protein_name,
            "exression in FPKM",
            picuture_file_extension,
            sep = "."
          ),
          plot = plot_FPKM,
          device = picuture_file_extension,
          path =  output_folder_protein_plots,
          width = 1920,
          height = 1080,
          units = "px"
        )
      }
      ##############################m START chromatin speciifc ########################
      
      plot_list_chromatin_state_specific <- list()
      
      for (state_counter in 1:length(chromatin_state_names)) {
        #dubug  state_counter <- 1
        
        c_state <- chromatin_state_names[state_counter]
        #print(c_state)
        # only look at one chromatine state
        df_state <-
          df[!is.na(df$chromatin_state == c_state) &
               df$chromatin_state == c_state, ]
        
        # if no data or random data
        if (nrow(df_state) < 1 | protein_name == "random") {
          empty_plot <- ggplot2::ggplot()
          
          plot_list_chromatin_state_specific <-
            append(plot_list_chromatin_state_specific,
                   list(empty_plot,
                        empty_plot))
          
          df_temp <- rbind(
            df_temp,
            data.frame(
              protein_name = protein_name,
              chromatin_state = c_state,
              n_state = nrow(df_state),
              p_love = NA,
              p_hate = NA,
              statistic  = NaN,
              n_Chip_both = 0,
              n_Chip_HepG2 = 0,
              n_Chip_K562 = 0
            )
          )
          next
        }
        
        sim_number <- 100000
        
        meth_statistic_vec <-
          function(ChiP_state,
                   fRead_HepG2,
                   fRead_K562) {
            mean_K562_statistic_value <-
              mean(fRead_HepG2[which(ChiP_state == "K562")] -  fRead_K562[which(ChiP_state == "K562")], na.rm = TRUE)
            mean_HepG2_statistic_value <-
              mean(fRead_HepG2[which(ChiP_state == "HepG2")] - fRead_K562[which(ChiP_state == "HepG2")], na.rm = TRUE)
            statistic <-
              mean_HepG2_statistic_value  - mean_K562_statistic_value
            return(statistic)
          }
        
        #real_stat <-  meth_statistic(df_state)
        real_stat <-
          meth_statistic_vec(df_state$ChiP_state,
                             df_state$fRead_HepG2,
                             df_state$fRead_K562)
        
        
        file_name_mote_carlo <-
          file.path(
            output_folder_protein_data,
            paste(protein_name,
                  c_state,
                  "monte_carlo_values",
                  "RDS", sep = ".")
          )
        
        if (file.exists(file_name_mote_carlo)) {
          monte_carlo_values <- readRDS(file = file_name_mote_carlo)
        } else{
          fRead_HepG2 <- df_state$fRead_HepG2
          fRead_K562 <- df_state$fRead_K562
          
          monte_carlo_values <- sapply(
            1:sim_number,
            FUN =  function(x) {
              return(
                meth_statistic_vec(
                  ChiP_state = sample(df_state$ChiP_state),
                  fRead_HepG2 = fRead_HepG2,
                  fRead_K562 = fRead_K562
                )
              )
            }
          )
          
          saveRDS(object = monte_carlo_values,
                  file = file.path(
                    output_folder_protein_data,
                    paste(
                      protein_name,
                      c_state,
                      "monte_carlo_values",
                      "RDS",
                      sep = "."
                    )
                  ))
        }
        
        p_love = sum(real_stat < monte_carlo_values) / length(monte_carlo_values)
        p_hate = sum(real_stat > monte_carlo_values) / length(monte_carlo_values)
        
        if (!cluster) {
          plot_meth_statistic <-
            ggplot2::ggplot(data = data.frame(statistic  = monte_carlo_values)) +
            ggplot2::geom_histogram(mapping = ggplot2::aes(x = statistic),
                                    bins = 1000) +
            ggplot2::geom_vline(xintercept = real_stat, color = "red") +
            # ggplot2::ggtitle(c_state)+
            #   paste(
            #     protein_name,
            #     ": " ,
            #     c_state,
            #     sim_number,
            #     "perm. p_affin:",
            #     format.pval(p_love, eps = 1 / sim_number),
            #     "  p_avers:",
            #     format.pval(p_hate, eps = 1 / sim_number)
            #   )
          # ) +
          ggplot2::xlab("methylation statistic") +
            ggplot2::theme(
              legend.key = ggplot2::element_rect(fill = NA),
              text = ggplot2::element_text(size = 10),
              axis.line = ggplot2::element_line(colour = "black"),
              axis.text = ggplot2::element_text(size = 15),
              axis.title  = ggplot2::element_text(size = 20),
              panel.background = ggplot2::element_blank(),
            )
          
          ggplot2::ggsave(
            filename = paste(
              protein_name,
              c_state,
              "permutation test",
              picuture_file_extension,
              sep = "."
            ),
            plot = plot_meth_statistic,
            device = picuture_file_extension,
            path =  output_folder_protein_plots,
            width = 1920,
            height = 1080,
            units = "px"
          )
          
          
          plot_fReads <- ggplot2::ggplot(data = df_state) +
            ggplot2::geom_point(mapping =
                                  ggplot2::aes(
                                    x = fRead_HepG2,
                                    y = fRead_K562,
                                    color = ChiP_state
                                  )) +
            ggplot2::scale_color_manual(name = "ChiP_state",
                                        values = group.colors) +
            ggplot2::ggtitle(c_state)+
            # ggplot2::ggtitle(label = paste(
            #   protein_name,
            #   c_state ,
            #   "HepG2:",
            #   sum(df_state$ChiP_state == "HepG2"),
            #   "K562:",
            #   sum(df_state$ChiP_state == "K562"),
            #   "both:",
            #   sum(df_state$ChiP_state == "HepG2")
            # )) +
            ggplot2::scale_x_continuous(limits = c(0, 1),
                                        name = expression('methylation_value'["Hepg2"]),) +
            ggplot2::scale_y_continuous(limits = c(0, 1),
                                        name = expression('methylation_value'["K562"]),) +
            ggplot2::theme(
              legend.key = ggplot2::element_rect(fill = NA),
              text = ggplot2::element_text(size = 12),
              axis.text.x = ggplot2::element_text(face = "bold"),
              #color = group.colors["Hepg2"]), #, size = 22),
              axis.title.y = ggplot2::element_text(face = "bold"),
              #color = group.colors["K562"]), # size = 32,
              axis.line = ggplot2::element_line(colour = "black"),
              panel.background = ggplot2::element_blank(),
              legend.position = "none"
            )
          
          ggplot2::ggsave(
            filename = paste(
              protein_name,
              c_state,
              "methylation values",
              picuture_file_extension,
              sep = "."
            ),
            plot = plot_fReads,
            device = picuture_file_extension,
            path =  output_folder_protein_plots,
            width = 1920,
            height = 1080,
            units = "px"
          )
          
          
          # plot_density_methylation_values <-
          #   ggplot2::ggplot(data = df_state,
          #                   mapping = ggplot2::aes(x = fRead_HepG2,
          #                                          y = fRead_K562)) +
          #   stat_density_2d(geom =  "polygon",#  "raster", # , "point"
          #                   aes(alpha = (..level..) ^ 2, fill = ChiP_state))+
          #   ggplot2::geom_point(mapping =
          #                         ggplot2::aes(x = fRead_HepG2,
          #                                      y = fRead_K562,
          #                                      color = ChiP_state))
          #
          #
          #   ggplot2::ggplot(data = df_state,
          #                   mapping = ggplot2::aes(x = fRead_HepG2,
          #                                          y = fRead_K562)) +
          
          #                   sum(df_state$ChiP_state == "K562"),
          #                   "both:",
          #                   sum(df_state$ChiP_state == "HepG2")
          #     ))+
          #   ggplot2::geom_point( mapping = ggplot2::aes(color = ChiP_state))+
          #   ggplot2::stat_density_2d(geom = "polygon",
          #   bins = 10,
          #   mapping = ggplot2::aes(alpha = ..level..,
          #                          fill = ChiP_state,
          #                          group = ChiP_state))
          #
          #
          #   # ggplot2::scale_color_manual(name = "ChiP_state",
          #   #                             values = group.colors)+
          #
          # set.seed(123)
          # plot_data <-
          #   data.frame(
          #     X = c(rnorm(300, 3, 2.5), rnorm(150, 7, 2)),
          #     Y = c(rnorm(300, 6, 2.5), rnorm(150, 2, 2)),
          #     Label = c(rep('A', 300), rep('B', 150))
          #   )
          # ggplot(data = df_state, aes(x = fRead_HepG2,
          #                             y = fRead_K562, group = ChiP_state)) +
          #   stat_density_2d(geom = "polygon",
          #                   aes(alpha = (..level..) ^ 2, fill = ChiP_state),
          #                   bins = 10)+
          #   geom_point( aes(color = ChiP_state))+
          #   xlim(0, 0.02)+
          #   ylim(0, 0.02)
          #
          #
          #
          #
          #
          #     ggplot2::stat_density_2d(geom = "polygon",
          #                   aes(alpha = ..level.., fill = ChiP_state),
          #                   bins = 4,group = df_state$ChiP_state)
          #
          # +
          # ggplot2::annotate(
          #   "text",
          #   x = 0.1,
          #   y = 1.05,
          #   label = paste(sum(df_state$ChiP_state == "K562"), "K562 datapoints"),
          #   size = 4,
          #   color = group.colors["K562"]
          # ) +
          # ggplot2::annotate(
          #   "text",
          #   x = 0.5,
          #   y = 1.05,
          #   label = paste(sum(df_state$ChiP_state == "both"), "both datapoints"),
          #   size = 4,
          #   color = group.colors["both"]
          # ) +
          # ggplot2::annotate(
          #   "text",
          #   x = 0.9,
          #   y = 1.05,
          #   label = paste(sum(df_state$ChiP_state == "HepG2"), "HepG2 datapoints"),
          #   size = 4,
          #   color = group.colors["HepG2"]
          #)
          
          # plot_islands_chip <-
          #   ggplot2::ggplot(data = df_state,
          #                   mapping = ggplot2::aes(x = name, fill = ChiP_state)) +
          #   ggplot2::geom_bar() +
          #   ggplot2::ggtitle(paste(protein_name,c_state, "CpG Island and Chip state"))
          #
          # plot_islands_hist <-
          #   ggplot2::ggplot(data = df_state,
          #                   mapping = ggplot2::aes(x = n_CPG , fill = name)) +
          #   ggplot2::geom_histogram(breaks = seq(from = 0, to = 100, by = 2)) +
          #   ggplot2::ggtitle(paste(protein_name,c_state, "CpG Island and nCpG"))
          #
          # plot_CpG_quality_K562 <- ggplot2::ggplot(data = df_state) +
          #   ggplot2::geom_point(mapping = ggplot2::aes(x = n_CPG, y = mean_n_K562, color = ChiP_state)) +
          #   ggplot2::ggtitle(paste(protein_name,c_state, "CpG_quality K562")) +
          #   ggplot2::ylim(0, 300)
          #
          #
          # plot_CpG_quality_HepG2 <- ggplot2::ggplot(data = df_state) +
          #   ggplot2::geom_point(mapping = ggplot2::aes(x = n_CPG, y = mean_n_HepG2, color = ChiP_state)) +
          #   ggplot2::ggtitle(paste(protein_name,c_state, "CpG_quality HepG2")) +
          #   ggplot2::ylim(0, 300)
          #
          # plot_methylation_read_depth <-  ggplot2::ggplot(data = df_state) +
          #   ggplot2::geom_point(mapping = ggplot2::aes(x = fRead_HepG2, y = fRead_K562, color = n_CPG)) +
          #   #ggplot2::scale_colour_gradient(low = "yellow", high = "red", na.value = NA)+
          #   ggplot2::scale_colour_gradient2(
          #     #breaks = c(0,5,10,20),
          #     low = "red",
          #     mid = "white",
          #     high = "blue",
          #     midpoint = 25,
          #     space = "Lab",
          #     na.value = NA,
          #     guide = "colourbar",
          #     aesthetics = "colour"
          #   ) +
          #   ggplot2::ggtitle(paste(protein_name,c_state, "methylation and read_depth"))
          #
          
          
          plot_list_chromatin_state_specific <-
            append(
              plot_list_chromatin_state_specific,
              list(plot_fReads, plot_meth_statistic)
            )
          
        } # end plotting
        
        
        df_temp <- rbind(
          df_temp,
          data.frame(
            protein_name = protein_name,
            chromatin_state = c_state,
            n_state = nrow(df_state),
            p_love = p_love,
            p_hate = p_hate,
            statistic  = real_stat,
            n_Chip_both = sum(df_state$ChiP_state == "both"),
            n_Chip_HepG2 = sum(df_state$ChiP_state == "HepG2"),
            n_Chip_K562 = sum(df_state$ChiP_state == "K562")
          )
        )
        
        saveRDS(object = df_state,
                file = file.path(
                  output_folder_protein_data,
                  paste(protein_name,
                        c_state,
                        "RDS",
                        sep = ".")
                ))
        
      } # end c_state loop
      
      saveRDS(object = df_temp,
              file = summery_file_name)
      
      # # lay <- rbind(
      # #   c(1, 2, 3, 4),
      # #   c(5, 6, 7, 8),
      # #   c(9, 9, 10, 11),
      # #   c(9, 9, 10, 11),
      # #   c(12, 12, 13, 14),
      # #   c(15, 15, 16, 17)
      # # )
      #
      # plot_list_general <- list(plot_ven, plot_bio_rep, plot_FPKM)
      #
      # # plot_bio_rep,
      # # plot_states,
      # # plot_islands_3corner,
      # # plot_islands_chip,
      # # plot_islands_hist,
      # #plot_CpG_quality_K562,
      # #plot_CpG_quality_HepG2,
      # #plot_methylation_read_depth,
      # #plot_methylation_state_read_depth
      # # plot_ven,
      # # plot_FPKM,
      #
      plot_list <- plot_list_chromatin_state_specific
      
      #append(plot_list_general,
      #       plot_list_chromatin_state_specific)
      
      final_plot <-
        gridExtra::grid.arrange(
          grobs = plot_list,
          shared_legend,
          top = protein_name, #paste(protein_name, Sys.Date()),
          ncol = 6
          # layout_matrix = lay
        )
      
      
      pic_file_name <-
        paste(protein_name,
              "overview",
              picuture_file_extension,
              sep = ".")
      
      ggplot2::ggsave(
        filename = pic_file_name,
        plot = final_plot,
        device = picuture_file_extension,
        path =  output_folder_protein_plots,
        width = 7680,
        height = 4320,
        units = "px"
      )
      
      print(paste("saved", pic_file_name))
    }
    print(sprintf(
      "%s Finished in %s on %s",
      protein_name,
      format(Sys.time() - start_time),
      sprintf("on %s.RDS",
              format(Sys.time(),
                     "%d-%b-%Y %H.%M"))
    ))
    return(df_temp)
    #processed_results_summery <- rbind(processed_results_summery,df_temp)
  }
#parallel::stopCluster(cl)
stopImplicitCluster()

saveRDS(object = processed_results_summery,
        file = file.path(
          this.dir,
          output_folder ,
          "meta",
          sprintf(
            "protein_sensetivity_summery %s.RDS",
            format(Sys.time(), "%d-%b-%Y %H.%M")
          )
        ))

print(sprintf(
  "%s Finished : %s %s",
  "Script",
  format(Sys.time() - start_script),
  sprintf(
    "protein_sensetivity_summery %s.RDS",
    format(Sys.time(),
           "%d-%b-%Y %H.%M")
  )
))
}