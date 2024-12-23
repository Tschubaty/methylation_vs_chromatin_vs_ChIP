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
library(ggalluvial)
library(stringr)
library(gridExtra)

##################################### INPUT ########################################
# Define the input directory for chromosome files
input_WGBS_dir <- file.path(this.dir, "WGBS/byChr/concat_methylation")


input_ChIP_dir <- file.path(
  this.dir,
  "meme",
  "fimo_single_experiments_on_known_motif_and_peaks_add_methylation_v3"
)

# Set up parallel processing
n_cores <- detectCores() - 1
registerDoParallel(cores = n_cores)


# Define the directory where you want to save the plots
output_dir_plots <- file.path(this.dir, "plots")

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
  "1_TssA",
  # Active TSS
  "2_TssFlnk",
  # Flanking Active TSS
  "3_TssFlnkU",
  # Flanking TSS Upstream
  "4_TssFlnkD",
  # Flanking TSS Downstream
  "5_Tx",
  # Strong Transcription
  "6_TxWk",
  # Weak Transcription
  "7_EnhG1",
  # Genic Enhancer 1
  "8_EnhG2",
  # Genic Enhancer 2
  "9_EnhA1",
  # Active Enhancer 1
  "10_EnhA2",
  # Active Enhancer 2
  "11_EnhWk",
  # Weak Enhancer
  "12_ZNF_Rpts",
  # ZNF Genes & Repeats
  "13_Het",
  # Heterochromatin
  "14_TssBiv",
  # Bivalent TSS
  "15_EnhBiv",
  # Bivalent Enhancer
  "16_ReprPC",
  # Repressed Polycomb
  "17_ReprPCWk",
  # Weak Repressed Polycomb
  "18_Quies"      # Quiescent/Low Activity
)

# Define corresponding hex color codes for each chromatin state
chromatin_state_colors <- c(
  "#FF0000",
  # 1_TssA       -> Red
  "#FF4500",
  # 2_TssFlnk    -> Orange-Red
  "#FF9900",
  # 3_TssFlnkU   -> Orange
  "#FFCC00",
  # 4_TssFlnkD   -> Yellow-Orange
  "#00CC00",
  # 5_Tx         -> Green
  "#006400",
  # 6_TxWk       -> Dark Green
  "#FFD700",
  # 7_EnhG1      -> Gold
  "#FFD700",
  # 8_EnhG2      -> Gold (same as EnhG1)
  "#FFFF00",
  # 9_EnhA1      -> Yellow
  "#FFDD00",
  # 10_EnhA2     -> Yellow-Orange
  "#FFEA73",
  # 11_EnhWk     -> Light Yellow
  "#9370DB",
  # 12_ZNF_Rpts  -> Purple
  "#C0C0C0",
  # 13_Het       -> Light Gray
  "#FF4500",
  # 14_TssBiv    -> Orange-Red (similar to TssFlnk)
  "#FFDD00",
  # 15_EnhBiv    -> Yellow-Orange (similar to EnhA2)
  "#808080",
  # 16_ReprPC    -> Dark Gray
  "#A9A9A9",
  # 17_ReprPCWk  -> Light Gray
  "#000000"   # 18_Quies     -> Black
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


# List of biosamples and their respective ChIP directories
biosamples <- c("A549", "GM12878", "HepG2", "K562")

# Iterate through all proteins and their respective experiments
proteins <- list.files(path = input_ChIP_dir)

#######################################################################################
# # reads WGBS data
if (FALSE) {
  combined_methylation_data_compact <- readRDS(file.path(input_WGBS_dir, "combined_methylation_data_compact.rds"))
}
####################### make WGBS plots and save in list ###########################
if (FALSE) {
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
  
  saveRDS(object = plot_list,
          file = file.path(out_dir, "plot_list,rds"))
  rm(plot_list)
  rm(ecdf_plot)
  rm(final_grid_plot)
  rm(histogram_plot)
  rm(combined_plot)
  rm(state_data)
  gc()
}
####################### make grid plots and save plots  ###########################
if (FALSE) {
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
if (FALSE) {
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
    if (protein == "log")
      next
    
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
      if (!(biosample %in% biosamples))
        next
      
      # Read ChIP data
      start_read <- Sys.time()
      df_ChIP <- read_bed_file(file_path)
      end_read <- Sys.time()
      print(paste("Data loaded in:", round(
        difftime(end_read, start_read, units = "secs"), 2
      ), "seconds"))
      
      # Identify matching CpGs in combined_methylation_data_compact
      start_match <- Sys.time()
      match_idx <- which(
        combined_methylation_data_compact$chr %in% df_ChIP$chr &
          combined_methylation_data_compact$Start %in% df_ChIP$start_cg &
          combined_methylation_data_compact$End %in% df_ChIP$end_cg
      )
      end_match <- Sys.time()
      print(paste(
        "Index matching completed in:",
        round(difftime(end_match, start_match, units = "secs"), 2),
        "seconds"
      ))
      
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
        combined_methylation_data_compact[match_idx, column_name] <- vapply(match_idx, function(i) {
          current_value <- combined_methylation_data_compact[i, column_name]
          new_value <- paste0(protein, "_", experiment_id, "_", motif)
          
          if (is.na(current_value)) {
            new_value  # Insert new value if NA
          } else {
            paste(current_value, new_value, sep = ";")  # Append new value
          }
        }, FUN.VALUE = character(1)  # Ensure a character vector is returned)
        
      }
      end_update <- Sys.time()
      print(paste(
        "Protein information added in:",
        round(difftime(end_update, start_update, units = "secs"), 2),
        "seconds"
      ))
      
      # Total time for this iteration
      end_total <- Sys.time()
      print(paste("Total time for this file:", round(
        difftime(end_total, start_total, units = "secs"), 2
      ), "seconds"))
    }
  }
  
  
  # Final output
  print(head(combined_methylation_data_compact))
  saveRDS(
    object = combined_methylation_data_compact,
    file = file.path(
      out_dir,
      "combined_methylation_data_compact_withpeaks.rds"
    )
  )
  test <- head(combined_methylation_data_compact[!is.na(combined_methylation_data_compact$protein_hits_GM12878), ])
}
#################################################################################

combined_methylation_data_compact <- readRDS(file = file.path(out_dir, "combined_methylation_data_compact_withpeaks.rds"))
#lot_list <- readRDS(file = file.path(out_dir,"plot_list,rds"))

# if (FALSE) {
#   # List of biosamples
#   biosamples <- c("A549", "GM12878", "HepG2", "K562")
#   
#   # Iterate through each biosample
#   for (biosample in biosamples) {
#     print(paste("Processing biosample:", biosample))
#     
#     # Identify relevant columns
#     protein_hits_col <- paste0("protein_hits_", biosample)
#     methylation_col <- paste0("fRead_", biosample)
#     chromatin_state_col <- paste0("Chromatin_State_", biosample)
#     
#     # Subset only relevant columns and create 'is_occupied' column
#     biosample_methylation_data <- combined_methylation_data_compact %>%
#       select(!!sym(methylation_col),
#              !!sym(protein_hits_col),
#              !!sym(chromatin_state_col)) %>%
#       rename(
#         Methylation = !!sym(methylation_col),
#         ProteinHits = !!sym(protein_hits_col),
#         ChromatinState = !!sym(chromatin_state_col)
#       ) %>%
#       mutate(is_occupied = !is.na(ProteinHits))  # Logical for protein occupancy
#     
#     # Print summary of the processed data
#     print(paste(
#       "Total CpGs in",
#       biosample,
#       ":",
#       nrow(biosample_methylation_data)
#     ))
#     print(paste(
#       "Occupied CpGs:",
#       sum(biosample_methylation_data$is_occupied)
#     ))
#     print(paste(
#       "Unoccupied CpGs:",
#       sum(!biosample_methylation_data$is_occupied)
#     ))
#     
#     # Histogram plot
#     histogram_plot <- ggplot(biosample_methylation_data,
#                              aes(x = Methylation, fill = is_occupied)) +
#       geom_histogram(bins = 100,
#                      alpha = 0.7,
#                      position = "identity") +
#       scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "grey")) +
#       labs(
#         title = paste("Methylation Histogram -", biosample),
#         x = "Methylation Level",
#         y = "Frequency",
#         fill = "Occupancy"
#       ) +
#       theme_minimal()
#     
#     # Save histogram plot
#     histogram_file <- paste0("occupancy_histogram_methylation_", biosample, ".png")
#     ggsave(
#       filename = file.path(output_dir_plots, histogram_file),
#       plot = histogram_plot,
#       width = 10,
#       height = 6
#     )
#     print(paste("Histogram plot saved:", histogram_file))
#     
#     # Summarize data for the bar plot
#     plot_data <- biosample_methylation_data %>%
#       group_by(ChromatinState, is_occupied) %>%
#       summarize(Count = n(), .groups = "drop") %>%
#       mutate(is_occupied = ifelse(is_occupied, "Occupied", "Unoccupied"))
#     
#     # Bar plot
#     bar_plot <- ggplot(plot_data,
#                        aes(x = ChromatinState, y = Count, fill = is_occupied)) +
#       geom_bar(stat = "identity",
#                position = "stack",
#                alpha = 0.8) +
#       scale_fill_manual(values = c(
#         "Occupied" = "blue",
#         "Unoccupied" = "grey"
#       )) +
#       labs(
#         title = paste("CpG Occupancy in Chromatin States -", biosample),
#         x = "Chromatin State",
#         y = "Number of CpGs",
#         fill = "Occupancy"
#       ) +
#       theme_minimal() +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1),
#             # Rotate x-axis labels for readability
#             panel.grid.major.x = element_blank()  # Remove vertical grid lines)
#             
#             # Save bar plot
#             bar_plot_file <- paste0("occupancy_barplot_chromatin_state_", biosample, ".png")
#             ggsave(
#               filename = file.path(output_dir_plots, bar_plot_file),
#               plot = bar_plot,
#               width = 10,
#               height = 6
#             )
#             print(paste("Bar plot saved:", bar_plot_file))
#   }
# }
############# make plots of each protein  is_occupied #########################


# Define the output directory for plots
output_plots_subfolder <- file.path(output_dir_plots, "protein_presence_plots")
if (!dir.exists(output_plots_subfolder)) {
  dir.create(output_plots_subfolder, recursive = TRUE)
}
if (FALSE) {
  for (protein in proteins) {
    print(paste("Processing protein:", protein))
    
    # Skip log folder
    if (protein == "log")
      next
    
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
        select(
          !!sym(methylation_col),
          !!sym(protein_hits_col),
          !!sym(chromatin_state_col)
        ) %>%
        rename(
          Methylation = !!sym(methylation_col),
          ProteinHits = !!sym(protein_hits_col),
          ChromatinState = !!sym(chromatin_state_col)
        ) %>%
        mutate(
          ProteinHits = case_when(
            is.na(ProteinHits) ~ "not occupied",
            # If NA, not occupied
            str_detect(ProteinHits, fixed(peak_name)) ~ "peak_name",
            # If peak_name matches
            TRUE ~ "other protein"  # Otherwise, other proteins
          ),
          is_occupied = ProteinHits != "not occupied"  # Logical for occupancy
        )
      
      # Define a color palette for ProteinHits
      protein_hits_colors <- c(
        "not occupied" = "grey",
        "peak_name" = "red",
        "other protein" = "blue"
      )
      
      # Create separate histograms for each ProteinHits category
      histograms <- lapply(unique(biosample_methylation_data$ProteinHits), function(hit_type) {
        ggplot(
          biosample_methylation_data %>% filter(ProteinHits == hit_type),
          aes(x = Methylation, fill = ProteinHits)
        ) +
          geom_histogram(bins = 100,
                         alpha = 0.7,
                         position = "identity") +
          scale_fill_manual(values = protein_hits_colors) +
          labs(
            title = paste(
              "Methylation Histogram -",
              hit_type,
              "\nProtein:",
              protein,
              "Biosample:",
              biosample,
              "Motif:",
              motif,
              "Experiment:",
              experiment_id
            ),
            x = "Methylation Level",
            y = "Frequency"
          ) +
          theme_minimal()
      })
      
      # Combine histograms into a grid layout
      grid_histograms <- grid.arrange(
        grobs = histograms,
        nrow = length(histograms),
        top = paste(
          "Methylation Histograms for Protein:",
          protein,
          "Biosample:",
          biosample,
          "Motif:",
          motif,
          "Experiment:",
          experiment_id
        )
      )
      
      # Save the grid plot of histograms
      histogram_file <- file.path(
        output_plots_subfolder,
        paste0(
          "histogram_grid_",
          protein,
          "_",
          biosample,
          "_",
          experiment_id,
          "_",
          motif,
          ".png"
        )
      )
      ggsave(
        histogram_file,
        plot = grid_histograms,
        width = 10,
        height = 15
      )
      print(paste("Saved histogram plot:", histogram_file))
      
      # Generate a combined ECDF plot
      ecdf_plot_combined <- ggplot(biosample_methylation_data,
                                   aes(x = Methylation, color = ProteinHits)) +
        stat_ecdf(geom = "step", size = 1.2) +
        scale_color_manual(values = protein_hits_colors) +
        labs(
          title = paste(
            "Combined ECDF of Methylation\nProtein:",
            protein,
            "Biosample:",
            biosample,
            "Motif:",
            motif,
            "Experiment:",
            experiment_id
          ),
          x = "Methylation Level",
          y = "ECDF",
          color = "Protein Hits"
        ) +
        theme_minimal()
      
      # Save the ECDF plot
      ecdf_file <- file.path(
        output_plots_subfolder,
        paste0(
          "ecdf_combined_",
          protein,
          "_",
          biosample,
          "_",
          experiment_id,
          "_",
          motif,
          ".png"
        )
      )
      ggsave(ecdf_file,
             plot = ecdf_plot_combined,
             width = 10,
             height = 6)
      print(paste("Saved ECDF plot:", ecdf_file))
      
      # Save the processed data as an RDS file
      output_file <- file.path(
        output_plots_subfolder,
        paste0(
          "biosample_methylation_data_",
          protein,
          "_",
          biosample,
          "_",
          experiment_id,
          "_",
          motif,
          ".rds"
        )
      )
      saveRDS(biosample_methylation_data, output_file)
      print(
        paste(
          "Data saved for biosample:",
          biosample,
          "experiment:",
          experiment_id,
          "motif:",
          motif,
          "at",
          output_file
        )
      )
      
      # Log total time for processing this file
      end_total <- Sys.time()
      print(paste("Total time for this file:", round(
        difftime(end_total, start_total, units = "secs"), 2
      ), "seconds"))
    }
  }
}
############# make plots of each state & protein  is_occupied #########################

if (FALSE) {
  for (protein in proteins) {
    print(paste("Processing protein:", protein))
    
    # Skip log folder
    if (protein == "log")
      next
    
    # Get all files for the protein
    protein_out_dir <- file.path(output_plots_subfolder, protein)
    
    
    experiment_files <- list.files(
      pattern = paste0(protein , ".*.rds"),
      path = output_plots_subfolder,
      full.names = TRUE
    )
    
    for (file_path in experiment_files) {
      cat("\nStarting experiment file:", file_path, "\n")
      start_total <- Sys.time()
      
      # Extract experiment details
      file_name <- basename(file_path)
      experiment_id <- strsplit(file_name, "_")[[1]][6]
      biosample <- strsplit(file_name, "_")[[1]][5]
      motif <- strsplit(file_name, "_")[[1]][7]
      
      # Identify relevant columns
      protein_hits_col <- paste0("protein_hits_", biosample)
      methylation_col <- paste0("fRead_", biosample)
      chromatin_state_col <- paste0("Chromatin_State_", biosample)
      
      # Name for the current peak
      peak_name <- paste0(protein, "_", experiment_id, "_", motif)
      
      # Process biosample-specific data
      biosample_methylation_data <- readRDS(file = file_path)
      
      # Define a color palette for ProteinHits
      protein_hits_colors <- c(
        "not occupied" = "grey",
        "peak_name" = "red",
        "other protein" = "blue"
      )
      
      for (state in levels(biosample_methylation_data$ChromatinState)) {
        state_biosample_methylation_data <- biosample_methylation_data[biosample_methylation_data$ChromatinState == state, ]
        
        
        
        
        # Create separate histograms for each ProteinHits category
        histograms <- lapply(unique(state_biosample_methylation_data$ProteinHits), function(hit_type) {
          ggplot(
            state_biosample_methylation_data %>% filter(ProteinHits == hit_type),
            aes(x = Methylation, fill = ProteinHits)
          ) +
            geom_histogram(bins = 100,
                           alpha = 0.7,
                           position = "identity") +
            scale_fill_manual(values = protein_hits_colors) +
            labs(
              title = paste(
                "Methylation Histogram -",
                hit_type,
                "\nProtein:",
                protein,
                "Biosample:",
                biosample,
                "Motif:",
                motif,
                "Experiment:",
                experiment_id
              ),
              x = "Methylation Level",
              y = "Frequency"
            ) +
            theme_minimal()
        })
        
        # Combine histograms into a grid layout
        grid_histograms_state <- grid.arrange(
          grobs = histograms,
          nrow = length(histograms),
          top = paste(
            "state",
            state,
            protein,
            "Biosample:",
            biosample,
            "Motif:",
            motif,
            "Experiment:",
            experiment_id
          )
        )
        
        
        # Generate a combined ECDF plot
        ecdf_plot_state <- ggplot(state_biosample_methylation_data,
                                  aes(x = Methylation, color = ProteinHits)) +
          stat_ecdf(geom = "step", size = 1.2) +
          scale_color_manual(values = protein_hits_colors) +
          labs(
            title = paste(
              "Combined ECDF of Methylation\nProtein:",
              protein,
              "Biosample:",
              biosample,
              "Motif:",
              motif,
              "Experiment:",
              experiment_id
            ),
            x = "Methylation Level",
            y = "ECDF",
            color = "Protein Hits"
          ) +
          theme_minimal()
        
        # Log total time for processing this file
        end_total <- Sys.time()
        print(paste(
          "Total time for this file:",
          round(difftime(end_total, start_total, units = "secs"), 2),
          "seconds"
        ))
      }
    }
  }
}

############## chat gpt temp ########################################


for (protein in proteins[142]) {
  print(paste("Processing protein:", protein))
  
  # Skip log folder
  if (protein == "log")
    next
  
  # Get all experiment files for the protein
  protein_out_dir <- file.path(output_plots_subfolder, protein)
  
  # Check if the directory exists, and create it if it doesn't
  if (!dir.exists(protein_out_dir)) {
    dir.create(protein_out_dir, recursive = TRUE)
  }
  
  experiment_files <- list.files(
    pattern = paste0(protein , ".*.rds"),
    path = output_plots_subfolder,
    full.names = TRUE
  )
  
  # List to store all plots for grid layout
  plot_list <- list()
  
  for (file_path in experiment_files) {
    cat("\nStarting experiment file:", file_path, "\n")
    start_total <- Sys.time()
    
    # Extract experiment details
    file_name <- basename(file_path)
    experiment_id <- strsplit(file_name, "_")[[1]][6]
    biosample <- strsplit(file_name, "_")[[1]][5]
    motif <- strsplit(file_name, "_")[[1]][7]
    
    # Process biosample-specific data
    biosample_methylation_data <- readRDS(file = file_path)
    
    # Define a color palette for ProteinHits
    protein_hits_colors <- c(
      "not occupied" = "grey",
      "peak_name" = "red",
      "other protein" = "blue"
    )
    
    # Iterate through each chromatin state
    for (state in levels(biosample_methylation_data$ChromatinState)) {
      state_data <- biosample_methylation_data %>%
        filter(ChromatinState == state)
      
      # Generate Histogram Plot
      histograms <- lapply(unique(state_data$ProteinHits), function(hit_type) {
        ggplot(
          state_data %>% filter(ProteinHits == hit_type),
          aes(x = Methylation, fill = ProteinHits)
        ) +
          geom_histogram(bins = 100,
                         alpha = 0.7,
                         position = "identity") +
          scale_fill_manual(values = protein_hits_colors) +
          labs(
            title = paste("Histogram -", hit_type),
            x = "Methylation Level",
            y = "Frequency"
          ) +
          theme_minimal()
      })
      
      # Combine histograms into a grid layout
      grid_histograms_state <- grid.arrange(
        grobs = histograms,
        nrow = length(histograms),
        top = paste(
          "Histograms - State:",
          state,
          "motif:",
          motif,
          "Biosample:",
          biosample,
          "experiment_id",
          experiment_id
        )
      )
      
      # Generate ECDF Plot
      ecdf_plot_state <- ggplot(state_data, aes(x = Methylation, color = ProteinHits)) +
        stat_ecdf(geom = "step", size = 1.2) +
        scale_color_manual(values = protein_hits_colors) +
        labs(
          title = paste(
            "ECDF - State:",
            state,
            "motif:",
            motif,
            "Biosample:",
            biosample,
            "experiment_id",
            experiment_id
          ),
          x = "Methylation Level",
          y = "ECDF",
          color = "Protein Hits"
        ) +
        theme_minimal()
      
      # Combine both plots (histogram + ECDF) for the state
      combined_state_plot <- grid.arrange(grid_histograms_state, ecdf_plot_state, nrow = 2,
                                          top = paste("State:", state, "Protein:", protein, "Biosample:", biosample,"motif",motif))
                                          
                                          # Add to the plot list for the experiment
                                          plot_list[[paste(state, experiment_id, motif, sep = "_")]] <- combined_state_plot
    }
    
    # Log total time for processing this experiment
    end_total <- Sys.time()
    print(paste("Total time for this file:", round(
      difftime(end_total, start_total, units = "secs"), 2
    ), "seconds"))
  }
  
  # Combine all plots into a grid layout for the protein
  combined_protein_plot <- grid.arrange(
    grobs = plot_list,
    ncol = 18,
    top = paste("Overview of Protein:", protein)
  )
  
  # Save the combined plot
  ggsave(
    filename = file.path(protein_out_dir, paste0("overview_", protein, ".svg")),
    plot = combined_protein_plot,
    width = 180,
    height = length(experiment_files) * 10,
    limitsize = FALSE
  )
  
  print(paste("Saved overview plot for protein:", protein))
}



############# Vhat gpt end

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

if (FALSE) {
  # Iterate through all proteins and their respective experiments
  proteins <- list.files(path = input_ChIP_dir)
  for (protein in proteins) {
    print(protein)
    # Skip log folder
    if (protein == "log")
      next
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
    mutate(Chromatin_State = factor(
      Chromatin_State,
      levels = c(
        names(chromatin_state_colors_short),
        "No motif in peak",
        "No CG in motif",
        "No State assignment for CG"
      )
    ))
  
  # Save the summary to a CSV file
  output_csv <- file.path(out_dir, "chromatin_state_tests.csv")
  write_csv(final_summary, output_csv)
}
###############################################################################
if (FALSE) {
  # Loop through each protein folder
  for (protein_file in protein_files[1:10]) {
    # protein_file in protein_files[1:10]
    protein_name <- gsub(
      pattern = ".fimo_methylation_peaks.bed",
      replacement = "",
      x = basename(protein_file)
    )  # Extract protein name
    df <- read.delim(
      header = FALSE,
      file = file.path(input_folder, protein_file),
      col.names = colnames_df
    )
    head(df)
    
    ChIP_files <- list.files(path = file.path(ChIP.seq_folder, protein_name),
                             pattern = "*.bed")
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

