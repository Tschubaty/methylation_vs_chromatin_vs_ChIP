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
output_folder <- file.path(this.dir, "permutation_test_2")

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
      # biosample_comb <- biosample_combinations[[1]]
      print(biosample_comb)
      
      merge_df_folder <-  file.path(output_folder,
      protein,
      motif)
      
      # Create output directory if it doesn't exist
      if (!dir.exists(merge_df_folder)) {
        dir.create(merge_df_folder, recursive = TRUE)
      }
 
      merge_df_file_name <- file.path(
        merge_df_folder,
        paste(
          unlist(biosample_comb)[1],
          unlist(biosample_comb)[2],
          "best_comparison_all_data_points.rds",
          sep = "_"
        )
      )
      if (exists(merge_df_file_name)) {
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
      
      # Discard peaKS WITHOUT cpG 
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
      
      
      # --- Full join on coordinates only, keep original names ---
      merged_df <- dplyr::full_join(
        df1,
        df2,
        by = c("chr", "start_cg", "end_cg", "strand"),
        suffix = c(".x", ".y")
      ) %>%
        dplyr::mutate(
          # Define combined sample label
          sample = dplyr::case_when(
            !is.na(sample.x) & !is.na(sample.y) ~ paste0(biosample_comb[1], "_", biosample_comb[2]),
            !is.na(sample.x) ~ biosample_comb[1],
            !is.na(sample.y) ~ biosample_comb[2],
            TRUE ~ NA_character_
          ),
          # Define combined experiment label
          experiment = dplyr::case_when(
            !is.na(experiment.x) & !is.na(experiment.y) ~ paste0(experiment.x, "_", experiment.y),
            !is.na(experiment.x) ~ experiment.x,
            !is.na(experiment.y) ~ experiment.y,
            TRUE ~ NA_character_
          )
        ) %>%
        # Select and rename the columns back to the expected names
        dplyr::transmute(
          chr,
          start_cg,
          end_cg,
          strand,
          sample,
          experiment,
          !!paste0("Chromatin_State_", biosample_comb[1]) := coalesce(
            .data[[paste0("Chromatin_State_", biosample_comb[1], ".x")]],
            .data[[paste0("Chromatin_State_", biosample_comb[1], ".y")]]
          ),
          !!paste0("fRead_", biosample_comb[1]) := coalesce(
            .data[[paste0("fRead_", biosample_comb[1], ".x")]],
            .data[[paste0("fRead_", biosample_comb[1], ".y")]]
          ),
          !!paste0("Chromatin_State_", biosample_comb[2]) := coalesce(
            .data[[paste0("Chromatin_State_", biosample_comb[2], ".x")]],
            .data[[paste0("Chromatin_State_", biosample_comb[2], ".y")]]
          ),
          !!paste0("fRead_", biosample_comb[2]) := coalesce(
            .data[[paste0("fRead_", biosample_comb[2], ".x")]],
            .data[[paste0("fRead_", biosample_comb[2], ".y")]]
          )
        ) %>%
        dplyr::distinct()
      
      # Check the structure of the cleaned merged data
      #str(merged_df)
      
      saveRDS(object = merged_df, file =  merge_df_file_name)
      
    }
  }
  
}

############################## calculate statistic #############################

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

start <- which("NR3C1" == protein_folders) -1
stop <- start +3
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
        motif,
        paste(
          biosample_both,
          "best_comparison_all_data_points.rds",
          sep = "_"
        )
      )
      
      merged_df <- readRDS(rds_file)
      merged_df <- merged_df[complete.cases(merged_df), ]
      
      print(table(merged_df$sample))
      
      # === GRID PLOT with per-state counts ===
      chromatin_state_col1 <- paste0("Chromatin_State_", biosample_comb[1])
      chromatin_state_col2 <- paste0("Chromatin_State_", biosample_comb[2])
      fRead_col1 <- paste0("fRead_", biosample_comb[1])
      fRead_col2 <- paste0("fRead_", biosample_comb[2])
      
      grid_df <- merged_df %>%
        dplyr::filter(
          !is.na(.data[[chromatin_state_col1]]),
          !is.na(.data[[chromatin_state_col2]]),
          .data[[chromatin_state_col1]] == .data[[chromatin_state_col2]]
        ) %>%
        dplyr::mutate(
          state = factor(as.character(.data[[chromatin_state_col1]]), levels = as.character(1:18)),
          sample_label = dplyr::case_when(
            sample == biosample_comb[1] ~ biosample_comb[1],
            sample == biosample_comb[2] ~ biosample_comb[2],
            TRUE ~ "present in both samples"
          )
        )
      
      if (nrow(grid_df) == 0L) {
        grid_df <- tibble::tibble(
          !!fRead_col1 := numeric(),
          !!fRead_col2 := numeric(),
          state = factor(levels = as.character(1:18)),
          sample_label = character()
        )
      }
      
      count_1    <- sum(grid_df$sample_label == biosample_comb[1], na.rm = TRUE)
      count_2    <- sum(grid_df$sample_label == biosample_comb[2], na.rm = TRUE)
      count_both <- sum(grid_df$sample_label == "present in both samples", na.rm = TRUE)
      
      # per-state counts (force all 18 states AND all 3 categories to exist)
      labels_df <- grid_df %>%
        dplyr::count(state, sample_label, name = "n") %>%
        tidyr::complete(
          state = factor(as.character(1:18), levels = as.character(1:18)),
          sample_label = factor(
            c(biosample_comb[1], biosample_comb[2], "present in both samples"),
            levels = c(biosample_comb[1], biosample_comb[2], "present in both samples")
          ),
          fill = list(n = 0L)
        ) %>%
        tidyr::pivot_wider(names_from = sample_label, values_from = n) %>%
        dplyr::mutate(
          # after complete()+wider, these columns are guaranteed to exist; still replace any NAs defensively
          !!biosample_comb[1] := tidyr::replace_na(.data[[biosample_comb[1]]], 0L),
          !!biosample_comb[2] := tidyr::replace_na(.data[[biosample_comb[2]]], 0L),
          `present in both samples` = tidyr::replace_na(.data[["present in both samples"]], 0L),
          label = paste0(
            biosample_comb[1], "=", .data[[biosample_comb[1]]], "\n",
            biosample_comb[2], "=", .data[[biosample_comb[2]]], "\n",
            "both=", `present in both samples`
          ),
          # position for the text inside each facet (adjust if your fRead range differs)
          x = 0.02, y = 0.98
        )
      
      
      palette_now <- c(group.colors, "present in both samples" = "gray")
      
      grid_title <- paste0(
        motif, " | ", biosample_comb[1], " vs ", biosample_comb[2],
        "  —  counts: ",
        biosample_comb[1], "=", count_1, " | ",
        biosample_comb[2], "=", count_2, " | both=", count_both
      )
      
      grid_plot <- ggplot(grid_df, aes(
        x = .data[[fRead_col1]],
        y = .data[[fRead_col2]],
        color = sample_label
      )) +
        geom_point(size = 0.9, alpha = 0.85, na.rm = TRUE) +
        scale_color_manual(
          values = palette_now,
          name = "Category",
          breaks = c(biosample_comb[1], biosample_comb[2], "present in both samples")
        ) +
        facet_wrap(~ state, ncol = 6, drop = FALSE) +
        labs(
          title = grid_title,
          x = paste0("fRead_", biosample_comb[1]),
          y = paste0("fRead_", biosample_comb[2])
        ) +
        # fix axis so empty facets can still show annotations; adjust if your fRead range differs
        coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
        geom_text(
          data = labels_df,
          aes(x = x, y = y, label = label),
          inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 11, face = "bold"),
          legend.position = "right"
        )
      
      ggsave(
        filename = file.path(
          sim_output_dir,
          paste0(
            motif, "_", biosample_comb[1], "_vs_", biosample_comb[2],
            "_states_1to18_grid.svg"
          )
        ),
        plot = grid_plot,
        width = 12, height = 8, device = "svg", limitsize = FALSE
      )
      # === END GRID PLOT ===
      
      
      
      # --- NEW: loop over each chromatin state (1..18) ---
      for (state_i in 1:18) {
        data_set <- paste0("state_", state_i)
        
        # Equal-state filter + select the specific state
        current_df <- merged_df %>%
          dplyr::filter(
            !!sym(chromatin_state_col1) == state_i,
            !!sym(chromatin_state_col2) == state_i
          )
        
        sample_counts <- table(current_df$sample)
        
        # Handle sparse/no data
        if (nrow(current_df) == 0 || length(sample_counts) < 3) {
          stratified_test[[length(stratified_test) + 1]] <- data.frame(
            protein = protein,
            biosample1 = biosample_comb[1],
            biosample2 = biosample_comb[2],
            motif = motif,
            data_set = data_set,
            n_1 = ifelse(
              biosample_comb[1] %in% names(sample_counts),
              as.numeric(sample_counts[biosample_comb[1]]), 0
            ),
            n_2 = ifelse(
              biosample_comb[2] %in% names(sample_counts),
              as.numeric(sample_counts[biosample_comb[2]]), 0
            ),
            n_both = ifelse(
              biosample_both %in% names(sample_counts),
              as.numeric(sample_counts[biosample_both]), 0
            ),
            p_value_S_bigger = NA,
            p_value_S_smaller = NA,
            observed_S_statistic = NA
          )
          next
        }
        
        # Observed S
        observed_S_statistic <- calculate_S_statistic(current_df, biosample_comb)
        
        sim_output_dir <- file.path(motif_out_dir, biosample_both)
        dir.create(sim_output_dir, recursive = TRUE, showWarnings = FALSE)
        permuted_file <- file.path(sim_output_dir, paste0("permuted_stratified_S_df_state_", state_i, ".csv"))
        
        
        
        # Permutation (reuse cache if exists)
        if (file.exists(permuted_file)) { # 
          permuted_stratified_S_df <- read.csv(permuted_file)
          permuted_S_values <- permuted_stratified_S_df$S_statistic
        } else {
          permuted_S_values <- foreach(
            i = 1:n_permutations,
            .combine = 'c',
            .packages = c('dplyr', 'rlang')
          ) %dopar% {
            current_df %>%
              dplyr::group_by(
                !!sym(chromatin_state_col1),
                !!sym(chromatin_state_col2)
              ) %>%
              dplyr::mutate(sample = sample(sample)) %>%
              dplyr::ungroup() %>%
              calculate_S_statistic(biosample_comb)
          }
        }
        
        p_value_S_bigger  <- mean(permuted_S_values >= observed_S_statistic)
        p_value_S_smaller <- mean(permuted_S_values <= observed_S_statistic)
        
        cat(sprintf(
          "[%s] Observed S: %.4f | %d perms: %d >= (p=%s), %d <= (p=%s), two-sided p=%s\n",
          data_set,
          observed_S_statistic,
          length(permuted_S_values),
          sum(permuted_S_values >= observed_S_statistic),
          ifelse(p_value_S_bigger == 0, "< 1/length(permuted_S_values)", sprintf("%.4f", p_value_S_bigger)),
          sum(permuted_S_values <= observed_S_statistic),
          ifelse(p_value_S_smaller == 0, "< 1/length(permuted_S_values)", sprintf("%.4f", p_value_S_smaller)),
          ifelse(2 * min(p_value_S_bigger, p_value_S_smaller) == 0,
                 "< 1/length(permuted_S_values)",
                 sprintf("%.4f", 2 * min(p_value_S_bigger, p_value_S_smaller)))
        ))
        
        # Save permutation distribution
        permuted_stratified_S_df <- data.frame(S_statistic = permuted_S_values)
        write.csv(permuted_stratified_S_df, permuted_file, row.names = FALSE)
        
        # Histogram
        histogram_plot <- ggplot(permuted_stratified_S_df, aes(x = S_statistic)) +
          geom_histogram(bins = 100, fill = "lightblue", color = "black") +
          geom_vline(aes(xintercept = observed_S_statistic),
                     color = "red", linetype = "dashed", linewidth = 1) +
          labs(
            title = paste("Histogram of Permuted S Values for", data_set),
            x = "S_statistic", y = "Frequency"
          ) +
          theme_minimal()
        
        ggsave(
          file = file.path(sim_output_dir,
                           paste0("S_stratified_statistic_histogram_state_", state_i, ".svg")),
          plot = histogram_plot,
          device = "svg",
          limitsize = FALSE
        )
        
        # Scatter (A vs B)
        current_df <- current_df %>%
          dplyr::mutate(
            sample_label = dplyr::case_when(
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
              motif, data_set, "Scatterplot:",
              "CpG Methylation in", biosample_comb[1], "vs", biosample_comb[2]
            ),
            x = paste("CpG Methylation in", biosample_comb[1]),
            y = paste("CpG Methylation in", biosample_comb[2])
          ) +
          theme_minimal()
        
        ggsave(
          filename = file.path(
            sim_output_dir,
            paste0(
              motif, "_state_", state_i,
              "_Scatterplot_CpG_Methylation_",
              biosample_comb[1], "_vs_", biosample_comb[2],
              "_stratified.svg"
            )
          ),
          plot = scatter_plot,
          width = 8, height = 6, device = "svg"
        )
        
        # Append result row
        stratified_test[[length(stratified_test) + 1]] <- data.frame(
          protein = protein,
          biosample1 = biosample_comb[1],
          biosample2 = biosample_comb[2],
          motif = motif,
          data_set = data_set,   # <-- writes chromatin state
          n_1 = as.numeric(sample_counts[biosample_comb[1]]),
          n_2 = as.numeric(sample_counts[biosample_comb[2]]),
          n_both = as.numeric(sample_counts[biosample_both]),
          p_value_S_bigger = p_value_S_bigger,
          p_value_S_smaller = p_value_S_smaller,
          observed_S_statistic = observed_S_statistic
        )
      } # end state loop
    }
  }
}
  # Stop the cluster after all loops are done
  stopCluster(cl)
  
  # Combine results into a single data frame and save
  stratified_test <- do.call(rbind, stratified_test)
  rownames(stratified_test) <- NULL
  saveRDS(stratified_test,
          file = file.path(output_folder, "stratified_test_100k_26.10.2025.rds"))
  
  
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