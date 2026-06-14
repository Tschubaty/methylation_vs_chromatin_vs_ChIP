#################################################################
##  Chip-seq vs Methylation vs Chromatin Analysis (V3)
##
##  Description: Permutation test for ALL experiment pairs within 
##               each protein+motif combination (not filtered).
##               Results annotated with antibody, experiment_id, file_id
##
##  Input: Methylation and protein peak data (HG38)
##         Located in Encode3/meme/fimo_single_experiments_on_known_motif_and_peaks_add_methylation_v3/
##
##  Output: Results with full metadata saved in:
##          Encode3/permutation_test_v3/
##
##  Version: 10.01.2025
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

#################################### Libraries ###################################
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rlang)
library(doParallel)
library(foreach)
library(tibble)

#################################### Constants ###################################
# Start time for the script
start_script <- Sys.time()

# Define input and output directories
input_folder <- file.path(
  this.dir,
  "meme",
  "fimo_single_experiments_on_known_motif_and_peaks_add_methylation_v3"
)
output_folder <- file.path(this.dir, "permutation_test_v3")

# Create output directory if it doesn't exist
if (!dir.exists(output_folder)) {
  dir.create(output_folder, recursive = TRUE)
}

group.colors <- c(
  HepG2    = "#F8766D",
  K562     = "#00BFC4",
  GM12878  = "#A3A500",
  A549     = "#E76BF3"
)

chromatin_state_colors_short <- c(
  "1" = "#FF0000", "2" = "#FF4500", "3" = "#FF9900", "4" = "#FFCC00",
  "5" = "#00CC00", "6" = "#006400", "7" = "#FFD700", "8" = "#FFD700",
  "9" = "#FFFF00", "10" = "#FFDD00", "11" = "#FFEA73", "12" = "#9370DB",
  "13" = "#C0C0C0", "14" = "#FF4500", "15" = "#FFDD00", "16" = "#808080",
  "17" = "#A9A9A9", "18" = "#000000"
)

# Helper function to extract column names from BED files
extract_colnames <- function(file_path) {
  first_row <- readLines(file_path, n = 1)
  cleaned_row <- gsub("#", "", first_row)
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
# Build comprehensive summary of all experiments
summary_df_path <- file.path(output_folder, "summary_df_all_experiments.rds")
if (!file.exists(summary_df_path)) {
  summary_list <- list()
  
  for (protein in protein_folders) {
    print(paste("Processing protein:", protein))
    
    files <- list.files(path = file.path(input_folder, protein))
    
    for (bed_file in files) {
      print(paste("Processing file:", bed_file))
      
      # Extract metadata from file name
      # Format: biosample_protein_experiment_id_motif_ChIP_with_methylation.bed
      parts <- strsplit(bed_file, "_")[[1]]
      biosample <- parts[1]
      protein_from_file <- parts[2]
      experiment_id <- parts[3]
      motif <- sub("\\.bed$", "", parts[4])
      
      # Read the BED file
      df <- read_bed_file(file_path = file.path(input_folder, protein, bed_file))
      
      # Create a data frame with only unique peaks
      df_unique_peaks <- df %>% distinct(chr, start_peak, end_peak, .keep_all = TRUE)
      
      # Calculate categories
      n_peak_with_CpG <- sum(df_unique_peaks$start_motif != -1 &
                               df_unique_peaks$start_cg != 0)
      n_peak_no_CpG <- sum(df_unique_peaks$start_motif != -1 &
                             df_unique_peaks$start_cg == 0)
      n_peak_no_motif_no_CpG <- sum(df_unique_peaks$start_motif == -1 &
                                      df_unique_peaks$start_cg == 0)
      n_peaks <- nrow(df_unique_peaks)
      
      # Append results to summary list (keeping all experiments)
      summary_list[[length(summary_list) + 1]] <- data.frame(
        protein = protein,
        biosample = biosample,
        experiment_id = experiment_id,
        motif = motif,
        n_peaks = n_peaks,
        n_peak_with_CpG = n_peak_with_CpG,
        n_peak_no_CpG = n_peak_no_CpG,
        n_peak_no_motif_no_CpG = n_peak_no_motif_no_CpG,
        n_CpG = sum(df$start_cg != -1),
        bed_file = bed_file,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Combine all summaries into a single data frame
  summary_df <- do.call(rbind, summary_list)
  rownames(summary_df) <- NULL
  
  # Save the summary data frame
  saveRDS(object = summary_df, file = summary_df_path)
  print(paste("Summary data saved at:", summary_df_path))
} else {
  # Load existing summary data
  summary_df <- readRDS(summary_df_path)
  print(paste("Loaded existing summary data from:", summary_df_path))
}

end_script <- Sys.time()
print(paste("Summary generation completed in", round(
  difftime(end_script, start_script, units = "mins"), 2
), "minutes"))

############################## Permutation Test #############################

# Start time for the script
start_script <- Sys.time()

# Fast S statistic using only the two exclusive groups.
#
# Input:
#   diff_excl:
#     Numeric methylation-difference vector for CpGs that belong to
#     biosample1-only or biosample2-only groups.
#
#   label_excl:
#     Labels for diff_excl. Must contain only biosample1 and biosample2 labels.
#
#   label1:
#     biosample1-only label, for example "HepG2"
#
#   label2:
#     biosample2-only label, for example "K562"
#
# Statistic:
#   S = mean(diff_excl[label1]) - mean(diff_excl[label2])
calculate_S_statistic_exclusive <- function(diff_excl, label_excl, label1, label2) {
  
  mean_1 <- mean(diff_excl[label_excl == label1], na.rm = TRUE)
  mean_2 <- mean(diff_excl[label_excl == label2], na.rm = TRUE)
  
  if (is.nan(mean_1) || is.nan(mean_2)) {
    return(NA_real_)
  }
  
  mean_1 - mean_2
}


# CURRENT (no progress):
# num_cores <- detectCores() - 1
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)

# NEW (with progress support):
library(doSNOW)
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores, type = "SOCK")
registerDoSNOW(cl)

stratified_test <- list()
n_permutations <- 100000 #
min_n_per_group <- 10

#protein_folders <- protein_folders[4:5]

# Loop over each protein
for (protein in protein_folders) {
  #protein <- protein_folders[1]
  print(paste("=== Processing protein:", protein, "==="))
  
  protein_data <- summary_df[summary_df$protein == protein, ]
  
  # Loop over unique motifs
  unique_motifs <- unique(protein_data$motif)
  protein_output_dir <- file.path(output_folder, protein)
  
  for (motif in unique_motifs) {
    # motif <- unique_motifs[1]
    print(paste("Motif:", motif))
    motif_data <- protein_data[protein_data$motif == motif, ]
    
    # Get ALL experiment pairs (not filtered by quality)
    # Generate all combinations of rows (experiments) within this motif
    if (nrow(motif_data) < 2) {
      print(paste("  Skip: only 1 experiment for", protein, motif))
      next
    }
    
    # Get all pairwise experiment combinations.
    # Important: only compare experiments from DIFFERENT biosamples.
    # We do not compare K562-vs-K562, HepG2-vs-HepG2, etc.,
    # because methylation values are sample-specific and would be identical
    # on both sides of the comparison.
    motif_rows <- 1:nrow(motif_data)
    
    experiment_pairs <- combn(motif_rows, 2, simplify = FALSE)
    
    experiment_pairs <- experiment_pairs[
      sapply(experiment_pairs, function(pair) {
        biosample_a <- motif_data$biosample[pair[1]]
        biosample_b <- motif_data$biosample[pair[2]]
        biosample_a != biosample_b
      })
    ]
    
    # If no valid cross-biosample pairs remain, skip this motif.
    if (length(experiment_pairs) == 0) {
      print(paste("  Skip: no cross-biosample experiment pairs for", protein, motif))
      next
    }
    
    motif_out_dir <- file.path(protein_output_dir, motif)
    
    # Loop through each experiment pair
    for (pair_idx in seq_along(experiment_pairs)) {
      # pair_idx <- seq_along(experiment_pairs)[1]
      pair <- experiment_pairs[[pair_idx]]
      idx1 <- pair[1]
      idx2 <- pair[2]
      
      # Get experiment metadata
      exp1_row <- motif_data[idx1, ]
      exp2_row <- motif_data[idx2, ]
      
      biosample1 <- exp1_row$biosample
      biosample2 <- exp2_row$biosample
      experiment_id1 <- exp1_row$experiment_id
      experiment_id2 <- exp2_row$experiment_id
      bed_file1 <- exp1_row$bed_file
      bed_file2 <- exp2_row$bed_file
      
      biosample_both <- paste(biosample1, biosample2, sep = "_")
      exp_both <- paste(experiment_id1, experiment_id2, sep = "_")
      
      print(paste("  Pair:", pair_idx, "/", length(experiment_pairs),
                  "-", biosample1, "[", experiment_id1, "] vs",
                  biosample2, "[", experiment_id2, "]"))
      
      # Create output directory
      pair_out_dir <- file.path(motif_out_dir, exp_both)
      dir.create(pair_out_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Read the BED files
      df1 <- read_bed_file(file_path = file.path(input_folder, protein, bed_file1))
      df2 <- read_bed_file(file_path = file.path(input_folder, protein, bed_file2))
      
      # Discard peaks without CpG
      df1 <- df1[df1$start_cg != 0, ]
      df2 <- df2[df2$start_cg != 0, ]
      
      # Convert all "." values to NA
      df1 <- df1 %>% mutate_all(~ replace(., . == ".", NA))
      df2 <- df2 %>% mutate_all(~ replace(., . == ".", NA))
      
      # Convert relevant columns to integers
      integer_columns <- c(
        "start_cg", "end_cg", "start_motif", "end_motif", "start_peak", "end_peak",
        "mRead_A549", "nRead_A549", "mRead_GM12878", "nRead_GM12878",
        "mRead_HepG2", "nRead_HepG2", "mRead_K562", "nRead_K562"
      )
      
      df1[integer_columns] <- lapply(df1[integer_columns], function(x)
        as.integer(as.character(x)))
      df2[integer_columns] <- lapply(df2[integer_columns], function(x)
        as.integer(as.character(x)))
      
      # Remove duplicates
      df1 <- df1 %>% distinct(chr, start_cg, end_cg, strand, .keep_all = TRUE)
      df2 <- df2 %>% distinct(chr, start_cg, end_cg, strand, .keep_all = TRUE)
      
      # Recalculate fRead columns
      for (biosample_col in c(biosample1, biosample2)) {
        for (df_name in c("df1", "df2")) {
          df_temp <- get(df_name)
          df_temp <- df_temp %>%
            mutate(across(
              starts_with("fRead_"),
              .fns = function(x) {
                base_col <- sub("fRead_", "", cur_column())
                mRead_col <- paste0("mRead_", base_col)
                nRead_col <- paste0("nRead_", base_col)
                ifelse(df_temp[[nRead_col]] == 0, NA, df_temp[[mRead_col]] / df_temp[[nRead_col]])
              },
              .names = "{col}"
            ))
          assign(df_name, df_temp, envir = parent.frame())
        }
      }
      
      # Full join on coordinates
      biosample_comb <- c(biosample1, biosample2)
      merged_df <- dplyr::full_join(
        df1, df2,
        by = c("chr", "start_cg", "end_cg", "strand"),
        suffix = c(".x", ".y")
      ) %>%
        dplyr::mutate(
          sample = dplyr::case_when(
            !is.na(sample.x) & !is.na(sample.y) ~ paste0(biosample1, "_", biosample2),
            !is.na(sample.x) ~ biosample1,
            !is.na(sample.y) ~ biosample2,
            TRUE ~ NA_character_
          ),
          experiment = dplyr::case_when(
            !is.na(experiment.x) & !is.na(experiment.y) ~ paste0(experiment.x, "_", experiment.y),
            !is.na(experiment.x) ~ experiment.x,
            !is.na(experiment.y) ~ experiment.y,
            TRUE ~ NA_character_
          )
        ) %>%
        dplyr::transmute(
          chr,
          start_cg, end_cg, strand,
          sample, experiment,
          !!paste0("Chromatin_State_", biosample1) := coalesce(
            .data[[paste0("Chromatin_State_", biosample1, ".x")]],
            .data[[paste0("Chromatin_State_", biosample1, ".y")]]
          ),
          !!paste0("fRead_", biosample1) := coalesce(
            .data[[paste0("fRead_", biosample1, ".x")]],
            .data[[paste0("fRead_", biosample1, ".y")]]
          ),
          !!paste0("Chromatin_State_", biosample2) := coalesce(
            .data[[paste0("Chromatin_State_", biosample2, ".x")]],
            .data[[paste0("Chromatin_State_", biosample2, ".y")]]
          ),
          !!paste0("fRead_", biosample2) := coalesce(
            .data[[paste0("fRead_", biosample2, ".x")]],
            .data[[paste0("fRead_", biosample2, ".y")]]
          )
        ) %>%
        dplyr::distinct()
      
      # Remove incomplete cases
      merged_df <- merged_df[complete.cases(merged_df), ]
      
      if (nrow(merged_df) == 0) {
        print(paste("    No complete data for this pair"))
        next
      }
      
      # ===== LOOP OVER CHROMATIN STATES =====
      chromatin_state_col1 <- paste0("Chromatin_State_", biosample1)
      chromatin_state_col2 <- paste0("Chromatin_State_", biosample2)
      
      for (state_i in 1:18) {
        # state_i <- 1
        # Filter CpGs where BOTH biosamples are in the same chromatin state.
        # Even if this gives zero rows, we will still write an output row.
        current_df <- merged_df %>%
          dplyr::filter(
            !!sym(chromatin_state_col1) == state_i,
            !!sym(chromatin_state_col2) == state_i
          )
        
        # Count how many CpGs belong to biosample1-only, biosample2-only, and both.
        # If current_df is empty, all counts will correctly become 0.
        sample_counts <- table(current_df$sample)
        
        n_sample1 <- as.numeric(ifelse(
          biosample1 %in% names(sample_counts),
          sample_counts[biosample1],
          0
        ))
        
        n_sample2 <- as.numeric(ifelse(
          biosample2 %in% names(sample_counts),
          sample_counts[biosample2],
          0
        ))
        
        n_both <- as.numeric(ifelse(
          biosample_both %in% names(sample_counts),
          sample_counts[biosample_both],
          0
        ))
        
        # Default values.
        # These stay NA unless the state has enough data to run the permutation test.
        observed_S_statistic <- NA_real_
        p_value_S_bigger <- NA_real_
        p_value_S_smaller <- NA_real_
        p_value_two_sided <- NA_real_
        
        # Add a status column so you can later see WHY a row has NA statistics.
        test_status <- NA_character_
        
        # Case 1: no CpGs in this chromatin state.
        if (nrow(current_df) == 0) {
          
          test_status <- "no_CpGs_in_state"
          
        } else if (n_sample1 < min_n_per_group ||
                   n_sample2 < min_n_per_group) {
          
          # Case 2: CpGs exist, but at least one of the two exclusive groups is either:
          #   1. completely missing, n = 0
          #   2. present but below the minimum CpG threshold
          #
          # Required groups for the current S statistic:
          #   - biosample1-only CpGs
          #   - biosample2-only CpGs
          #
          # The shared/both group is saved as annotation but is not required,
          # because the current S statistic ignores it.
          
          low_groups <- c()
          
          if (n_sample1 == 0) {
            low_groups <- c(
              low_groups,
              paste0("missing_", biosample1, "_only_n0")
            )
          } else if (n_sample1 < min_n_per_group) {
            low_groups <- c(
              low_groups,
              paste0("too_few_", biosample1, "_only_n", n_sample1)
            )
          }
          
          if (n_sample2 == 0) {
            low_groups <- c(
              low_groups,
              paste0("missing_", biosample2, "_only_n0")
            )
          } else if (n_sample2 < min_n_per_group) {
            low_groups <- c(
              low_groups,
              paste0("too_few_", biosample2, "_only_n", n_sample2)
            )
          }
          
          test_status <- paste(low_groups, collapse = ";")
          
        } else {
          
          # Case 3: enough data exists, so calculate observed statistic
          # and run/read the permutation test.
          test_status <- "tested"
          
          # Precompute methylation difference once for this chromatin state.
          # One value per CpG:
          #   fRead_biosample1 - fRead_biosample2
          diff_fRead <- current_df[[paste0("fRead_", biosample1)]] -
            current_df[[paste0("fRead_", biosample2)]]
          
          # Store ChIP-binding group labels once.
          sample_labels <- current_df$sample
          
          # Keep only the two exclusive groups used by the S statistic.
          # The shared/both group is ignored by the current statistic.
          exclusive_idx <- sample_labels %in% c(biosample1, biosample2)
          
          # Remove shared/both CpGs before the permutation.
          diff_excl <- diff_fRead[exclusive_idx]
          label_excl <- sample_labels[exclusive_idx]
          
          # Remove rows with NA methylation differences.
          # This keeps the permutation vector clean and avoids repeated na.rm work.
          valid_excl <- !is.na(diff_excl)
          
          diff_excl <- diff_excl[valid_excl]
          label_excl <- label_excl[valid_excl]
          
          # Count the two exclusive groups.
          n1_perm <- sum(label_excl == biosample1)
          n2_perm <- sum(label_excl == biosample2)
          
          # Observed S statistic.
          observed_S_statistic <- calculate_S_statistic_exclusive(
            diff_excl = diff_excl,
            label_excl = label_excl,
            label1 = biosample1,
            label2 = biosample2
          )
          
          # Include the number of permutations in the file name,
          # so old 100-permutation files are not reused for 100k runs.
          permuted_file <- file.path(
            pair_out_dir,
            paste0("permuted_stratified_S_df_state_", state_i, "_nperm_", n_permutations, ".csv")
          )
          
          if (file.exists(permuted_file)) {
            
            # Read existing permutation distribution if already calculated
            permuted_stratified_S_df <- read.csv(permuted_file)
            permuted_S_values <- permuted_stratified_S_df$S_statistic
            
          } else {
            
            cat(sprintf("      State %d: Running %d permutations...\n", state_i, n_permutations))
            p_start_time <- Sys.time()
            
            pb <- txtProgressBar(max = n_permutations, style = 3)
            opts <- list(progress = function(n) setTxtProgressBar(pb, n))
            
            sum_total <- sum(diff_excl)
            n_excl <- length(diff_excl)
            
            permuted_S_values <- foreach(
              i = 1:n_permutations,
              .combine = 'c',
              .options.snow = opts
            ) %dopar% {
              
              # Randomly choose n1_perm CpGs for the permuted biosample1-only group.
              idx_group1 <- sample.int(n_excl, n1_perm)
              
              # Sum for permuted biosample1-only group.
              sum_group1 <- sum(diff_excl[idx_group1])
              
              # Sum for permuted biosample2-only group.
              # This avoids creating diff_excl[-idx_group1].
              sum_group2 <- sum_total - sum_group1
              
              # Permuted S statistic.
              (sum_group1 / n1_perm) - (sum_group2 / n2_perm)
            }
            
            close(pb)
            elapsed <- difftime(Sys.time(), p_start_time, units = "secs")
            cat(sprintf("Completed in %.1f seconds\n", elapsed))
            
            # Save permutation distribution for later reuse
            permuted_stratified_S_df <- data.frame(S_statistic = permuted_S_values)
            write.csv(permuted_stratified_S_df, permuted_file, row.names = FALSE)
          }
          
          # Calculate empirical p-values
          p_value_S_bigger <- mean(permuted_S_values >= observed_S_statistic)
          p_value_S_smaller <- mean(permuted_S_values <= observed_S_statistic)
          
          # Two-sided empirical p-value.
          # min(1, ...) prevents values above 1.
          p_value_two_sided <- min(1, 2 * min(p_value_S_bigger, p_value_S_smaller))
        }
        
        # Always append one row per chromatin state.
        # If the test could not be run, counts are still saved and statistics stay NA.
        stratified_test[[length(stratified_test) + 1]] <- data.frame(
          protein = protein,
          motif = motif,
          chromatin_state = state_i,
          
          # Experiment 1 metadata
          biosample1 = biosample1,
          experiment_id1 = experiment_id1,
          bed_file1 = bed_file1,
          
          # Experiment 2 metadata
          biosample2 = biosample2,
          experiment_id2 = experiment_id2,
          bed_file2 = bed_file2,
          
          # Sample counts
          n_sample1 = n_sample1,
          n_sample2 = n_sample2,
          n_both = n_both,
          
          # Test results
          observed_S_statistic = observed_S_statistic,
          p_value_S_bigger = p_value_S_bigger,
          p_value_S_smaller = p_value_S_smaller,
          p_value_two_sided = p_value_two_sided,
          
          # Reason why test was / was not run
          test_status = test_status,
          
          stringsAsFactors = FALSE
        )
        
      } # end state loop
    } # end experiment pair loop
  } # end motif loop
} # end protein loop

# Stop the cluster
stopCluster(cl)

# Combine results into a single data frame and save
stratified_test_df <- do.call(rbind, stratified_test)
rownames(stratified_test_df) <- NULL

# Create dynamic output file names based on the number of permutations
output_prefix <- paste0("stratified_test_all_pairs_", n_permutations, "perm")

# Save as RDS
saveRDS(
  stratified_test_df,
  file = file.path(output_folder, paste0(output_prefix, ".rds"))
)

# Save as Excel
if (require("openxlsx")) {
  openxlsx::write.xlsx(
    stratified_test_df,
    file = file.path(output_folder, paste0(output_prefix, ".xlsx")),
    overwrite = TRUE
  )
}

print(paste("Permutation test complete. Results saved to:", output_folder))
print(paste("Total test rows:", nrow(stratified_test_df)))

end_script <- Sys.time()
print(paste("Total script time:", round(
  difftime(end_script, start_script, units = "mins"), 2
), "minutes"))
