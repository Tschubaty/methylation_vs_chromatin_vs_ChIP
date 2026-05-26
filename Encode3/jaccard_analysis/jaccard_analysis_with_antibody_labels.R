#################################################################
##  Calculate Jaccard Index for ChIP-seq Peaks by Antibody
##
##  Compares peak overlap between experiments per protein-motif
##  Uses interval overlap (any 1bp overlap = match)
##
##  Input: Encode3/meme/fimo_single_experiments_on_known_motif_and_peaks_add_methylation_v3
##  Output: Encode3/jaccard_analysis/
##
##  Author: Daniel Batyrev
##  Date: 2026-05-14
#################################################################

# Clear R working environment
rm(list = ls())

cluster <- FALSE
if (cluster) {
  this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/Encode3/"
  picture_file_extension <- "pdf"
} else {
  this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
  picture_file_extension <- "png"
}

setwd(this.dir)

#################################### Libs ########################################
library(stringr)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(GenomicRanges)
library(IRanges)

################################## Constants #####################################
start_script <- Sys.time()

# Set input and output directories
input_folder <- file.path(
  dirname(this.dir),
  "meme",
  "fimo_single_experiments_on_known_motif_and_peaks_add_methylation_v3"
)
output_folder <- file.path(this.dir)

############################ Main-loop cache settings ############################

# Set TRUE only if you want to redo the expensive protein/motif loop from scratch.
force_run_main_jaccard_loop <- FALSE

# If this summary already exists, the whole protein/motif loop can be skipped.
main_jaccard_summary_csv <- file.path(
  output_folder,
  "ALL_PROTEINS_MOTIFS_jaccard_summary.csv"
)

# RDS version is better for later reuse inside R.
main_jaccard_summary_rds <- file.path(
  output_folder,
  "ALL_PROTEINS_MOTIFS_jaccard_summary.rds"
)

main_jaccard_cache_available <- file.exists(main_jaccard_summary_rds) ||
  file.exists(main_jaccard_summary_csv)

# Antibody metadata file
antibody_mapping_file <- file.path(
  dirname(this.dir),
  "Antibody_Info",
  "ENCODE_antibody_complete_mapping.csv"
)

# Suggested colors for the biosamples
group.colors <- c(
  HepG2   = "#F8766D",     # A warm reddish color
  K562    = "#00BFC4",     # A cool cyan color
  GM12878 = "#A3A500",     # A yellow-green color
  A549    = "#E76BF3"      # A vibrant purple color
)

################################ Functions ######################################

# Function to extract and clean column names from the first row of a text file
extract_colnames <- function(file_path) {
  first_row <- readLines(file_path, n = 1)
  cleaned_row <- gsub("#", "", first_row)          # Remove leading #
  colnames <- unlist(strsplit(cleaned_row, "\t"))
  return(colnames)
}

# Function to read .bed files with dynamically extracted column names
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

# Calculate Jaccard Index with interval overlap (any 1bp overlap = match)
calculate_jaccard_peaks_overlap <- function(peaks_A, peaks_B) {
  # Convert to GRanges (genomic ranges) for interval overlap
  gr_A <- GRanges(
    seqnames = peaks_A$chr,
    ranges = IRanges(start = peaks_A$start_peak, end = peaks_A$end_peak)
  )
  
  gr_B <- GRanges(
    seqnames = peaks_B$chr,
    ranges = IRanges(start = peaks_B$start_peak, end = peaks_B$end_peak)
  )
  
  # Remove duplicate peaks within each set
  gr_A <- unique(gr_A)
  gr_B <- unique(gr_B)
  
  # Find overlapping peaks (any 1 bp overlap = match)
  overlaps <- findOverlaps(gr_A, gr_B, minoverlap = 1L)
  
  # Count unique peaks from A that overlap with B
  n_intersection <- length(unique(queryHits(overlaps)))
  
  # Union = total unique peaks across both sets
  # (peaks in A + peaks in B - overlapping peaks)
  n_union <- length(gr_A) + length(gr_B) - n_intersection
  
  # Jaccard Index
  if (n_union == 0) {
    jaccard <- 0
  } else {
    jaccard <- n_intersection / n_union
  }
  
  return(jaccard)
}


# Make multiline labels for heatmap rows/columns
make_heatmap_label <- function(biosample, Experiment_Accession, antibody_accession,experiment_id) {
  antibody_accession <- ifelse(
    is.na(antibody_accession) | antibody_accession == "" | antibody_accession == "N/A",
    "no_antibody_ID",
    antibody_accession
  )
  
  paste0(
    "Cell line: ", biosample, "\n",
    "ID: ", Experiment_Accession, "\n",
    "Antibody ID: ", antibody_accession
  )
}

################################ Antibody Mapping ###############################

if (!file.exists(antibody_mapping_file)) {
  stop(paste("Antibody mapping file not found:", antibody_mapping_file))
}

antibody_mapping_raw <- readr::read_csv(
  antibody_mapping_file,
  show_col_types = FALSE
)

antibody_mapping <- antibody_mapping_raw %>%
  mutate(
    File_Accession = as.character(File_Accession),
    Antibody_Accession = as.character(Antibody_Accession),
    Antibody_ID_clean = stringr::str_remove_all(as.character(Antibody_ID), "^/antibodies/|/$"),
    Antibody_Label = dplyr::case_when(
      !is.na(Antibody_Accession) & Antibody_Accession != "" & Antibody_Accession != "N/A" ~ Antibody_Accession,
      !is.na(Antibody_ID_clean) & Antibody_ID_clean != "" & Antibody_ID_clean != "N/A" ~ Antibody_ID_clean,
      TRUE ~ "no_antibody_ID"
    )
  ) %>%
  group_by(File_Accession) %>%
  summarise(
    Experiment_Accession = paste(unique(na.omit(Experiment_Accession)), collapse = ";"),
    Antibody_Title = paste(unique(na.omit(Antibody_Title)), collapse = ";"),
    Antibody_ID = paste(unique(na.omit(Antibody_ID_clean)), collapse = ";"),
    Antibody_Accession = paste(unique(na.omit(Antibody_Accession)), collapse = ";"),
    Antibody_Label = paste(unique(na.omit(Antibody_Label)), collapse = ";"),
    Antibody_Target = paste(unique(na.omit(Antibody_Target)), collapse = ";"),
    Antibody_Lot = paste(unique(na.omit(Antibody_Lot)), collapse = ";"),
    Antibody_Catalog = paste(unique(na.omit(Antibody_Catalog)), collapse = ";"),
    Antibody_Lab = paste(unique(na.omit(Antibody_Lab)), collapse = ";"),
    Experiment_Lab = paste(unique(na.omit(Experiment_Lab)), collapse = ";"),
    .groups = "drop"
  )

################################ Main Loop ######################################

# Get the list of all protein folders
protein_folders <- list.dirs(
  path = input_folder,
  recursive = FALSE,
  full.names = FALSE
)

# Exclude sbatch_scripts and log folders
protein_folders <- protein_folders[!grepl("sbatch_scripts|log", protein_folders)]

# Sort for consistency
protein_folders <- sort(protein_folders)

############################ Main-loop cache check ###############################

if (!force_run_main_jaccard_loop && main_jaccard_cache_available) {
  
  print("Global Jaccard summary already exists.")
  print("Skipping entire protein/motif loop.")
  print("BED files will NOT be read.")
  print("Heatmaps will NOT be recreated.")
  
  if (file.exists(main_jaccard_summary_rds)) {
    all_results_df <- readRDS(main_jaccard_summary_rds)
  } else {
    all_results_df <- readr::read_csv(
      main_jaccard_summary_csv,
      show_col_types = FALSE
    ) %>%
      as.data.frame()
    
    saveRDS(
      all_results_df,
      file = main_jaccard_summary_rds
    )
  }
  
  # Recreate this variable in case later code expects it.
  all_jaccard_results <- split(all_results_df, seq_len(nrow(all_results_df)))
  
} else {
  
  print("No global Jaccard summary cache found, or recalculation forced.")
  print("Running full protein/motif loop.")
  
  # Initialize list to store all results
  all_jaccard_results <- list()
  
  # Loop over each protein (ChIP target)
  for (protein in protein_folders) {
  # protein <- protein_folders[1]  # For debugging: run single protein
  
  print(paste("Processing protein:", protein))
  
  # Get all files for the current protein
  files <- list.files(path = file.path(input_folder, protein))
  
  if (length(files) == 0) {
    print(paste("  No files found for", protein))
    next
  }
  
  # Extract unique motifs and biosamples from file names
  file_metadata <- data.frame()
  
  for (bed_file in files) {
    if (!grepl("\\.bed$", bed_file)) next  # Skip non-BED files
    
    # Split the file name by "_"
    parts <- strsplit(bed_file, "_")[[1]]
    
    # Extract the components based on position
    biosample <- parts[1]            # Example: GM12878
    experiment_id <- parts[3]        # Example: ENCFF320KXO
    motif <- sub("\\.bed$", "", parts[4])  # Extract motif
    
    file_metadata <- rbind(
      file_metadata,
      data.frame(
        bed_file = bed_file,
        biosample = biosample,
        experiment_id = experiment_id,
        motif = motif,
        experiment_label = paste(biosample, experiment_id, sep = "_"),
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (nrow(file_metadata) == 0) {
    print(paste("  Could not parse any files for", protein))
    next
  }
  
  file_metadata <- file_metadata %>%
    left_join(
      antibody_mapping,
      by = c("experiment_id" = "File_Accession")
    ) %>%
    mutate(
      Antibody_Label = ifelse(
        is.na(Antibody_Label) | Antibody_Label == "" | Antibody_Label == "N/A",
        "no_antibody_ID",
        Antibody_Label
      ),
      
      # Keep this one-line label for joins, CSV tables, and matrix rownames
      experiment_label = paste(biosample, experiment_id, Antibody_Label, sep = "_"),
      
      # Use this only for the heatmap display
      heatmap_label = make_heatmap_label(
        biosample = biosample,
        Experiment_Accession = Experiment_Accession,
        antibody_accession = Antibody_Label,
        experiment_id = experiment_id
      )
    )
  
  # Get unique motifs for this protein
  unique_motifs <- unique(file_metadata$motif)
  
  # Loop over each motif
  for (motif in unique_motifs) {
    # motif <- unique_motifs[1]  # For debugging: run single motif
    
    print(paste("  Processing motif:", motif))
    
    # Filter metadata for this protein-motif combination
    motif_files <- file_metadata[file_metadata$motif == motif, ]
    
    if (nrow(motif_files) < 2) {
      print(paste("    Only", nrow(motif_files), "experiment(s) found. Skipping."))
      next
    }
    
    print(paste("    Found", nrow(motif_files), "experiments"))
    
    # Create output directory for this protein-motif
    protein_motif_dir <- file.path(output_folder, protein, motif)
    if (!dir.exists(protein_motif_dir)) {
      dir.create(protein_motif_dir, recursive = TRUE)
    }
    
    # Save experiment-to-antibody mapping for this protein-motif
    write.csv(
      motif_files,
      file = file.path(protein_motif_dir, "experiment_antibody_metadata.csv"),
      row.names = FALSE
    )
    
    # Initialize Jaccard matrix
    ######################## Load cached Jaccard or calculate ########################
    
    jaccard_matrix_file <- file.path(protein_motif_dir, "jaccard_matrix.rds")
    jaccard_table_file <- file.path(protein_motif_dir, "jaccard_table.csv")
    jaccard_summary_file <- file.path(protein_motif_dir, "jaccard_summary.csv")
    
    cached_jaccard_available <- file.exists(jaccard_matrix_file) &&
      file.exists(jaccard_table_file) &&
      file.exists(jaccard_summary_file)
    
    if (cached_jaccard_available) {
      
      print("    Cached Jaccard files found. Loading saved results; BED files will NOT be read.")
      
      jaccard_matrix <- readRDS(jaccard_matrix_file)
      
      jaccard_long <- readr::read_csv(
        jaccard_table_file,
        show_col_types = FALSE,
        col_types = readr::cols(
          .default = readr::col_character(),
          jaccard_index = readr::col_double(),
          same_biosample = readr::col_logical(),
          same_antibody = readr::col_logical()
        )
      )
      
      summary_stats <- readr::read_csv(
        jaccard_summary_file,
        show_col_types = FALSE
      ) %>%
        as.data.frame()
      
    } else {
      
      print("    Cached Jaccard files not found, or recalculation forced. Reading BED files now.")
      
      # Initialize Jaccard matrix
      n_exp <- nrow(motif_files)
      jaccard_matrix <- matrix(NA, nrow = n_exp, ncol = n_exp)
      rownames(jaccard_matrix) <- motif_files$experiment_label
      colnames(jaccard_matrix) <- motif_files$experiment_label
      
      # Read all peak data for this motif
      peak_data <- list()
      for (i in 1:nrow(motif_files)) {
        bed_file_path <- file.path(input_folder, protein, motif_files$bed_file[i])
        print(paste("    Reading:", motif_files$bed_file[i]))
        
        peak_data[[i]] <- read_bed_file(bed_file_path)
      }
      
      # Calculate pairwise Jaccard indices
      print("    Calculating Jaccard indices...")
      for (i in 1:n_exp) {
        for (j in i:n_exp) {
          
          peaks_i <- peak_data[[i]][, c("chr", "start_peak", "end_peak")]
          peaks_j <- peak_data[[j]][, c("chr", "start_peak", "end_peak")]
          
          peaks_i <- peaks_i[complete.cases(peaks_i), ]
          peaks_j <- peaks_j[complete.cases(peaks_j), ]
          
          jaccard_value <- calculate_jaccard_peaks_overlap(peaks_i, peaks_j)
          
          jaccard_matrix[i, j] <- jaccard_value
          jaccard_matrix[j, i] <- jaccard_value
        }
      }
      
      # Save Jaccard matrix as RDS
      saveRDS(
        jaccard_matrix,
        file = jaccard_matrix_file
      )
      
      # Convert matrix to long format data frame
      jaccard_df <- as.data.frame(jaccard_matrix)
      jaccard_df$experiment_1 <- rownames(jaccard_matrix)
      
      experiment_lookup <- motif_files %>%
        select(
          experiment_label,
          biosample,
          experiment_id,
          Antibody_Label,
          Antibody_Accession,
          Antibody_Title,
          Antibody_Lot,
          Antibody_Catalog,
          Antibody_Lab,
          Experiment_Lab
        )
      
      jaccard_long <- jaccard_df %>%
        tidyr::pivot_longer(
          cols = -experiment_1,
          names_to = "experiment_2",
          values_to = "jaccard_index"
        ) %>%
        dplyr::filter(experiment_1 < experiment_2) %>%
        dplyr::left_join(
          experiment_lookup,
          by = c("experiment_1" = "experiment_label")
        ) %>%
        dplyr::rename(
          biosample_1 = biosample,
          experiment_id_1 = experiment_id,
          antibody_label_1 = Antibody_Label,
          antibody_accession_1 = Antibody_Accession,
          antibody_title_1 = Antibody_Title,
          antibody_lot_1 = Antibody_Lot,
          antibody_catalog_1 = Antibody_Catalog,
          antibody_lab_1 = Antibody_Lab,
          experiment_lab_1 = Experiment_Lab
        ) %>%
        dplyr::left_join(
          experiment_lookup,
          by = c("experiment_2" = "experiment_label")
        ) %>%
        dplyr::rename(
          biosample_2 = biosample,
          experiment_id_2 = experiment_id,
          antibody_label_2 = Antibody_Label,
          antibody_accession_2 = Antibody_Accession,
          antibody_title_2 = Antibody_Title,
          antibody_lot_2 = Antibody_Lot,
          antibody_catalog_2 = Antibody_Catalog,
          antibody_lab_2 = Antibody_Lab,
          experiment_lab_2 = Experiment_Lab
        ) %>%
        dplyr::mutate(
          same_biosample = biosample_1 == biosample_2,
          same_antibody = antibody_label_1 == antibody_label_2
        ) %>%
        dplyr::arrange(dplyr::desc(jaccard_index))
      
      write.csv(
        jaccard_long,
        file = jaccard_table_file,
        row.names = FALSE
      )
      
      # Calculate summary statistics
      jaccard_values <- jaccard_matrix[upper.tri(jaccard_matrix)]
      jaccard_values <- jaccard_values[!is.na(jaccard_values)]
      
      same_antibody_values <- jaccard_long$jaccard_index[jaccard_long$same_antibody]
      different_antibody_values <- jaccard_long$jaccard_index[!jaccard_long$same_antibody]
      
      summary_stats <- data.frame(
        protein = protein,
        motif = motif,
        n_experiments = n_exp,
        n_antibodies = length(unique(motif_files$Antibody_Label)),
        n_pairwise_comparisons = length(jaccard_values),
        mean_jaccard = mean(jaccard_values, na.rm = TRUE),
        median_jaccard = median(jaccard_values, na.rm = TRUE),
        min_jaccard = min(jaccard_values, na.rm = TRUE),
        max_jaccard = max(jaccard_values, na.rm = TRUE),
        sd_jaccard = sd(jaccard_values, na.rm = TRUE),
        mean_jaccard_same_antibody = ifelse(
          length(same_antibody_values) > 0,
          mean(same_antibody_values, na.rm = TRUE),
          NA
        ),
        mean_jaccard_different_antibody = ifelse(
          length(different_antibody_values) > 0,
          mean(different_antibody_values, na.rm = TRUE),
          NA
        ),
        stringsAsFactors = FALSE
      )
      
      write.csv(
        summary_stats,
        file = jaccard_summary_file,
        row.names = FALSE
      )
    }
    
    # Generate heatmap
    png(
      file = file.path(
        protein_motif_dir,
        paste0("jaccard_heatmap_", protein, "_", motif, ".png")
      ),
      width = 3200,
      height = 2600,
      res = 180
    )
    
    pheatmap(
      jaccard_matrix,
      main = paste("Jaccard Index -", protein, "(", motif, ")"),
      display_numbers = TRUE,
      number_format = "%.3f",
      breaks = seq(0, 1, length.out = 101),
      color = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
      
      labels_row = motif_files$heatmap_label,
      labels_col = motif_files$heatmap_label,
      angle_col = 0,
      
      cellwidth = 220,
      cellheight = 150,
      
      fontsize = 14,
      fontsize_row = 13,
      fontsize_col = 13,
      fontsize_number = 18,
      
      border_color = "grey40"
    )
    
    dev.off()
    
    print(paste("    Saved outputs to:", protein_motif_dir))
    print(paste("    Mean Jaccard Index:", round(summary_stats$mean_jaccard, 4)))
    print(paste("    Range:", round(summary_stats$min_jaccard, 4), "-", 
                round(summary_stats$max_jaccard, 4)))
    
    # Add to all results
    all_jaccard_results[[length(all_jaccard_results) + 1]] <- summary_stats
  }
}

################################ Summary Report #################################

if (length(all_jaccard_results) > 0) {
  # Combine all results
  all_results_df <- do.call(rbind, all_jaccard_results)
  rownames(all_results_df) <- NULL
  
  # Save comprehensive summary
  write.csv(
    all_results_df,
    file = file.path(output_folder, "ALL_PROTEINS_MOTIFS_jaccard_summary.csv"),
    row.names = FALSE
  )
  
  # Create visualization of all results
  png(
    file = file.path(
      protein_motif_dir,
      paste0("jaccard_heatmap_", protein, "_", motif, ".png")
    ),
    width = 2200,
    height = 1800,
    res = 150
  )
  
  ggplot(all_results_df, aes(x = reorder(paste(protein, motif, sep = "_"), mean_jaccard),
                             y = mean_jaccard,
                             fill = protein)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = mean_jaccard - sd_jaccard,
                      ymax = mean_jaccard + sd_jaccard),
                  width = 0.3) +
    coord_flip() +
    labs(
      title = "Mean Jaccard Index by Protein-Motif",
      x = "Protein - Motif",
      y = "Mean Jaccard Index",
      fill = "Protein"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  dev.off()
  
  print("\n========================================")
  print("SUMMARY RESULTS")
  print("========================================")
  print(all_results_df)
  print("\n")
  print(paste("Total protein-motif combinations analyzed:", nrow(all_results_df)))
  print(paste("Overall mean Jaccard Index:", round(mean(all_results_df$mean_jaccard), 4)))
  print(paste("Overall range:", 
              round(min(all_results_df$min_jaccard), 4), "-",
              round(max(all_results_df$max_jaccard), 4)))
  
  # Save RDS cache so future runs can skip the entire main loop
  saveRDS(
    all_results_df,
    file = main_jaccard_summary_rds
  )
}
  
} # closes: if cached summary exists, skip loop; else run full main loop

end_script <- Sys.time()
print(paste("Total runtime:", round(difftime(end_script, start_script, units = "mins"), 2), "minutes"))
print(paste("Results saved to:", output_folder))

################################ Jaccard Grouped Analysis ########################
##################################################################################


jaccard_grouped_output_dir <- file.path(output_folder, "Jaccard_grouped_analysis")

if (!dir.exists(jaccard_grouped_output_dir)) {
  dir.create(jaccard_grouped_output_dir, recursive = TRUE)
}

print("Starting grouped Jaccard analysis...")

# Find all jaccard_table.csv files created in the protein/motif folders
jaccard_table_files <- list.files(
  path = output_folder,
  pattern = "^jaccard_table\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(jaccard_table_files) == 0) {
  warning("No jaccard_table.csv files found. Grouped Jaccard analysis skipped.")
}
  
  # Read all pairwise Jaccard tables and add protein/motif from folder names
all_jaccard_pairs <- lapply(jaccard_table_files, function(f) {
  
  motif_name <- basename(dirname(f))
  protein_name <- basename(dirname(dirname(f)))
  
  readr::read_csv(
    f,
    show_col_types = FALSE,
    col_types = readr::cols(
      .default = readr::col_character(),
      jaccard_index = readr::col_double(),
      same_biosample = readr::col_logical(),
      same_antibody = readr::col_logical()
    )
  ) %>%
    dplyr::mutate(
      protein = protein_name,
      motif = motif_name
    )
}) %>%
  dplyr::bind_rows()
  
  # Clean and create grouping variables
  all_jaccard_pairs <- all_jaccard_pairs %>%
    dplyr::mutate(
      antibody_label_1 = dplyr::if_else(
        is.na(antibody_label_1) | antibody_label_1 == "" | antibody_label_1 == "N/A",
        "no_antibody_ID",
        antibody_label_1
      ),
      antibody_label_2 = dplyr::if_else(
        is.na(antibody_label_2) | antibody_label_2 == "" | antibody_label_2 == "N/A",
        "no_antibody_ID",
        antibody_label_2
      ),
      
      same_biosample = biosample_1 == biosample_2,
      same_antibody = antibody_label_1 == antibody_label_2,
      
      antibody_known_1 = antibody_label_1 != "no_antibody_ID",
      antibody_known_2 = antibody_label_2 != "no_antibody_ID",
      both_antibodies_known = antibody_known_1 & antibody_known_2,
      
      same_antibody_strict = same_antibody & both_antibodies_known,
      different_antibody_strict = !same_antibody & both_antibodies_known,
      
      biosample_group = dplyr::if_else(
        same_biosample,
        "Same cell line",
        "Different cell line"
      ),
      
      antibody_group = dplyr::case_when(
        same_antibody_strict ~ "Same antibody",
        different_antibody_strict ~ "Different antibody",
        TRUE ~ "Unknown antibody"
      ),
      
      combined_group = dplyr::case_when(
        same_biosample & same_antibody_strict ~ "Same cell line + same antibody",
        same_biosample & different_antibody_strict ~ "Same cell line + different antibody",
        !same_biosample & same_antibody_strict ~ "Different cell line + same antibody",
        !same_biosample & different_antibody_strict ~ "Different cell line + different antibody",
        TRUE ~ "Unknown antibody"
      ),
      
      biosample_pair = paste(pmin(biosample_1, biosample_2), pmax(biosample_1, biosample_2), sep = " vs "),
      antibody_pair = paste(pmin(antibody_label_1, antibody_label_2), pmax(antibody_label_1, antibody_label_2), sep = " vs "),
      
      protein_motif = paste(protein, motif, sep = " - ")
    )
  
  # Save complete combined pairwise table
  write.csv(
    all_jaccard_pairs,
    file = file.path(jaccard_grouped_output_dir, "ALL_pairwise_jaccard_values_with_groups.csv"),
    row.names = FALSE
  )
  
  ############################ Summary tables ###################################
  
  jaccard_summary_by_biosample_group <- all_jaccard_pairs %>%
    dplyr::group_by(biosample_group) %>%
    dplyr::summarise(
      n_pairs = dplyr::n(),
      mean_jaccard = mean(jaccard_index, na.rm = TRUE),
      median_jaccard = median(jaccard_index, na.rm = TRUE),
      sd_jaccard = sd(jaccard_index, na.rm = TRUE),
      min_jaccard = min(jaccard_index, na.rm = TRUE),
      max_jaccard = max(jaccard_index, na.rm = TRUE),
      .groups = "drop"
    )
  
  jaccard_summary_by_antibody_group <- all_jaccard_pairs %>%
    dplyr::group_by(antibody_group) %>%
    dplyr::summarise(
      n_pairs = dplyr::n(),
      mean_jaccard = mean(jaccard_index, na.rm = TRUE),
      median_jaccard = median(jaccard_index, na.rm = TRUE),
      sd_jaccard = sd(jaccard_index, na.rm = TRUE),
      min_jaccard = min(jaccard_index, na.rm = TRUE),
      max_jaccard = max(jaccard_index, na.rm = TRUE),
      .groups = "drop"
    )
  
  jaccard_summary_by_combined_group <- all_jaccard_pairs %>%
    dplyr::group_by(combined_group) %>%
    dplyr::summarise(
      n_pairs = dplyr::n(),
      mean_jaccard = mean(jaccard_index, na.rm = TRUE),
      median_jaccard = median(jaccard_index, na.rm = TRUE),
      sd_jaccard = sd(jaccard_index, na.rm = TRUE),
      min_jaccard = min(jaccard_index, na.rm = TRUE),
      max_jaccard = max(jaccard_index, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(mean_jaccard))
  
  jaccard_summary_by_biosample_pair <- all_jaccard_pairs %>%
    dplyr::group_by(biosample_pair) %>%
    dplyr::summarise(
      n_pairs = dplyr::n(),
      mean_jaccard = mean(jaccard_index, na.rm = TRUE),
      median_jaccard = median(jaccard_index, na.rm = TRUE),
      sd_jaccard = sd(jaccard_index, na.rm = TRUE),
      min_jaccard = min(jaccard_index, na.rm = TRUE),
      max_jaccard = max(jaccard_index, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(mean_jaccard))
  
  jaccard_summary_by_antibody_pair <- all_jaccard_pairs %>%
    dplyr::filter(both_antibodies_known) %>%
    dplyr::group_by(antibody_pair) %>%
    dplyr::summarise(
      n_pairs = dplyr::n(),
      mean_jaccard = mean(jaccard_index, na.rm = TRUE),
      median_jaccard = median(jaccard_index, na.rm = TRUE),
      sd_jaccard = sd(jaccard_index, na.rm = TRUE),
      min_jaccard = min(jaccard_index, na.rm = TRUE),
      max_jaccard = max(jaccard_index, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(mean_jaccard))
  
  jaccard_summary_by_protein_motif_groups <- all_jaccard_pairs %>%
    dplyr::group_by(protein, motif, biosample_group, antibody_group, combined_group) %>%
    dplyr::summarise(
      n_pairs = dplyr::n(),
      mean_jaccard = mean(jaccard_index, na.rm = TRUE),
      median_jaccard = median(jaccard_index, na.rm = TRUE),
      sd_jaccard = sd(jaccard_index, na.rm = TRUE),
      min_jaccard = min(jaccard_index, na.rm = TRUE),
      max_jaccard = max(jaccard_index, na.rm = TRUE),
      .groups = "drop"
    )
  
  write.csv(
    jaccard_summary_by_biosample_group,
    file = file.path(jaccard_grouped_output_dir, "summary_by_biosample_group.csv"),
    row.names = FALSE
  )
  
  write.csv(
    jaccard_summary_by_antibody_group,
    file = file.path(jaccard_grouped_output_dir, "summary_by_antibody_group.csv"),
    row.names = FALSE
  )
  
  write.csv(
    jaccard_summary_by_combined_group,
    file = file.path(jaccard_grouped_output_dir, "summary_by_combined_group.csv"),
    row.names = FALSE
  )
  
  write.csv(
    jaccard_summary_by_biosample_pair,
    file = file.path(jaccard_grouped_output_dir, "summary_by_biosample_pair.csv"),
    row.names = FALSE
  )
  
  write.csv(
    jaccard_summary_by_antibody_pair,
    file = file.path(jaccard_grouped_output_dir, "summary_by_antibody_pair.csv"),
    row.names = FALSE
  )
  
  write.csv(
    jaccard_summary_by_protein_motif_groups,
    file = file.path(jaccard_grouped_output_dir, "summary_by_protein_motif_groups.csv"),
    row.names = FALSE
  )
  
  ############################ Basic statistical tests ###########################
  
  statistical_tests <- list()
  
  # Same cell line vs different cell line
  biosample_test_df <- all_jaccard_pairs %>%
    dplyr::filter(!is.na(jaccard_index))
  
  if (length(unique(biosample_test_df$biosample_group)) == 2) {
    statistical_tests$biosample_group_wilcox <- wilcox.test(
      jaccard_index ~ biosample_group,
      data = biosample_test_df
    )
  }
  
  # Same antibody vs different antibody, only where antibody is known
  antibody_test_df <- all_jaccard_pairs %>%
    dplyr::filter(
      antibody_group %in% c("Same antibody", "Different antibody"),
      !is.na(jaccard_index)
    )
  
  if (length(unique(antibody_test_df$antibody_group)) == 2) {
    statistical_tests$antibody_group_wilcox <- wilcox.test(
      jaccard_index ~ antibody_group,
      data = antibody_test_df
    )
  }
  
  sink(file.path(jaccard_grouped_output_dir, "jaccard_group_statistical_tests.txt"))
  cat("Jaccard grouped statistical tests\n")
  cat("================================\n\n")
  
  cat("Same cell line vs different cell line\n")
  cat("------------------------------------\n")
  if (!is.null(statistical_tests$biosample_group_wilcox)) {
    print(statistical_tests$biosample_group_wilcox)
  } else {
    cat("Test not possible: fewer than two groups available.\n")
  }
  
  cat("\n\nSame antibody vs different antibody\n")
  cat("----------------------------------\n")
  if (!is.null(statistical_tests$antibody_group_wilcox)) {
    print(statistical_tests$antibody_group_wilcox)
  } else {
    cat("Test not possible: fewer than two groups available.\n")
  }
  
  sink()
  
  ############################ Plots ############################################
  
  # 1. Overall histogram
  png(
    file = file.path(jaccard_grouped_output_dir, "histogram_all_jaccard_values.png"),
    width = 2200,
    height = 1600,
    res = 180
  )
  
  print(
    ggplot(all_jaccard_pairs, aes(x = jaccard_index)) +
      geom_histogram(bins = 50, color = "black") +
      labs(
        title = "Distribution of all pairwise Jaccard values",
        x = "Jaccard index",
        y = "Number of pairwise comparisons"
      ) +
      theme_minimal(base_size = 18)
  )
  
  dev.off()
  
  # 2. Histogram by same/different cell line
  png(
    file = file.path(jaccard_grouped_output_dir, "histogram_jaccard_by_biosample_group.png"),
    width = 2400,
    height = 1700,
    res = 180
  )
  
  print(
    ggplot(all_jaccard_pairs, aes(x = jaccard_index, fill = biosample_group)) +
      geom_histogram(bins = 50, alpha = 0.65, position = "identity", color = "black") +
      labs(
        title = "Jaccard distribution: same vs different cell line",
        x = "Jaccard index",
        y = "Number of pairwise comparisons",
        fill = "Comparison group"
      ) +
      theme_minimal(base_size = 18)
  )
  
  dev.off()
  
  # 3. Histogram by same/different antibody
  png(
    file = file.path(jaccard_grouped_output_dir, "histogram_jaccard_by_antibody_group.png"),
    width = 2400,
    height = 1700,
    res = 180
  )
  
  print(
    ggplot(
      all_jaccard_pairs %>% dplyr::filter(antibody_group != "Unknown antibody"),
      aes(x = jaccard_index, fill = antibody_group)
    ) +
      geom_histogram(bins = 50, alpha = 0.65, position = "identity", color = "black") +
      labs(
        title = "Jaccard distribution: same vs different antibody",
        x = "Jaccard index",
        y = "Number of pairwise comparisons",
        fill = "Comparison group"
      ) +
      theme_minimal(base_size = 18)
  )
  
  dev.off()
  
  # 4. Boxplot: same vs different cell line
  png(
    file = file.path(jaccard_grouped_output_dir, "boxplot_jaccard_same_vs_different_cell_line.png"),
    width = 2200,
    height = 1600,
    res = 180
  )
  
  print(
    ggplot(all_jaccard_pairs, aes(x = biosample_group, y = jaccard_index, fill = biosample_group)) +
      geom_boxplot(outlier.alpha = 0.35) +
      geom_jitter(width = 0.15, alpha = 0.35, size = 1.8) +
      labs(
        title = "Jaccard overlap: same vs different cell line",
        x = "",
        y = "Jaccard index",
        fill = "Comparison group"
      ) +
      theme_minimal(base_size = 18) +
      theme(legend.position = "none")
  )
  
  dev.off()
  
  # 5. Boxplot: same vs different antibody
  png(
    file = file.path(jaccard_grouped_output_dir, "boxplot_jaccard_same_vs_different_antibody.png"),
    width = 2200,
    height = 1600,
    res = 180
  )
  
  print(
    ggplot(
      all_jaccard_pairs %>% dplyr::filter(antibody_group != "Unknown antibody"),
      aes(x = antibody_group, y = jaccard_index, fill = antibody_group)
    ) +
      geom_boxplot(outlier.alpha = 0.35) +
      geom_jitter(width = 0.15, alpha = 0.35, size = 1.8) +
      labs(
        title = "Jaccard overlap: same vs different antibody",
        x = "",
        y = "Jaccard index",
        fill = "Comparison group"
      ) +
      theme_minimal(base_size = 18) +
      theme(legend.position = "none")
  )
  
  dev.off()
  
  # 6. Combined group boxplot
  png(
    file = file.path(jaccard_grouped_output_dir, "boxplot_jaccard_combined_cell_line_and_antibody.png"),
    width = 2800,
    height = 1800,
    res = 180
  )
  
  print(
    ggplot(
      all_jaccard_pairs %>% dplyr::filter(combined_group != "Unknown antibody"),
      aes(x = reorder(combined_group, jaccard_index, FUN = median), y = jaccard_index, fill = combined_group)
    ) +
      geom_boxplot(outlier.alpha = 0.35) +
      geom_jitter(width = 0.15, alpha = 0.30, size = 1.6) +
      coord_flip() +
      labs(
        title = "Jaccard overlap grouped by cell line and antibody",
        x = "",
        y = "Jaccard index",
        fill = "Comparison group"
      ) +
      theme_minimal(base_size = 18) +
      theme(legend.position = "none")
  )
  
  dev.off()
  
  # 7. Mean Jaccard by biosample pair
  png(
    file = file.path(jaccard_grouped_output_dir, "barplot_mean_jaccard_by_biosample_pair.png"),
    width = 2400,
    height = 1800,
    res = 180
  )
  
  png(
    file = file.path(jaccard_grouped_output_dir, "jitter_boxplot_jaccard_by_biosample_pair_colored_by_antibody.png"),
    width = 2600,
    height = 1900,
    res = 180
  )
  
  print(
    ggplot(
      all_jaccard_pairs %>%
        dplyr::filter(antibody_group != "Unknown antibody"),
      aes(
        x = reorder(biosample_pair, jaccard_index, FUN = median),
        y = jaccard_index
      )
    ) +
      geom_boxplot(
        outlier.shape = NA,
        fill = "grey85",
        color = "black"
      ) +
      geom_jitter(
        aes(color = antibody_group),
        width = 0.18,
        height = 0,
        alpha = 0.65,
        size = 2.3
      ) +
      stat_summary(
        fun = mean,
        geom = "point",
        shape = 23,
        size = 4,
        fill = "white",
        color = "black"
      ) +
      scale_color_manual(
        values = c(
          "Same antibody" = "#0072B2",
          "Different antibody" = "#D55E00"
        )
      ) +
      coord_flip() +
      labs(
        title = "Jaccard value distribution by cell-line pair",
        subtitle = "Each point is one pairwise experiment comparison; diamond = mean",
        x = "Cell-line pair",
        y = "Jaccard index",
        color = "Antibody comparison"
      ) +
      theme_minimal(base_size = 18)
  )
  
  dev.off()
  
  # 8. Mean Jaccard by antibody pair
  png(
    file = file.path(jaccard_grouped_output_dir, "barplot_mean_jaccard_by_antibody_pair.png"),
    width = 3000,
    height = 2200,
    res = 180
  )
  
  print(
    ggplot(
      jaccard_summary_by_antibody_pair %>% dplyr::filter(n_pairs >= 2),
      aes(x = reorder(antibody_pair, mean_jaccard), y = mean_jaccard)
    ) +
      geom_col(color = "black") +
      geom_errorbar(
        aes(
          ymin = pmax(mean_jaccard - sd_jaccard, 0),
          ymax = mean_jaccard + sd_jaccard
        ),
        width = 0.25
      ) +
      coord_flip() +
      labs(
        title = "Mean Jaccard index by antibody pair",
        subtitle = "Only antibody pairs with at least 2 pairwise comparisons shown",
        x = "Antibody pair",
        y = "Mean Jaccard index"
      ) +
      theme_minimal(base_size = 16)
  )
  
  dev.off()
  
  # 9. Per protein-motif: mean Jaccard by combined group
  png(
    file = file.path(jaccard_grouped_output_dir, "barplot_protein_motif_mean_jaccard_by_combined_group.png"),
    width = 3600,
    height = 2600,
    res = 180
  )
  
  print(
    ggplot(
      jaccard_summary_by_protein_motif_groups %>%
        dplyr::filter(combined_group != "Unknown antibody"),
      aes(
        x = reorder(paste(protein, motif, sep = " - "), mean_jaccard),
        y = mean_jaccard,
        fill = combined_group
      )
    ) +
      geom_col(position = "dodge", color = "black") +
      coord_flip() +
      labs(
        title = "Mean Jaccard index per protein-motif by cell-line/antibody grouping",
        x = "Protein - motif",
        y = "Mean Jaccard index",
        fill = "Group"
      ) +
      theme_minimal(base_size = 15)
  )
  
  dev.off()
  
  ############################ Console summary ##################################
  
  print("Grouped Jaccard analysis finished.")
  print(paste("Grouped analysis outputs saved to:", jaccard_grouped_output_dir))
  
  print("Summary by cell-line group:")
  print(jaccard_summary_by_biosample_group)
  
  print("Summary by antibody group:")
  print(jaccard_summary_by_antibody_group)
  
  print("Summary by combined group:")
  print(jaccard_summary_by_combined_group)
  
  
