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
        # Extract only chr, start_peak, end_peak columns
        peaks_i <- peak_data[[i]][, c("chr", "start_peak", "end_peak")]
        peaks_j <- peak_data[[j]][, c("chr", "start_peak", "end_peak")]
        
        # Remove rows with missing values
        peaks_i <- peaks_i[complete.cases(peaks_i), ]
        peaks_j <- peaks_j[complete.cases(peaks_j), ]
        
        # Calculate Jaccard
        jaccard_value <- calculate_jaccard_peaks_overlap(peaks_i, peaks_j)
        
        jaccard_matrix[i, j] <- jaccard_value
        jaccard_matrix[j, i] <- jaccard_value  # Symmetric matrix
      }
    }
    
    # Save Jaccard matrix as RDS
    saveRDS(
      jaccard_matrix,
      file = file.path(protein_motif_dir, "jaccard_matrix.rds")
    )
    
    # Convert matrix to long format data frame
    jaccard_df <- as.data.frame(jaccard_matrix)
    jaccard_df$experiment_1 <- rownames(jaccard_matrix)
    
    jaccard_long <- jaccard_df %>%
      pivot_longer(
        cols = -experiment_1,
        names_to = "experiment_2",
        values_to = "jaccard_index"
      ) %>%
      filter(experiment_1 < experiment_2) %>%  # Only keep upper triangle
      arrange(desc(jaccard_index))
    
    # Save as CSV
    write.csv(
      jaccard_long,
      file = file.path(protein_motif_dir, "jaccard_table.csv"),
      row.names = FALSE
    )
    
    # Calculate summary statistics
    jaccard_values <- jaccard_matrix[upper.tri(jaccard_matrix)]
    jaccard_values <- jaccard_values[!is.na(jaccard_values)]
    
    summary_stats <- data.frame(
      protein = protein,
      motif = motif,
      n_experiments = n_exp,
      n_pairwise_comparisons = length(jaccard_values),
      mean_jaccard = mean(jaccard_values, na.rm = TRUE),
      median_jaccard = median(jaccard_values, na.rm = TRUE),
      min_jaccard = min(jaccard_values, na.rm = TRUE),
      max_jaccard = max(jaccard_values, na.rm = TRUE),
      sd_jaccard = sd(jaccard_values, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    write.csv(
      summary_stats,
      file = file.path(protein_motif_dir, "jaccard_summary.csv"),
      row.names = FALSE
    )
    
    # Generate heatmap
    png(
      file = file.path(
        protein_motif_dir,
        paste0("jaccard_heatmap_", protein, "_", motif, ".png")
      ),
      width = 800,
      height = 700
    )
    
    pheatmap(
      jaccard_matrix,
      main = paste("Jaccard Index -", protein, "(", motif, ")"),
      display_numbers = TRUE,
      number_format = "%.3f",
      breaks = seq(0, 1, length.out = 101),
      color = colorRampPalette(c("white", "yellow", "orange", "red"))(100),
      cellwidth = 100,
      cellheight = 100,
      fontsize = 10,
      fontsize_number = 9
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
    file = file.path(output_folder, "ALL_PROTEINS_MOTIFS_jaccard_summary.png"),
    width = 1200,
    height = 800
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
}

end_script <- Sys.time()
print(paste("Total runtime:", round(difftime(end_script, start_script, units = "mins"), 2), "minutes"))
print(paste("Results saved to:", output_folder))
