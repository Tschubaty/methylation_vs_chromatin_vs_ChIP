#################################################################
##  Chip-seq vs methylation vs Chromatin
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

##################################### INPUT ########################################
# Define the input directory for chromosome files
input_dir <- "C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/WGBS/byChr/concat_methylation"

# Set up parallel processing
n_cores <- detectCores() - 1
registerDoParallel(cores = n_cores)


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

# Create a named vector to map chromatin states to their colors
chromatin_state_map <- setNames(chromatin_state_colors, chromatin_state_names)

# View the mapping (optional)
chromatin_state_map

##################################### Functions #################################

combined_df_chr <- NULL
chr_files <- file.path(input_dir,paste0(CHR_NAMES,"_cpgs_island_all_methylation.bed"))

for(chr_file in chr_files){   #chr_file
  print(chr_file)
  df_chr <- readr::read_delim(chr_file, delim = "\t", col_names = TRUE, show_col_types = FALSE)
  
  # Rename the `#Chromosome` column to `chr`
  df_chr <- df_chr %>%
    rename(chr = `#Chromosome`)  # Rename the column
  
  #  Replace "." with NA in all columns (character and numeric-like)
  df_chr <- df_chr %>%
    mutate(across(everything(), ~ ifelse(. == "." | .== "Not_in_CpG_Island", NA, .)))  # Replace "." with NA for all columns
  
  # Assuming `df_chr` is your data frame
  
  # Convert specific columns to factors and the rest to integers
  df_chr <- df_chr %>%
    # Convert Chromosome and Chromatin_State columns to factors
    mutate(
      across(chr, as.factor),  # Convert `#Chromosome` to factor
      across(starts_with("Chromatin_State"), as.factor)  # Convert all Chromatin_State columns to factors
    ) %>%
    # Convert remaining numeric columns to integers
    mutate(
      across(where(is.numeric), as.integer)
    )

  if (is.null(combined_df_chr)) {
    combined_df_chr <- df_chr
  } else {
    combined_df_chr <- bind_rows(combined_df_chr, df_chr)
  }
  
  
}


# Save the concatenated dataframe to a file
output_file <- file.path(input_dir, "combined_methylation_data.bed")
write_delim(combined_df_chr, output_file, delim = "\t")

cat("Final concatenated data saved to:", output_file, "\n")

# End script timing
end_script <- Sys.time()
cat("Script completed in:", end_script - start_script, "\n")

##################################### Main ######################################
combined_df_chr <- read_delim("C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/WGBS/byChr/concat_methylation/combined_methylation_data.bed")


# Define the directory where you want to save the plots
output_dir <- "C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/plots/"

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
################################### Figure 2a) #####################

# Number of chunks
num_chunks <- 1000

# Define the number of rows per chunk
total_rows <- nrow(combined_df_chr)
chunk_size <- ceiling(total_rows / num_chunks)

# Loop over each chunk to process sequentially
for (chunk_idx in 1:num_chunks) {
  cat("Processing chunk", chunk_idx, "of", num_chunks, "\n")
  
  # Define the row indices for the current chunk
  start_idx <- (chunk_idx - 1) * chunk_size + 1
  end_idx <- min(chunk_idx * chunk_size, total_rows)
  
  # Subset the data for the current chunk
  df_chunk <- combined_df_chr[start_idx:end_idx, ]
  
  # Recalculate fRead for each biosample in the chunk
  df_chunk <- df_chunk %>%
    mutate(
      fRead_A549 = ifelse(nRead_A549 == 0, NA, as.double(mRead_A549 / nRead_A549)),
      fRead_GM12878 = ifelse(nRead_GM12878 == 0, NA, as.double(mRead_GM12878 / nRead_GM12878)),
      fRead_HepG2 = ifelse(nRead_HepG2 == 0, NA, as.double(mRead_HepG2 / nRead_HepG2)),
      fRead_K562 = ifelse(nRead_K562 == 0, NA, as.double(mRead_K562 / nRead_K562))
    )
  
  # Reshape the data: gather fRead, mRead, nRead, and Chromatin_State into long format
  df_long_fRead_chunk <- df_chunk %>%
    pivot_longer(
      cols = starts_with("fRead"),  # Gather all fRead columns
      names_to = "Biosample",        # Create a new column called "Biosample"
      values_to = "fRead"            # Store the recalculated fRead values
    ) %>%
    mutate(
      # Remove the "fRead_" prefix from the biosample names
      Biosample = gsub("fRead_", "", Biosample)
    ) %>%
    pivot_longer(
      cols = c(starts_with("mRead"), starts_with("nRead"), starts_with("Chromatin_State")),
      names_pattern = "(.*)_(.*)",  # Extract variable name and biosample name
      names_to = c(".value", "BiosampleMatch")  # `.value` splits into the respective columns
    ) %>%
    filter(Biosample == BiosampleMatch) %>%  # Keep rows where Biosample matches the specific columns
    select(chr, Start, End, CpG_Island_Status, Strand, fRead, mRead, nRead, Chromatin_State, Biosample)  # Keep only the required columns
  
  # Save each chunk directly to disk
  chunk_filename <- file.path(output_dir, paste0("df_long_fRead_chunk_", chunk_idx, ".csv"))
  write.csv(df_long_fRead_chunk, chunk_filename, row.names = FALSE)
  
  # Remove the temporary variables and run garbage collection
  rm(df_chunk, df_long_fRead_chunk)
  gc()
  
  # Optional: Print progress
  cat("Processed and saved chunk", chunk_idx, "of", num_chunks, "\n")
}

# Now, after all chunks are saved, combine the files from disk

# Get a list of all chunk files
chunk_files <- list.files(output_dir, pattern = "df_long_fRead_chunk_.*\\.csv", full.names = TRUE)

# Combine the files incrementally
df_long_combined <- do.call(rbind, lapply(chunk_files, read.csv))

# Save the combined data to disk
write.csv(df_long_combined, file = file.path(output_dir, "df_long_combined.csv"), row.names = FALSE)

# Cleanup: Remove only the temporary variables created in this script, not everything
rm(chunk_idx, chunk_size, total_rows, num_chunks, chunk_files)
gc()

cat("Processing complete and combined dataframe saved to", file.path(output_dir, "df_long_combined.csv"), "\n")


# Create violin plot using the combined dataframe
violin_plot <- ggplot(df_long_fRead_combined, aes(x = Biosample, y = fRead, fill = Biosample)) +
  geom_violin(trim = TRUE) +
  labs(
    title = "Distribution of fRead by Biosample",
    x = "Biosample",
    y = "fRead"
  ) +
  theme_minimal()

# Display the violin plot
print(violin_plot)

# Save the violin plot to disk
ggsave(
  filename = file.path(output_dir, "fRead_violin_plot.png"),
  plot = violin_plot,
  width = 10, height = 6, units = "in"
)


#################################################################

# Reshape the data: gather chromatin state columns into long format
df_long <- combined_df_chr %>%
  pivot_longer(
    cols = starts_with("Chromatin_State"), 
    names_to = "Biosample", 
    values_to = "Chromatin_State"
  ) %>%
  mutate(
    # Remove the "Chromatin_State_" prefix from the biosample names
    Biosample = gsub("Chromatin_State_", "", Biosample)
  )

# Plot the distribution of chromatin states, with dodge by biosample
chromatin_plot <- ggplot(df_long, aes(x = Chromatin_State, fill = Biosample)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = group.colors) +
  labs(
    title = "Distribution of Chromatin States by Biosample",
    x = "Chromatin State",
    y = "Count",
    fill = "Biosample"
  ) +
  theme_minimal()

# Save the chromatin state distribution plot
ggsave(
  filename = file.path(output_dir, "chromatin_state_distribution.png"),
  plot = chromatin_plot,
  width = 10, height = 6, units = "in"
)

rm(chromatin_plot, df_long)
gc()  # Run garbage collection to free up memory

# Reshape the data: gather fRead columns into long format
df_long_fRead <- combined_df_chr %>%
  pivot_longer(
    cols = starts_with("fRead"),  # Gather all fRead columns
    names_to = "Biosample",        # Create a new column called "Biosample"
    values_to = "fRead"            # Store the values in a column called "fRead"
  ) %>%
  mutate(
    # Remove the "fRead_" prefix from the biosample names
    Biosample = gsub("fRead_", "", Biosample),
    
    # Convert Chromatin_State columns (for each biosample) to factors
    across(starts_with("Chromatin_State"), as.factor)
  )

# Assuming you want to use one chromatin state column for plotting (e.g., Chromatin_State_A549)
fread_boxplot <- ggplot(df_long_fRead, aes(x = Chromatin_State_A549, y = fRead, fill = Biosample)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +
  scale_fill_manual(values = group.colors) +  # Use predefined colors for biosamples
  labs(
    title = "Distribution of fRead by Chromatin State and Biosample",
    x = "Chromatin State",
    y = "fRead",
    fill = "Biosample"
  ) +
  theme_minimal()

# Save the fRead boxplot
ggsave(
  filename = file.path(output_dir, "fRead_boxplot.png"),
  plot = fread_boxplot,
  width = 10, height = 6, units = "in"
)

rm(fread_boxplot, df_long_fRead)
gc()  # Run garbage collection to free up memory