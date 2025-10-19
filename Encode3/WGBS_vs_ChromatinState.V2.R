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
input_dir <- file.path(this.dir,"WGBS/byChr/concat_methylation")

# Set up parallel processing
n_cores <- detectCores() - 1
registerDoParallel(cores = n_cores)


# Define the directory where you want to save the plots
output_dir <- file.path(this.dir,"plots")

# Ensure the output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


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


######################## preprocessing long format start ###########################
if (FALSE) {
  
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
  
  #combined_df_chr <- read_delim(file = file.path(input_dir,"combined_methylation_data.bed"))
  #combined_df_chr <- readRDS(file = file.path(input_dir,"combined_methylation_data_compact.rds"))

  # ----- BEGIN: CpG chromatin state 1 summaries -----
  cell_lines <- c("A549", "GM12878", "HepG2", "K562")
  state_cols <- paste0("Chromatin_State_", cell_lines)
  total_cpgs <- nrow(combined_df_chr)
  percent_state1 <- sapply(state_cols, function(col) {
    sum(combined_df_chr[[col]] == 1, na.rm = TRUE) / total_cpgs * 100
  })
  print("Percentage of CpGs in state 1 in each cell line:")
  print(percent_state1)
  
  pairwise_stable <- matrix(NA, nrow=length(cell_lines), ncol=length(cell_lines), dimnames=list(cell_lines, cell_lines))
  for (i in seq_along(cell_lines)) {
    for (j in seq_along(cell_lines)) {
      if (i != j) {
        col1 <- state_cols[i]
        col2 <- state_cols[j]
        idx1 <- which(combined_df_chr[[col1]] == 1)
        stable_count <- sum(combined_df_chr[idx1, col2] == 1, na.rm = TRUE)
        pairwise_stable[i, j] <- stable_count / length(idx1) * 100
      }
    }
  }
  print("Pairwise percent of 'state 1' CpGs in cell line 1 that are also 'state 1' in cell line 2:")
  print(pairwise_stable)
  # ----- END: CpG chromatin state 1 summaries -----
  
  # bar plto of concervation fo states:
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  
  cell_lines <- c("A549", "GM12878", "HepG2", "K562")
  state_cols <- paste0("Chromatin_State_", cell_lines)
  states <- as.character(1:18)
  
  results <- data.frame()
  
  for (i in seq_along(cell_lines)) {
    for (j in seq_along(cell_lines)) {
      if (i != j) {
        col1 <- state_cols[i]
        col2 <- state_cols[j]
        for (state in states) {
          idx_state1 <- which(combined_df_chr[[col1]] == state)
          if (length(idx_state1) > 0) {
            n_conserved <- sum(combined_df_chr[idx_state1, col2] == state, na.rm=TRUE)
            percent_conserved <- n_conserved / length(idx_state1) * 100
          } else {
            percent_conserved <- NA  # avoid division by zero
          }
          results <- rbind(results, data.frame(
            State = state,
            CellLine1 = cell_lines[i],
            CellLine2 = cell_lines[j],
            PercentConserved = percent_conserved
          ))
        }
      }
    }
  }
  
  # Ensure your chromatin_state_colors vector is named (if not already)
  names(chromatin_state_colors) <- as.character(1:18)
  
  results$State <- factor(results$State, levels = as.character(1:18))
  
  ggplot(results, aes(x = State, y = PercentConserved, fill = State)) +
    geom_bar(stat = "identity") +
    facet_grid(CellLine1 ~ CellLine2) +
    theme_bw() +
    scale_fill_manual(
      values = chromatin_state_colors,
      labels = chromatin_state_names
    ) +
    labs(
      title = "Percentage of CpGs with conserved chromatin state (by state, between cell line pairs)",
      x = "Chromatin State",
      y = "% Conserved"
    ) +
    theme(axis.text.x = element_blank())
    #theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
  
  #save as : Barplot of CpG Chromatin State Conservation Across Cell Line Pairs
  
  # barplots finished 

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
        fRead_GM12878 = ifelse(
          nRead_GM12878 == 0,
          NA,
          as.double(mRead_GM12878 / nRead_GM12878)
        ),
        fRead_HepG2 = ifelse(nRead_HepG2 == 0, NA, as.double(mRead_HepG2 / nRead_HepG2)),
        fRead_K562 = ifelse(nRead_K562 == 0, NA, as.double(mRead_K562 / nRead_K562))
      )
    
    # Reshape the data: gather fRead, mRead, nRead, and Chromatin_State into long format
    df_long_fRead_chunk <- df_chunk %>%
      pivot_longer(
        cols = starts_with("fRead"),
        # Gather all fRead columns
        names_to = "Biosample",
        # Create a new column called "Biosample"
        values_to = "fRead"            # Store the recalculated fRead values
      ) %>%
      mutate(# Remove the "fRead_" prefix from the biosample names
        Biosample = gsub("fRead_", "", Biosample)) %>%
      pivot_longer(
        cols = c(
          starts_with("mRead"),
          starts_with("nRead"),
          starts_with("Chromatin_State")
        ),
        names_pattern = "(.*)_(.*)",
        # Extract variable name and biosample name
        names_to = c(".value", "BiosampleMatch")  # `.value` splits into the respective columns
      ) %>%
      filter(Biosample == BiosampleMatch) %>%  # Keep rows where Biosample matches the specific columns
      select(
        chr,
        Start,
        End,
        CpG_Island_Status,
        Strand,
        fRead,
        mRead,
        nRead,
        Chromatin_State,
        Biosample
      )  # Keep only the required columns
    
    # Save each chunk directly to disk
    chunk_filename <- file.path(output_dir,
                                paste0("df_long_fRead_chunk_", chunk_idx, ".csv"))
    write.csv(df_long_fRead_chunk, chunk_filename, row.names = FALSE)
    
    # Remove the temporary variables and run garbage collection
    rm(df_chunk, df_long_fRead_chunk)
    gc()
    
    # Optional: Print progress
    cat("Processed and saved chunk",
        chunk_idx,
        "of",
        num_chunks,
        "\n")
  }
  
  # Now, after all chunks are saved, combine the files from disk
  
  # Get a list of all chunk files
  chunk_files <- list.files(output_dir, pattern = "df_long_fRead_chunk_.*\\.csv", full.names = TRUE)
  
  # Combine the files incrementally
  df_long_combined <- do.call(rbind, lapply(chunk_files, read.csv))
  
  # Save the combined data to disk
  #write.csv(df_long_combined, file = file.path(output_dir, "df_long_combined.csv"), row.names = FALSE)
  
  # Cleanup: Remove only the temporary variables created in this script, not everything
  rm(chunk_idx, chunk_size, total_rows, num_chunks, chunk_files)
  gc()
  
  cat(
    "Processing complete and combined dataframe saved to",
    file.path(output_dir, "df_long_combined.csv"),
    "\n"
  )
  
  # Convert character columns to factors
  # Specify level order for chr, Strand, and Biosample
  df_long_combined$chr <- factor(df_long_combined$chr, levels = CHR_NAMES)
  df_long_combined$Strand <- factor(df_long_combined$Strand, levels = c("+", "-"))
  df_long_combined$Biosample <- factor(df_long_combined$Biosample, levels = names(group.colors))
  
  saveRDS(object = df_long_combined,
          file = file.path(output_dir, "df_long_combined.RDS"))
}
######################## preprocessing long format ended #######################
if(FALSE){
  df_long_combined <- readRDS(file = file.path(output_dir, "df_long_combined.RDS"))
  # Create violin plot using the combined dataframe
  
  df_long_Chromatin_State_Biosample_fRead <- df_long_combined[,c("Chromatin_State","Biosample","fRead")]
  rm(df_long_combined)
  gc()
  
  saveRDS(object = df_long_Chromatin_State_Biosample_fRead,
          file = file.path(output_dir, "df_long_Chromatin_State_Biosample_fRead.RDS"))
}
##################################### Main ######################################
df_long_Chromatin_State_Biosample_fRead <- readRDS(file = file.path(output_dir, 
                                                                    "df_long_Chromatin_State_Biosample_fRead.RDS"))
df_long_Chromatin_State_Biosample_fRead$Chromatin_State <- factor(df_long_Chromatin_State_Biosample_fRead$Chromatin_State, levels =seq(1,18))
#

library(dplyr)

# Summarize data with additional columns
genome_summary_table_with_meth_data <- df_long_Chromatin_State_Biosample_fRead %>%
  filter(!is.na(fRead)) %>% # Include only non-NA entries
  group_by(Biosample) %>%
  mutate(total_entries = n()) %>% # Total non-NA entries per biosample
  group_by(Chromatin_State, Biosample) %>%
  summarise(
    total_entries = first(total_entries), # Consistent total non-NA entries across states
    count_entries = n(), # Non-NA entries for the current state
    mean_fRead = mean(fRead, na.rm = TRUE),
    percentage_entries = count_entries / total_entries * 100
  )

# View the result
print(genome_summary_table_with_meth_data)

saveRDS(object = genome_summary_table_with_meth_data,
        file = file.path(output_dir, "genome_summary_table_with_meth_data.RDS"))


genome_summary_table <- readRDS(file = file.path(output_dir,"genome_summary_table.rds"))

################################### Figure 2a) ################################
# Load necessary libraries
library(ggplot2)
library(ggridges)

# Ensure no NA values in fRead
df_filtered <- df_long_Chromatin_State_Biosample_fRead[!is.na(df_long_Chromatin_State_Biosample_fRead$fRead), ]
rm("df_long_Chromatin_State_Biosample_fRead")
gc()


# Create the ridge plot
ggplot(df_filtered, aes(x = fRead, y = Biosample, fill = Biosample)) +
  geom_density_ridges(alpha = 0.7) +
  scale_fill_manual(values = group.colors) +
  labs(
    title = "Distribution of fRead per Biosample",
    x = "fRead",
    y = "Biosample"
  ) +
  theme_minimal()

##############Histogram vertsion ########################
library(patchwork)

# Define variables
n_breaks <- 50  # Number of bins for histograms
x_limits <- c(0, 1)  # Consistent x-axis limits

# Calculate y-axis limits based on the highest count across all histograms
y_limits <- c(0, max(sapply(names(group.colors), function(sample) {
  max(ggplot_build(
    ggplot(df_filtered[df_filtered$Biosample == sample, ], aes(x = fRead)) +
      geom_histogram(breaks = seq(x_limits[1], x_limits[2], length.out = n_breaks + 1))
  )$data[[1]]$count)
})))



# Create individual histograms for each sample
plots <- lapply(names(group.colors), function(sample) {
  
  # Filter data for the sample
  sample_data <- df_filtered[df_filtered$Biosample == sample, ]
  
  
  sample_quantiles <- quantile(sample_data$fRead, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  
  ggplot(sample_data, aes(x = fRead, fill = Biosample)) +
    geom_histogram(
      breaks = seq(x_limits[1], x_limits[2], length.out = n_breaks + 1),
      alpha = 0.7,
    ) +
    geom_vline(aes(xintercept = sample_quantiles[1]), color = "red", linetype = "dotted", size = 0.8) +  # 25th percentile
    geom_vline(aes(xintercept = sample_quantiles[2]), color = "green", linetype = "dotted", size = 0.8) +  # Median (50th percentile)
    geom_vline(aes(xintercept = sample_quantiles[3]), color = "blue", linetype = "dotted", size = 0.8) +  # 75th percentile
    scale_fill_manual(values = group.colors[sample]) +
    labs(
      title = paste("Histogram of CpG methylation -", sample),
      x = "CpG methylation",
      y = "CpG Count"
    ) +
    theme_minimal() +
    xlim(x_limits) +  # Set consistent x-axis limits
    ylim(y_limits)    # Set consistent y-axis limits
})

# Combine plots into one layout using wrap_plots
combined_plot <- wrap_plots(plots) +
  plot_annotation(
    title = "Histogram of CpG methylation per Biosample",
    theme = theme_minimal()
  )

# Print the combined plot
print(combined_plot)

ggsave(
  filename = file.path(output_dir, "combined_plot_CpG_methylation.png"),
  plot = combined_plot,
  width = 12,   # Set the desired width
  height = 8,   # Set the desired height
  dpi = 300     # Set the resolution for the image
)


######################### log scale idea 

library(scales)  # make sure it's loaded

hist_plot <- ggplot(df_filtered, aes(x = fRead)) +
  stat_bin(bins = 200, fill = "gray80", color = NA) +
  scale_y_continuous(
    trans = "log10",
    labels = label_number(scale_cut = cut_si("")),  # <-- replacement
    expand = expansion(mult = c(0, 0.08))           # a bit of headroom for labels
  ) +
  coord_cartesian(xlim = c(0,1)) +
  facet_wrap(~ Biosample, ncol = 2) +
  geom_vline(data = qdf, aes(xintercept = xint, color = q),
             linetype = "dashed", linewidth = 0.5, show.legend = FALSE) +
  geom_text(data = qdf, aes(x = xint, y = Inf, label = qlab, color = q),
            vjust = 1.2, size = 2.8, show.legend = FALSE) +
  scale_color_manual(values = c(q25="#D55E00", q50="#009E73", q75="#0072B2")) +
  labs(title = "CpG methylation per biosample",
       x = "CpG methylation (fRead)", y = "CpG count (log10)") +
  theme_minimal(base_size = 11) +
  theme(panel.spacing = unit(8, "pt"),
        strip.text = element_text(face = "bold"),
        axis.title.y = element_text(margin = margin(r = 8)),
        axis.title.x = element_text(margin = margin(t = 6)))



################### barplot states #######################################################
rm("data")
rm("barplot")
rm("combined_plot")
rm("plots")
gc()

# Create a bar plot
barplot <- ggplot(df_filtered, aes(x = Chromatin_State, fill = Biosample)) +
  geom_bar(position = "dodge") +  # Dodge to separate bars by Biosample
  scale_fill_manual(values = group.colors) +  # Apply the custom colors for Biosamples
  labs(
    title = "Counts of Chromatin States for Each Biosample",
    x = "Chromatin State",
    y = "CpG Count",
    fill = "Biosample"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(output_dir, "barplot_chromatin_state.png"),
  plot = barplot,
  width = 12,   # Set the desired width
  height = 8,   # Set the desired height
  dpi = 300     # Set the resolution for the image
)
rm(barplot)
gc()


library(scales)

# order states
dfp <- df_filtered %>%
  dplyr::mutate(Chromatin_State = factor(Chromatin_State, levels = as.character(1:18)))

barplot_pub <- ggplot(dfp, aes(x = Chromatin_State, fill = Biosample)) +
  geom_bar(position = position_dodge(width = 0.85), width = 0.72) +
  scale_fill_manual(values = group.colors, name = "Biosample") +
  scale_y_continuous(labels = label_number(big.mark = ","), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Counts of Chromatin States by Biosample",
    x = "ChromHMM state",
    y = "CpG count"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(face = "bold"),          # <-- make ALL text bold
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    plot.title  = element_text(size = 14)
  )

ggsave(
  filename = file.path(output_dir, " "),
  plot = barplot_pub,
  width = 11, height = 6, dpi = 300, limitsize = FALSE
)
rm(barplot_pub)
gc()


###############################################################################
library(dplyr)
library(ggplot2)
library(ggalluvial)
library(tidyr)
library(readr)

# Load data
combined_df_chr <- read_delim(file = file.path(input_dir, "combined_methylation_data.bed"))

# Identify columns with "Chromatin_State" in their names
chromatin_cols <- names(combined_df_chr)[grepl("Chromatin_State", names(combined_df_chr))]

# Convert all chromatin state columns to factors
combined_df_chr <- combined_df_chr %>%
  mutate(across(all_of(chromatin_cols), as.factor))

# Determine a unified set of levels, sorted numerically
unified_levels <- sort(as.numeric(unique(unlist(lapply(combined_df_chr[chromatin_cols], levels)))))

# Standardize all Chromatin_State columns to the unified levels
combined_df_chr <- combined_df_chr %>%
  mutate(across(all_of(chromatin_cols), ~ factor(.x, levels = as.character(unified_levels))))

# Define chromatin state colors
chromatin_state_colors <- c(
  "1" = "#FF0000",  "2" = "#FF4500",  "3" = "#FF9900",  "4" = "#FFCC00",
  "5" = "#00CC00",  "6" = "#006400",  "7" = "#FFD700",  "8" = "#FFD700",
  "9" = "#FFFF00",  "10" = "#FFDD00", "11" = "#FFEA73", "12" = "#9370DB",
  "13" = "#C0C0C0", "14" = "#FF4500", "15" = "#FFDD00", "16" = "#808080",
  "17" = "#A9A9A9", "18" = "#000000"
)

# Recalculate fRead as mRead / nRead for each biosample and drop mRead/nRead columns
combined_df_chr <- combined_df_chr %>%
  mutate(
    fRead_A549 = ifelse(nRead_A549 > 0, mRead_A549 / nRead_A549, NA),
    fRead_GM12878 = ifelse(nRead_GM12878 > 0, mRead_GM12878 / nRead_GM12878, NA),
    fRead_HepG2 = ifelse(nRead_HepG2 > 0, mRead_HepG2 / nRead_HepG2, NA),
    fRead_K562 = ifelse(nRead_K562 > 0, mRead_K562 / nRead_K562, NA)
  ) %>%
  select(-mRead_A549, -nRead_A549, 
         -mRead_GM12878, -nRead_GM12878, 
         -mRead_HepG2, -nRead_HepG2, 
         -mRead_K562, -nRead_K562)

# Convert Start and End columns to integers
combined_df_chr <- combined_df_chr %>%
  mutate(
    Start = as.integer(Start),
    End = as.integer(End)
  )

# Drop the CpG_Island_Status column
combined_df_chr <- combined_df_chr %>%
  select(-CpG_Island_Status)

# View the updated dataset
head(combined_df_chr)


# Define the output file path
output_file <- file.path(input_dir, "combined_methylation_data_compact.rds")

# Save the processed data as an RDS file
saveRDS(combined_df_chr, file = output_file)

# Confirm the file was saved
message("Data saved to: ", output_file)

# Loop through all combinations of chromatin state columns
combinations <- combn(chromatin_cols, 2, simplify = FALSE)
for (combo in combinations) {
  source_col <- combo[1]
  target_col <- combo[2]
  
  # Sankey diagram
  sankey_data <- combined_df_chr %>%
    count(!!sym(source_col), !!sym(target_col)) %>%
    rename(Source = !!sym(source_col), Target = !!sym(target_col), Count = n)
  
  sankey_plot <- ggplot(sankey_data, aes(axis1 = Source, axis2 = Target, y = Count)) +
    geom_alluvium(aes(fill = Source)) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_fill_manual(values = chromatin_state_colors) +
    labs(
      title = paste("Chromatin State Transitions Between", source_col, "and", target_col),
      x = source_col,
      y = target_col,
      fill = "Chromatin State"
    ) +
    theme_minimal()
  
  ggsave(
    filename = file.path(output_dir, paste0("sankey_", source_col, "_to_", target_col, ".png")),
    plot = sankey_plot,
    width = 12,
    height = 8,
    dpi = 300
  )
  
  # Heatmap
  total_counts <- combined_df_chr %>%
    count(!!sym(source_col)) %>%
    rename(Total = n, Source = !!sym(source_col))
  
  transition_matrix <- combined_df_chr %>%
    count(!!sym(source_col), !!sym(target_col)) %>%
    rename(Source = !!sym(source_col), Target = !!sym(target_col), Count = n) %>%
    left_join(total_counts, by = "Source") %>%
    mutate(Percentage = Count / Total) %>%
    complete(Source = as.factor(unified_levels), Target = as.factor(unified_levels), fill = list(Percentage = 0))
  
  heatmap_plot <- ggplot(transition_matrix, aes(x = Source, y = Target, fill = Percentage)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "blue", labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = paste("Heatmap of Percentage-Based Transitions Between", source_col, "and", target_col),
      x = source_col,
      y = target_col,
      fill = "Percentage"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    filename = file.path(output_dir, paste0("heatmap_", source_col, "_to_", target_col, ".png")),
    plot = heatmap_plot,
    width = 12,
    height = 8,
    dpi = 300
  )
}
###############################################################################
library(ggplot2)
library(gridExtra)

# Extract all biosample columns with chromatin states
chromatin_cols <- names(combined_df_chr)[grepl("Chromatin_State", names(combined_df_chr))]

# Extract all unique chromatin states across biosamples
all_states <- sort(unique(unlist(lapply(chromatin_cols, function(col) unique(combined_df_chr[[col]])))))

# Loop through all chromatin state combinations
for (stateX in all_states) {
  for (stateY in all_states) {
    # Initialize a list to store plots
    plot_list <- list()
    for (biosample_A in chromatin_cols) {
      for (biosample_B in chromatin_cols) {
        if (biosample_A == biosample_B) next  # Skip self-comparisons
        print(paste(stateX ,stateY,biosample_A,biosample_B))
        # Debug example
        # stateX <- all_states[1]
        # stateY <- all_states[2]
        # biosample_A <- chromatin_cols[1]
        # biosample_B <- chromatin_cols[2]

        # Filter data for the current chromatin state and biosample combination
        filtered_data <- combined_df_chr %>%
          filter(
            !!sym(biosample_A) == stateX & !!sym(biosample_B) == stateY
          ) %>%
          select(
            chr, Start, End,
            !!sym(biosample_A), !!sym(biosample_B),
            paste0("fRead_", gsub("Chromatin_State_", "", biosample_A)),
            paste0("fRead_", gsub("Chromatin_State_", "", biosample_B))
          ) %>%
          rename(
            fRead_A = paste0("fRead_", gsub("Chromatin_State_", "", biosample_A)),
            fRead_B = paste0("fRead_", gsub("Chromatin_State_", "", biosample_B))
          ) %>%
          mutate(
            Methylation_Difference = fRead_A - fRead_B,
            Biosample_A = biosample_A,
            Biosample_B = biosample_B
          )
        
        # Skip if no data for this combination
        if (nrow(filtered_data) == 0) next
        
        # Create the histogram of methylation differences
        p <- ggplot(filtered_data, aes(x = Methylation_Difference)) +
          geom_histogram(binwidth = 0.01, fill = "steelblue", color = "black", alpha = 0.7) +
          labs(
            title = paste("Methylation Difference:", stateX, "->", stateY, "(", biosample_A, "vs", biosample_B, ")"),
            x = "Methylation Difference",
            y = "counts"
          ) +
          theme_minimal() +
          theme(
            plot.title = element_text(size = 10),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8)
          )
        
        # Save the plot in the plot list
        plot_list[[paste(stateX, stateY, biosample_A, biosample_B, sep = "_")]] <- p
      }
    }
    # Combine all plots into a single panel
    combined_plots <- do.call(grid.arrange, c(plot_list, ncol = 3))
    
    # Save the combined plot
    ggsave(
      filename = file.path(output_dir, paste0("methylation_difference_histograms_state_",stateX, "_into_", stateY,".png")),
      plot = combined_plots,
      width = 24,
      height = 18,
      dpi = 300
    )
  }
}


########################################################################################

library(ggplot2)
library(dplyr)

# Inputs
biosample_A <- "Chromatin_State_A549"  # Replace with desired biosample A column
biosample_B <- "Chromatin_State_K562"  # Replace with desired biosample B column
stateX <- "1"  # Replace with desired state in A
stateY <- "2"  # Replace with desired state in B

# Filter data for the specific states
filtered_data <- combined_df_chr %>%
  filter(!!sym(biosample_A) == stateX & !!sym(biosample_B) == stateY) %>%
  select(
    !!sym(biosample_A),
    !!sym(biosample_B),
    paste0("fRead_", gsub("Chromatin_State_", "", biosample_A)),
    paste0("fRead_", gsub("Chromatin_State_", "", biosample_B))
  )


# Rename columns for clarity
filtered_data <- filtered_data %>%
  rename(
    Chromatin_State_A = !!sym(biosample_A),
    Chromatin_State_B = !!sym(biosample_B),
    fRead_A = paste0("fRead_", gsub("Chromatin_State_", "", biosample_A)),
    fRead_B = paste0("fRead_", gsub("Chromatin_State_", "", biosample_B))
  )



p <- ggplot(filtered_data, aes(x = fRead_A, y = fRead_B)) +
  geom_hex() +
  scale_fill_gradient(
    trans = "log",  # Apply log scale
    low = "white", 
    high = "blue", 
    name = "Log Count"
  ) +
  labs(
    title = paste("Hexbin Plot (Log Scale): fRead Levels (State", stateX, "in", biosample_A, "vs. State", stateY, "in", biosample_B, ")"),
    x = paste("fRead in", biosample_A, "(State", stateX, ")"),
    y = paste("fRead in", biosample_B, "(State", stateY, ")")
  ) +
  theme_minimal()

# Save the plot
ggsave(
  filename = file.path(output_dir, paste0("hexbin_plot_logscale_", biosample_A, "_", stateX, "_to_", biosample_B, "_", stateY, ".png")),
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)
######################### loop all state combos ################################
library(ggplot2)
library(dplyr)
library(tidyr)

# Extract all biosample columns with chromatin states
chromatin_cols <- names(combined_df_chr)[grepl("Chromatin_State", names(combined_df_chr))]

# Extract all unique chromatin states across biosamples
all_states <- sort(unique(unlist(lapply(chromatin_cols, function(col) unique(combined_df_chr[[col]])))))

# Loop through all chromatin state combinations
for (stateX in all_states) {
  for (stateY in all_states) {
    for (biosample_A in chromatin_cols) {
      for (biosample_B in chromatin_cols) {
        if (biosample_A == biosample_B) next  # Skip self-comparisons
        
        # # debug example 
        # stateX <- all_states[1]
        # stateY <- all_states[2]
        # biosample_A <- chromatin_cols[1]
        # biosample_B <- chromatin_cols[2]

        
        
        # Filter data for the current chromatin state and biosample combination
        filtered_data <- combined_df_chr %>%
          filter(!!sym(biosample_A) == stateX & !!sym(biosample_B) == stateY) %>%
          select(
            "chr",
            "Start",
            "End",
            !!sym(biosample_A),
            !!sym(biosample_B),
            paste0("fRead_", gsub("Chromatin_State_", "", biosample_A)),
            paste0("fRead_", gsub("Chromatin_State_", "", biosample_B))
          ) %>%
          rename(
            fRead_A = paste0("fRead_", gsub("Chromatin_State_", "", biosample_A)),
            fRead_B = paste0("fRead_", gsub("Chromatin_State_", "", biosample_B))
          ) %>%
          mutate(
            Biosample_A = biosample_A,
            Biosample_B = biosample_B
          )
        
        # Skip if no data for this combination
        if (nrow(filtered_data) == 0) next
        
        # Create the plot
        p <- ggplot(filtered_data, aes(x = fRead_A, y = fRead_B)) +
          geom_hex() +
          scale_fill_gradient(
            trans = "log",  # Apply log scale
            low = "white", 
            high = "blue", 
            name = "Log Count"
          ) +
          labs(
            title = paste("Hexbin Plot (Log Scale): State", stateX,"methylation", 
                          "in", gsub(pattern = "Chromatin_State_",replacement = "",x = biosample_A),
                          "vs. State", stateY, "methylation",
                          "in", gsub(pattern = "Chromatin_State_",replacement = "",x = biosample_B)),
            x = paste("CpG methylation", gsub(pattern = "Chromatin_State_",replacement = "",x = biosample_A), "(State", stateX, ")"),
            y = paste("CpG methylation", gsub(pattern = "Chromatin_State_",replacement = "",x = biosample_B), "(State", stateY, ")")
          ) +
          theme_minimal()
        
        
        # Save the plot
        plot_filename <- paste0("state_", stateX, "_to_", stateY, "_", 
                                gsub("Chromatin_State_", "", biosample_A), 
                                "_vs_", 
                                gsub("Chromatin_State_", "", biosample_B), 
                                ".png")
        
        ggsave(
          filename = file.path(output_dir, plot_filename),
          plot = p,
          width = 8,
          height = 6,
          dpi = 300
        )
        
        # Optionally print progress
        print(paste("Generated plot for State", stateX, "to", stateY, ":", biosample_A, "vs.", biosample_B))
      }
    }
  }
}

# ################### library(ggalluvial) #####################################
# 
# library(ggalluvial)
# 
# #load data
# combined_df_chr <- read_delim(file = file.path(input_dir,"combined_methylation_data.bed"))
# 
# 
# # make states as factors 
# 
# # Step 1: Identify columns with "Chromatin_State" in their names
# chromatin_cols <- names(combined_df_chr)[grepl("Chromatin_State", names(combined_df_chr))]
# 
# # Step 2: Convert all columns to factors if they are not already factors
# combined_df_chr <- combined_df_chr %>%
#   mutate(across(all_of(chromatin_cols), as.factor))
# 
# # Step 3: Determine a unified set of levels, sorted numerically
# unified_levels <- sort(as.numeric(unique(unlist(lapply(combined_df_chr[chromatin_cols], levels)))))
# 
# # Step 4: Standardize all Chromatin_State columns to the unified levels
# combined_df_chr <- combined_df_chr %>%
#   mutate(across(all_of(chromatin_cols), ~ factor(.x, levels = as.character(unified_levels))))
# 
# # Define corresponding hex color codes for each chromatin state
# chromatin_state_colors <- c(
#   "1" = "#FF0000",  # TssA       -> Red
#   "2" = "#FF4500",  # TssFlnk    -> Orange-Red
#   "3" = "#FF9900",  # TssFlnkU   -> Orange
#   "4" = "#FFCC00",  # TssFlnkD   -> Yellow-Orange
#   "5" = "#00CC00",  # Tx         -> Green
#   "6" = "#006400",  # TxWk       -> Dark Green
#   "7" = "#FFD700",  # EnhG1      -> Gold
#   "8" = "#FFD700",  # EnhG2      -> Gold
#   "9" = "#FFFF00",  # EnhA1      -> Yellow
#   "10" = "#FFDD00", # EnhA2      -> Yellow-Orange
#   "11" = "#FFEA73", # EnhWk      -> Light Yellow
#   "12" = "#9370DB", # ZNF_Rpts   -> Purple
#   "13" = "#C0C0C0", # Het        -> Light Gray
#   "14" = "#FF4500", # TssBiv     -> Orange-Red (similar to TssFlnk)
#   "15" = "#FFDD00", # EnhBiv     -> Yellow-Orange (similar to EnhA2)
#   "16" = "#808080", # ReprPC     -> Dark Gray
#   "17" = "#A9A9A9", # ReprPCWk   -> Light Gray
#   "18" = "#000000"  # Quies      -> Black
# )
# 
# # Step 5: Prepare the data for Sankey plotting
# sankey_data <- combined_df_chr %>%
#   count(Chromatin_State_HepG2, Chromatin_State_K562)
# 
# # Step 6: Plot the Sankey diagram with custom colors
# p <- ggplot(sankey_data, aes(axis1 = Chromatin_State_HepG2, axis2 = Chromatin_State_K562, y = n)) +
#   geom_alluvium(aes(fill = Chromatin_State_HepG2)) +
#   geom_stratum() +
#   geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
#   scale_fill_manual(values = chromatin_state_colors) +
#   labs(title = "Chromatin State Transitions Between HepG2 and K562",
#        x = "Chromatin State",
#        y = "Count",
#        fill = "Chromatin State") +
#   theme_minimal()
# 
# 
# ggsave(
#   filename = file.path(output_dir, "sankey.png"),
#   plot = p,
#   width = 12,   # Set the desired width
#   height = 8,   # Set the desired height
#   dpi = 300     # Set the resolution for the image
# )
# 
# ############# heatmaop tranissiton matrrix ####################################
# # Step 1: Calculate total transitions for each starting state
# total_counts <- combined_df_chr %>%
#   count(Chromatin_State_HepG2) %>%
#   rename(Total = n)
# 
# # Step 2: Calculate percentage of transitions for each pair
# transition_matrix <- combined_df_chr %>%
#   count(Chromatin_State_HepG2, Chromatin_State_K562) %>%
#   left_join(total_counts, by = "Chromatin_State_HepG2") %>%
#   mutate(Percentage = n / Total) %>%
#   complete(Chromatin_State_HepG2 = as.factor(unified_levels), 
#            Chromatin_State_K562 = as.factor(unified_levels), 
#            fill = list(Percentage = 0))  # Fill missing transitions with 0%
# 
# # Step 3: Convert to long format
# transition_long <- transition_matrix %>%
#   mutate(
#     Chromatin_State_HepG2 = factor(Chromatin_State_HepG2, levels = as.character(unified_levels)),
#     Chromatin_State_K562 = factor(Chromatin_State_K562, levels = as.character(unified_levels))
#   )
# 
# # Step 4: Plot the heatmap
# p <- ggplot(transition_long, aes(x = Chromatin_State_HepG2, y = Chromatin_State_K562, fill = Percentage)) +
#   geom_tile() +
#   scale_fill_gradient(low = "white", high = "blue", labels = scales::percent_format(accuracy = 1)) +  # Display as percentages
#   labs(
#     title = "Heatmap of Percentage-Based Chromatin State Transitions Between HepG2 and K562",
#     x = "HepG2 Chromatin State",
#     y = "K562 Chromatin State",
#     fill = "Percentage"
#   ) +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     axis.title.x = element_text(size = 12),
#     axis.title.y = element_text(size = 12)
#   )
# 
# ggsave(
#   filename = file.path(output_dir, "heatmap.png"),
#   plot = p,
#   width = 12,   # Set the desired width
#   height = 8,   # Set the desired height
#   dpi = 300     # Set the resolution for the image
# )
# ################################################################################
violin_plot <- ggplot(df_filtered, aes(x = Biosample, y = fRead, fill = Biosample)) +
  geom_violin(trim = TRUE) +
  scale_fill_manual(values = group.colors) +  # Use the group.colors for fill
  labs(
    title = "Distribution of fRead by Biosample",
    x = "Biosample",
    y = "fRead"
  ) +
  theme_minimal()

# Display the violin plot
#print(violin_plot)lui

# Save the violin plot to disk
ggsave(
  filename = file.path(output_dir, "fRead_violin_plot.png"),
  plot = violin_plot,
  width = 10, height = 6, units = "in"
)

rm(violin_plot)
gc()

box_plot <- ggplot(df_long_Chromatin_State_Biosample_fRead, aes(x = Biosample, y = fRead, fill = Biosample)) +
  scale_fill_manual(values = group.colors) +  # Use the group.colors for fill
  geom_boxplot(outliers = FALSE) +
  labs(
    title = "Distribution of fRead by Biosample",
    x = "Biosample",
    y = "fRead"
  ) +
  theme_minimal()

# Display the violin plot
#print(violin_plot)

# Save the violin plot to disk
ggsave(
  filename = file.path(output_dir, "box_plot2a.png"),
  plot = box_plot,
  width = 10, height = 6, units = "in"
)

rm(box_plot)
gc()


box_plot_states <- ggplot(df_filtered, aes(x = Chromatin_State, y = fRead, fill = Biosample)) +
  scale_fill_manual(values = group.colors) +  # Use the group.colors for fill
  geom_boxplot(outliers = FALSE) +
  labs(
    title = "Distribution of fRead by Biosample",
    x = "ChromHMM  states",
    y = "genome wide methylation values"
  ) +
  theme_minimal()

# Save the violin plot to disk
ggsave(
  filename = file.path(output_dir, "box_plot2bChromatin_State.png"),
  plot = box_plot_states,
  width = 10, height = 6, units = "in"
)

rm(box_plot_states)
gc()

# Facet by biosample
dfv <- df_filtered %>%
  dplyr::mutate(Chromatin_State = factor(Chromatin_State, levels = as.character(1:18)))

state_cols <- setNames(chromatin_state_colors, as.character(1:18))

p_violin_facets <- ggplot(dfv, aes(x = Chromatin_State, y = fRead, fill = Chromatin_State)) +
  geom_violin(trim = TRUE, scale = "width", color = NA) +
  stat_summary(fun = median, geom = "point", size = 0.8, color = "black") +
  facet_wrap(~ Biosample, ncol = 2) +
  scale_fill_manual(values = state_cols, guide = "none") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "CpG methylation by ChromHMM state",
       x = "ChromHMM state",
       y = "CpG methylation (fRead)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.spacing = unit(8, "pt"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 8),
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "violin_states_faceted.png"),
       p_violin_facets, width = 12, height = 8, dpi = 300, limitsize = FALSE)
rm(p_violin_facets)
gc()



# helper for IQR bars
iqr_fn <- function(y) data.frame(y = median(y, na.rm = TRUE),
                                 ymin = quantile(y, 0.25, na.rm = TRUE),
                                 ymax = quantile(y, 0.75, na.rm = TRUE))

p_violin_dodged <- ggplot(dfv, aes(x = Chromatin_State, y = fRead, fill = Biosample)) +
  geom_violin(trim = TRUE, scale = "width",
              position = position_dodge(width = 0.85),
              width = 0.85, color = NA, alpha = 0.85) +
  stat_summary(fun.data = iqr_fn, geom = "errorbar",
               position = position_dodge(width = 0.85),
               width = 0.18, linewidth = 0.35, color = "black") +
  stat_summary(fun = median, geom = "point",
               position = position_dodge(width = 0.85),
               size = 0.8, color = "black") +
  scale_fill_manual(values = group.colors) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Distribution of CpG methylation by ChromHMM state",
       x = "ChromHMM state", y = "CpG methylation (fRead)", fill = "Biosample") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 8),
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(output_dir, "violin_states_dodged.png"),
       p_violin_dodged, width = 12, height = 7, dpi = 300, limitsize = FALSE)
rm(p_violin_dodged)
gc()

# Create the bar plot
bar_plot_state_count <- ggplot(df_long_Chromatin_State_Biosample_fRead, aes(x = Chromatin_State, fill = Biosample)) +
  geom_bar(position = "dodge") +  # Use 'dodge' to place bars side by side
  scale_fill_manual(values = group.colors) +  # Use the custom colors for biosamples
  labs(
    x = "Chromatin State",
    y = "Count of Entries"
  ) +
  theme_minimal()

# Save the violin plot to disk
ggsave(
  filename = file.path(output_dir, "box_plot2dChromatin_State_count.png"),
  plot = bar_plot_state_count ,
  width = 10, height = 6, units = "in"
)
rm(bar_plot_state_count)
gc()


# nicer boxplot (dodged, slim boxes, median dots, tidy theme)
box_plot_states <- ggplot(
  df_filtered,
  aes(x = factor(Chromatin_State, levels = as.character(1:18)),
      y = fRead,
      fill = Biosample)
) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = 0.85),
    width = 0.70,
    linewidth = 0.3
  ) +
  # small black median dots to guide the eye
  stat_summary(
    fun = median, geom = "point",
    position = position_dodge(width = 0.85),
    size = 0.8, color = "black"
  ) +
  scale_fill_manual(values = group.colors) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "CpG methylation by ChromHMM state",
    x = "ChromHMM state",
    y = "CpG methylation"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(6, "pt"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(size = 8),
    plot.title = element_text(face = "bold")
  )


# Save the violin plot to disk
ggsave(
  filename = file.path(output_dir, paste0("box_plot_states.", picuture_file_extension)),
  plot = box_plot_states,
  width = 12, height = 7, units = "in",
  dpi = 300, limitsize = FALSE
)

rm(box_plot_states)
gc()


#################################################################


################## summer table ##############################################

# Calculate the summary table
genome_summary_table <- df_long_Chromatin_State_Biosample_fRead %>%
  # Calculate total entries per Biosample
  group_by(Biosample) %>%
  mutate(total_entries = n()) %>%
  # Group by Biosample and Chromatin_State
  group_by(Biosample, Chromatin_State, total_entries) %>%
  summarise(
    mean_fRead = mean(fRead, na.rm = TRUE),
    count_entries = n()
  ) %>%
  # Calculate percentage of entries
  mutate(percentage_entries = (count_entries / total_entries) * 100) %>%
  ungroup()

#sum(summary_table$percentage_entries[summary_table$Biosample == "HepG2"])
saveRDS(object = genome_summary_table,file = file.path(output_dir,"genome_summary_table.rds"))

################################################################################

# # Reshape the data: gather chromatin state columns into long format
# df_long <- combined_df_chr %>%
#   pivot_longer(
#     cols = starts_with("Chromatin_State"), 
#     names_to = "Biosample", 
#     values_to = "Chromatin_State"
#   ) %>%
#   mutate(
#     # Remove the "Chromatin_State_" prefix from the biosample names
#     Biosample = gsub("Chromatin_State_", "", Biosample)
#   )
# 
# # Plot the distribution of chromatin states, with dodge by biosample
# chromatin_plot <- ggplot(df_long, aes(x = Chromatin_State, fill = Biosample)) +
#   geom_bar(position = "dodge") +
#   scale_fill_manual(values = group.colors) +
#   labs(
#     title = "Distribution of Chromatin States by Biosample",
#     x = "Chromatin State",
#     y = "Count",
#     fill = "Biosample"
#   ) +
#   theme_minimal()
# 
# # Save the chromatin state distribution plot
# ggsave(
#   filename = file.path(output_dir, "chromatin_state_distribution.png"),
#   plot = chromatin_plot,
#   width = 10, height = 6, units = "in"
# )
# 
# rm(chromatin_plot, df_long)
# gc()  # Run garbage collection to free up memory
# 
# # Reshape the data: gather fRead columns into long format
# df_long_fRead <- combined_df_chr %>%
#   pivot_longer(
#     cols = starts_with("fRead"),  # Gather all fRead columns
#     names_to = "Biosample",        # Create a new column called "Biosample"
#     values_to = "fRead"            # Store the values in a column called "fRead"
#   ) %>%
#   mutate(
#     # Remove the "fRead_" prefix from the biosample names
#     Biosample = gsub("fRead_", "", Biosample),
#     
#     # Convert Chromatin_State columns (for each biosample) to factors
#     across(starts_with("Chromatin_State"), as.factor)
#   )
# 
# # Assuming you want to use one chromatin state column for plotting (e.g., Chromatin_State_A549)
# fread_boxplot <- ggplot(df_long_fRead, aes(x = Chromatin_State_A549, y = fRead, fill = Biosample)) +
#   geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +
#   scale_fill_manual(values = group.colors) +  # Use predefined colors for biosamples
#   labs(
#     title = "Distribution of fRead by Chromatin State and Biosample",
#     x = "Chromatin State",
#     y = "fRead",
#     fill = "Biosample"
#   ) +
#   theme_minimal()
# 
# # Save the fRead boxplot
# ggsave(
#   filename = file.path(output_dir, "fRead_boxplot.png"),
#   plot = fread_boxplot,
#   width = 10, height = 6, units = "in"
# )
# 
# rm(fread_boxplot, df_long_fRead)
# gc()  # Run garbage collection to free up memory