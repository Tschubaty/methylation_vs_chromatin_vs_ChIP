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
library(dplyr)
library(ggplot2)
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

##################################### Functions #################################

# Load data
combined_df_chr <- read_delim(file = file.path(input_WGBS_dir, "combined_methylation_data.bed"))

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

head(combined_df_chr)

#######################################################################################

# Load necessary libraries
library(readr)
library(dplyr)
library(betareg)
library(ggplot2)

protein <- "CTCF"

# File path
file_path <- list.files(path = file.path(input_ChIP_dir,protein),full.names = TRUE)[1]


# Load the ChIP data
df_ChIP <- read_bed_file(file_path)

# Add a column to indicate peak presence in the whole genome dataset
combined_df_chr <- combined_df_chr %>%
  mutate(Peak_Present = ifelse(paste(chr, Start, End) %in% paste(df_ChIP$chr, df_ChIP$start_cg, df_ChIP$end_cg), 1, 0))

# Define a small epsilon
epsilon <- 1e-6

# Apply epsilon adjustment to transformed values
combined_df_chr <- combined_df_chr %>%
  mutate(
    # Scale from percentages (0-100) to proportions (0-1)
    fRead_A549_scaled = fRead_GM12878 / 100,
    
    # Add epsilon adjustment to ensure values are strictly within (0, 1)
    fRead_A549_transformed = pmax(epsilon, pmin(1 - epsilon, fRead_A549_scaled))
  )

# Verify the range of the transformed values
range(combined_df_chr$fRead_A549_transformed, na.rm = TRUE)


# Initialize an empty list to store models and results
beta_models <- list()
results <- data.frame()

for (state in names(chromatin_state_colors)) {
  print(state)
  # debug state <-  names(chromatin_state_colors)[1]
  # Filter the data for the current chromatin state
  state_data <- combined_df_chr %>%
    filter(Chromatin_State_A549 == state)
  
  
  state_data <- state_data %>%
    mutate(fRead_A549_scaled = fRead_A549 / 100)
  
  # Filter the data for the current chromatin state
  state_data <- combined_df_chr %>%
    filter(
      Chromatin_State_A549 == state,                       # Include only the current chromatin state
      !is.na(fRead_A549_transformed),                     # Remove NA values in the transformed methylation column
      fRead_A549_transformed > 0,                         # Ensure values are greater than 0
      fRead_A549_transformed < 1                          # Ensure values are less than 1
    )
  
  # Fit the beta regression model
  model <- betareg(
    fRead_A549_transformed ~ Peak_Present,
    data = state_data
  )
  

  
  # Save the model summary
  beta_models[[as.character(state)]] <- summary(model)
  
  # Extract coefficients and p-values for results
  coefficients <- as.data.frame(coef(summary(model)))
  coefficients$Chromatin_State <- state
  coefficients$Variable <- rownames(coefficients)
  
  # Append to results
  results <- rbind(results, coefficients)
}

# Save all model summaries to a file
capture.output(
  lapply(names(beta_models), function(state) {
    cat(paste("\nModel Summary for Chromatin State:", state, "\n"))
    print(beta_models[[state]])
  }),
  file = "separate_beta_model_summaries.txt"
)

# Save results as a table
write_delim(results, "separate_beta_model_results.tsv", delim = "\t")

# Visualization for each chromatin state
ggplot(results, aes(x = Chromatin_State, y = Estimate, fill = Variable)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(
    title = "Beta Regression Coefficients by Chromatin State",
    x = "Chromatin State",
    y = "Coefficient Estimate"
  ) +
  theme_minimal()
