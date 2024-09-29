# Load necessary libraries
library(dplyr)
library(stringr)
library(readr)
library(rlang) 

# Define the path to the directory where the files are stored
path <- "C:/Users/Daniel Batyrev/Documents/GitHub/methylation_vs_chromatin_vs_ChIP/Encode3/fimo_methylation/"

# Function to extract and clean column names from the first row of a text file
extract_colnames <- function(file_path) {
  # Read the first row
  first_row <- readLines(file_path, n = 1)
  
  # Remove the leading "#" and replace double tabs with single tab
  cleaned_row <- gsub("#", "", first_row)          # Remove leading #
  cleaned_row <- gsub("\t\t", "\t", cleaned_row)   # Replace double tabs with single tab
  
  # Split the cleaned string into individual column names
  colnames <- unlist(strsplit(cleaned_row, "\t"))
  
  return(colnames)
}

# Function to read .bed files with dynamically extracted column names and correct column types
read_bed_file <- function(file_path) {
  # Extract the correct column names from the first row
  colnames <- extract_colnames(file_path)
  
  # # Define column types: integer for all except the specified character columns
  # col_types <- cols(
  #   Chromosome = col_character(),  # Chromosome as character
  #   Sample_Info = col_character(), # Sample info as character
  #   Strand = col_character(),      # Strand as character (e.g., "+" or "-")
  #   CpG_Island_Status = col_character(), # CpG Island status as character
  #   
  #   # Assume all other columns are integers, dynamically handle columns from the colnames
  #   .default = col_integer()       # All other columns are integers
  # )
  
  # Read the BED file, skipping the first line (header row)
  df <- readr::read_delim(file_path, delim = "\t", skip = 1, col_names = colnames)# , col_types = col_types
  
  return(df)
}

# Function to process each file and extract relevant data
process_file <- function(file) {
  # Extract the protein and biosample from the file name
  filename <- basename(file)
  protein_biosample <- str_extract(filename, "^[^_]+_[^_]+")
  protein <- str_split(protein_biosample, "_")[[1]][1]
  biosample <- str_split(protein_biosample, "_")[[1]][2]
  
  # Read the data from the file
  data <- read_bed_file(file_path = file)
  
  # Group the data by Chromatin State (8th column) and count the occurrences (CpG rows)
  chromatin_state_counts <- data %>%
    group_by(ChromatinState = Chromatin_State) %>%
    summarise(number_CpG = n()) %>%
    mutate(Protein = protein, Biosample = biosample) %>%
    select(Protein, Biosample, ChromatinState, number_CpG)
  
  return(chromatin_state_counts)
}

# Get the list of all files matching the pattern *_methylation.chromatinstate.fimo.bed
files <- list.files(path, recursive = TRUE, full.names = TRUE, pattern = "_methylation\\.chromatinstate\\.fimo\\.bed$")

print(head(process_file(files[1])))

# Process all files and combine the results into a single data frame
results <- bind_rows(lapply(files, process_file))

# Print the result to the console
print(results)

# Optionally, save the result to a CSV file
write.csv(results, file = "chromatin_state_summary.csv", row.names = FALSE)
