# Documentation -----------------------------------------------------------
##
##
##  input: 
##
##
##  output: 
##
##
##  v_01 21.07.2024
##  Author: Daniel Batyrev (HUJI 777634015)
##
# Set up Work Environment --------------------------------------------------
# Clear R working environment
rm(list = ls())
cluster <- FALSE
if (cluster) {
  this.dir <- "/ems/elsc-labs/meshorer-e/daniel.batyrev/HumanEvo/HumanEvo/"
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
#################################################################################

# Load necessary libraries
library(dplyr)
library(purrr)
library(readr)
library(stringr)

# Base directory path
base_dir <- "D:/Users/Daniel Batyrev/Documents/GitHub/meme/single_zip"

# Function to find all fimo.tsv files in directories containing the specified protein
find_fimo_files <- function(base_dir, protein) {
  all_dirs <- list.dirs(base_dir, recursive = TRUE)
  # Corrected grep pattern with underscore
  protein_dirs <- all_dirs[grepl(paste0("_", protein, "-"), all_dirs)]
  fimo_files <- file.path(protein_dirs, "fimo_out_1", "fimo.tsv")
  existing_fimo_files <- fimo_files[file.exists(fimo_files)]
  return(existing_fimo_files)
}

# Function to extract folder name
extract_folder_name <- function(file_path) {
  folder_name <- basename(dirname(dirname(file_path)))
  return(folder_name)
}

# Function to load a fimo.tsv file into a dataframe with preserved column names
load_fimo_file <- function(file_path) {
  df <- read_tsv(file_path, col_types = cols())
  return(df)
}

# Example protein name (you can change this to any protein)
protein_name <- "MAX"

# Find all relevant fimo.tsv files for the specified protein
fimo_files <- find_fimo_files(base_dir, protein_name)

# Load all fimo.tsv files into a named list of dataframes
fimo_data_list <- map(fimo_files, load_fimo_file) %>% 
  set_names(map_chr(fimo_files, extract_folder_name))

# Check the named list
print(names(fimo_data_list))  # Should print folder names corresponding to the protein name

# Print the first motif_id from each dataframe in fimo_data_list
for (name in names(fimo_data_list)) {
  cat("Folder:", name, "\n")
  cat("First motif_id:", fimo_data_list[[name]]$motif_id[1], "\n\n")
}
