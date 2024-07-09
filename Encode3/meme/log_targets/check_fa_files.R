# Documentation -----------------------------------------------------------
##
##
##  input: 
##
##
##  output: 
##
##
##  v_01 07.07.2024
##  Author: Daniel Batyrev (HUJI 777634015)
##
# Set up Work Environment --------------------------------------------------
#Clear R working environment
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

# Load necessary library
library(dplyr)

# Set file paths
chip_targets_file <- "chip_targets.txt"
idmapping_file <- "idmapping_reviewed_true_AND_model_organ_2024_06_16.tsv"

# Read the chip targets file
chip_targets <- readLines(chip_targets_file)

# Read the idmapping file
idmapping_df <- read.delim(idmapping_file, header = TRUE, sep = "\t")

# Extract the 'HGNC symbol' column (assuming the column name is 'HGNC.symbol')
targets_in_table <- unique(idmapping_df$From)

# Find targets in idmapping table that are not in chip targets
missing_in_chip_targets <- setdiff(targets_in_table, chip_targets)

# Extract the 'From' column
targets_in_table <- idmapping_df$From

# Find duplicates in the 'From' column
multiple_id <- duplicated(idmapping_df$From) | duplicated(idmapping_df$From, fromLast = TRUE)
#idmapping_df[multiple_id, 1:4]

non_trivial_match <- paste0(idmapping_df$From, "_HUMAN") != idmapping_df$Entry.Name
#idmapping_df[non_trivial_match,1:4 ]

idmapping_df[multiple_id & non_trivial_match, 1:4]

multiple_annotation <- idmapping_df[multiple_id, 1:4]
multiple_annotation[order(multiple_annotation$From),1:4]


idmapping_df_reviewed <- idmapping_df[!((duplicated(idmapping_df$From) | duplicated(idmapping_df$From, fromLast = TRUE)) & (paste0(idmapping_df$From, "_HUMAN") != idmapping_df$Entry.Name)), ]


idmapping_df_non_reviewed <- idmapping_df[!(!((duplicated(idmapping_df$From) | duplicated(idmapping_df$From, fromLast = TRUE)) & (paste0(idmapping_df$From, "_HUMAN") != idmapping_df$Entry.Name))),1:4 ]

idmapping_df_non_reviewed[! idmapping_df_non_reviewed$From %in% idmapping_df_reviewed$From,]